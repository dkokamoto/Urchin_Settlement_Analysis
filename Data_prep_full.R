######################################################
###  Script to format data for analysis and plots  ###
###  Author:  D.K. Okamoto                         ###
######################################################

### load necessary packages 
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan",
                "dplyr")

lapply(package.list,library,character.only = T)

### load urchin data ###
settlement <- read.csv("1_Data/Invertebrate_Settlement_All_Years.csv",header=T)

### format dates, etc.

settlement <-  mutate(settlement, 
                      DATE_DEPLOYED=parse_date_time(DATE_DEPLOYED, "%m/%d/%y"),
                      DATE_RETRIEVED=parse_date_time(DATE_RETRIEVED, "%m/%d/%y"))


settlement <- mutate(settlement,
                     julian= as.numeric(julian(DATE_RETRIEVED)),
                     ID = paste(DATE_DEPLOYED,DATE_RETRIEVED,SITE,sep= "_"),
                     ID2 =paste(DATE_DEPLOYED,SITE,sep= "_"),
                     ID3 =paste(DATE_RETRIEVED,SITE,sep= "_"),
                     month_dep= month(DATE_DEPLOYED),
                     month_ret = month(DATE_RETRIEVED),
                     year_dep = year(DATE_DEPLOYED),
                     year_ret = year(DATE_RETRIEVED),
                     day_dep = day(DATE_DEPLOYED),
                     day_ret = day(DATE_RETRIEVED))


settlement <- settlement%>%
  mutate(diff = as.numeric(Duration),
         S_PURPURATUS = ifelse(S_PURPURATUS==-99999,NA,S_PURPURATUS),
         M_FRANCISCANUS = ifelse(M_FRANCISCANUS==-99999,NA, M_FRANCISCANUS),
         TOTAL_URCHINS = ifelse(TOTAL_URCHINS==-99999,NA, TOTAL_URCHINS))%>%
  filter(!is.na(M_FRANCISCANUS)&
           !is.na(S_PURPURATUS)&
           !is.na(TOTAL_URCHINS)&
           Duration>0) %>%
  arrange(SITE,julian) %>%
  drop.levels()


### format data to estimate collection period means per brush per day

set.sum <- settlement %>% group_by(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration) %>%
  mutate(SP= S_PURPURATUS,
         SF= M_FRANCISCANUS,
         TOT= TOTAL_URCHINS) %>%
  select(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration,SP,SF,TOT)%>%
  mutate(SP_EM= as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SP/(SP+SF)*TOT)/Duration)),
         SF_EM = as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration)),
         julian = julian(DATE_RETRIEVED, origin = "1990-01-01"),
         SITE_NUM = as.numeric(SITE),
         yday = yday(DATE_RETRIEVED),
         biweek  = floor(yday/14)+1) %>%
  mutate(biweek =ifelse(biweek>26,26,biweek),
         biweek_year  = biweek+(year_ret-1990)*26) %>%
  filter(!is.na(TOT)&Duration<96) 

site_levels <- c("Anacapa[SB]","Fort Bragg[NorCal]","Gaviota[SB]","Ocean Beach[SD]","Ellwood[SB]","Stearns Wharf[SB]","Scripps[SD]")

set.monyr <- set.sum %>%
  group_by(biweek_year,SITE)%>% 
  summarise(mean_SP= mean( SP_EM),mean_SF= mean( SF_EM)) %>%
  mutate(site= factor(SITE,levels=levels(SITE),labels= site_levels))

site <- 1:nlevels(factor(set.monyr$SITE))
data_mat <- expand.grid(site=site_levels,biweek=1:26, YEAR= 1990:2017) %>%
  mutate(site=factor(site,levels= site_levels),
         biweek_year  = biweek+(YEAR-1990)*26,
         M1 = cos(2*pi*biweek/26),
         M2 = sin(2*pi*biweek/26),
         M3 = cos(2*pi*biweek/26/2),
         M4 = sin(2*pi*biweek/26/2),
         M3 = cos(2*pi*biweek/26/4),
         M4 = sin(2*pi*biweek/26/4),
         MID= 1:length(site))

data_mat2 <- expand.grid(site=site_levels,biweek=1:26) %>%
  mutate(site=factor(site,levels= site_levels),
         biweek  = biweek,
         M1 = cos(2*pi*biweek/26),
         M2 = sin(2*pi*biweek/26),
         M3 = cos(2*pi*biweek/26/2),
         M4 = sin(2*pi*biweek/26/2),
         M3 = cos(2*pi*biweek/26/4),
         M4 = sin(2*pi*biweek/26/4),
         MID= 1:length(site))

D_site_star <-as.matrix(dist(unique(data_mat$biweek_year),upper= T,diag= T))
yearly <- seq(10,max(data_mat$biweek_year),by = 26)
m= length(yearly)
n= length(unique(data_mat$biweek_year))
D_star2 <-as.matrix(dist(yearly,upper= T,diag= T))
D_site_star2 <- as.matrix(dist(rbind(data.frame(x=unique(data_mat$biweek_year)),
                                     data.frame(x=yearly))))[1:n, (n + 1):(n +m)]

set.sum2 <- join(data_mat,set.monyr) 

### algebraic observed settlement 
obs_SF <- matrix(set.sum2$mean_SF, ncol= 7, byrow= T)
obs_SP <- matrix(set.sum2$mean_SP, ncol= 7, byrow= T)

### define GP covariance matrices

phi  =1.5
phi2 = 1

NYM= 26
Cstar = array(NA,dim= c(NYM,NYM))

NT= 728
Cstar2 = array(NA,dim= c(NT,NT))

for (i in 1:(NYM-1)) {
  for (j in (i + 1):(NYM)) {
    Cstar[i, j] =  exp(-2*(0.5-0.5*cos(2*pi*abs(D_site_star[i, j]/26)))/((phi)^2));
    Cstar[j, i] = Cstar[i, j];
  }
}
for (k in 1:(NYM)) Cstar[k, k] = 1+0.0001;


w <- rnorm(NYM,0,1)
plot(c(t(chol(Cstar))%*%w,t(chol(Cstar))%*%w),type= "l")


for (i in 1:(NT-1)) {
  for (j in (i + 1):(NT)) {
    Cstar2[i, j] =   (1+(3^(0.5)*abs(D_site_star[i,j]/26)/(phi2)))*exp(-(3^(0.5)*abs(D_site_star[i, j]/26)/(phi2)));
    Cstar2[j, i] = Cstar2[i, j];
  }
}
for (k in 1:(NT)) Cstar2[k, k] = 1+0.0001;

z <- rnorm(NT,0,1)
plot(t(chol(Cstar2))%*%z,type= "l")
abline(v= seq(0,702,by= 26))

### list of data for model 
data <- with(set.sum,list(
  NO= nrow(set.sum),
  YSP = SP,
  NB = rep(1,nrow(set.sum)),
  NS = length(unique(SITE)),
  NM = length(unique(data_mat$biweek+(data_mat$YEAR)*26)),
  NT= m,
  NYM = 26,
  OBS_MONTH= biweek_year,
  PRED_MONTH= rep(unique(data_mat$biweek),length(unique(data_mat$YEAR))),
  OBS_SITE= as.numeric(SITE),
  D_seas = t(chol(Cstar)),
  D_ann = t(chol(Cstar2)),
  D= Duration,
  n= SP+SF,
  N= TOT
))