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
settlement <- read.csv("./1_Data/Invertebrate_Settlement_All_Years.csv",header=T)

### format dates, etc.
settlement <-  mutate(settlement, 
                       DATE_RETRIEVED=parse_date_time(DATE_RETRIEVED, "%Y-%m-%d"),
                       DATE_RETRIEVED=parse_date_time(DATE_RETRIEVED, "%Y-%m-%d"))

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

settlement2 <-  settlement %>% mutate(S_FRANCISCANUS = replace(S_FRANCISCANUS, 
      year_ret==2009&
      month_ret==04&
        day_ret==07&
        is.na(S_FRANCISCANUS)&
        !is.na(S_PURPURATUS)
        , 0))  

settlement <- settlement%>%
              mutate(diff = as.numeric(Duration)) %>%
              filter(!is.na(S_FRANCISCANUS)&
                       !is.na(S_PURPURATUS)&
                       !is.na(TOTAL_URCHINS)&
                       Duration>0) %>%
              arrange(SITE,julian) %>%
              drop.levels


### format data to estimate collection period means per brush per day

set.sum <- settlement %>% group_by(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration) %>%
  summarise(SP= sum(S_PURPURATUS),
            SF= sum(S_FRANCISCANUS),
            TOT= sum(TOTAL_URCHINS),
            NB= length(TOTAL_URCHINS)) %>%
  mutate(SP_EM= as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SP/(SP+SF)*TOT)/Duration/NB)),
         SF_EM = as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration/NB)),
         julian = julian(DATE_RETRIEVED, origin = "1990-01-01"),
         SITE_NUM = as.numeric(SITE),
         M1 = cos(2*pi*biweek/26),
         M2 = sin(2*pi*biweek/26),
         yday = yday(DATE_RETRIEVED),
         biweek  = floor(yday/14)+1) %>%
  mutate(biweek =ifelse(biweek>26,26,biweek),
          biweek_year  = biweek+(year_ret-1990)*26) %>%
  filter(!is.na(TOT)&Duration<45) 

site_levels <- c("Anacapa[SB]","Fort Bragg[NorCal]","Gaviota[SB]","Ocean Beach[SD]","Ellwood[SB]","Stearns Wharf[SB]","Scripps[SD]")

set.monyr <- set.sum %>%
              group_by(biweek_year,SITE)%>% 
                summarise(mean_SP= mean( SP_EM),mean_SF= mean( SF_EM)) %>%
              mutate(site= factor(SITE,levels=levels(SITE),labels= site_levels))

site <- 1:nlevels(factor(set.monyr$SITE))
data_mat <- expand.grid(site=site_levels,biweek=1:26, YEAR= 1990:2015) %>%
              mutate(site=factor(site,levels= site_levels),
                     biweek_year  = biweek+(YEAR-1990)*26,
                     M1 = cos(2*pi*biweek/26),
                     M2 = sin(2*pi*biweek/26),
                     M3 = sin(2*pi*biweek/13),
                     M4 = sin(2*pi*biweek/13),
                     MID= 1:length(site))

set.sum2 <- join(data_mat,set.monyr) 

### algebraic observed settlement 
obs_SF <- matrix(set.sum2$mean_SF, ncol= 7, byrow= T)
obs_SP <- matrix(set.sum2$mean_SP, ncol= 7, byrow= T)

### linear model matrix
MM <- model.matrix(MID~M1:factor(site)+M2:factor(site)+M3:factor(site)+M4:factor(site)+factor(site)-1,data= data_mat)
Mt <- t(MM)
dim(Mt) <- c(ncol(MM),7,nrow(MM)/7)
MA <- aperm(Mt,c(3,2,1))

### list of data for model 
data <- with(set.sum,list(
  NO= nrow(set.sum),
  YSP = SP,
  NB = NB,
  NS = length(unique(SITE)),
  NM = length(unique(data_mat$biweek+(data_mat$YEAR)*26)),
  OBS_MONTH= biweek_year,
  OBS_SITE= as.numeric(SITE),
  MM =MA,
  D= Duration,
  n= SP+SF,
  N= TOT,
  NP = dim(MM)[2]
))