##############################################################
###  Script to run GLM and Bayesian multiple regression    ###
###  Author:  D.K. Okamoto                                 ###
##############################################################

package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan",
                "dplyr","glmmTMB","parallel")

lapply(package.list,library,character.only=T)

### load urchin data ###
set.sum <- read.csv("1_Data/settlement_data_prepped.csv")

### global covariates 
covariates <- read.csv("1_Data/bayes_covariates.csv")[,-1]%>%
  mutate(region= factor(ifelse(SITE=="FBPC","NC",ifelse(SITE%in%c("OCNBCH","SIO"),"SD","SB"))))%>%
  filter(species=="SP")%>%
  dplyr::select(SITE,region,biweek,biweek_year,MEI_runlag2,NPGOlag,PDOlag)%>%
  data.frame()

set.agg <-  set.sum%>%
  group_by(SITE,biweek,biweek_year)%>%
  dplyr::summarize(P1= sum(SP)+0.5,
                   P2= sum(SF)+0.5,
                   TOT= sum(TOT))%>%
  filter((P1+P2-1)<TOT&P1>0.5)%>%
  data.frame()%>%
  mutate(ID = 1:nrow(.))


set.subset <- set.agg%>%
  select(ID,SITE,biweek,biweek_year)%>%
  right_join(set.sum)%>%
  mutate(ID= ifelse(is.na(ID),0,ID))

D_site_star <-as.matrix(dist(unique(covariates$biweek_year),upper= T,diag= T))
yearly <- seq(10,max(unique(covariates$biweek_year)),by = 26)
m= length(yearly)
n= length(unique(unique(covariates$biweek_year)))
D_star2 <-as.matrix(dist(yearly,upper= T,diag= T))
D_site_star2 <- as.matrix(dist(rbind(data.frame(x=unique(covariates$biweek_year)),
                                     data.frame(x=yearly))))[1:n, (n + 1):(n +m)]

D_site_star <-as.matrix(dist(1:702,upper= T,diag= T))
yearly <- seq(10,702,by = 26)
m= length(yearly)
n= length(unique(1:702))
NT= 702
phi  =1.5
phi2 = 1

NYM= 26
Cstar = array(NA,dim= c(NYM,NYM))
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
for (k in 1:(NT)) Cstar2[k, k] = 1;

w <- rnorm(NYM,0,1)
w2 <- rnorm(NT,0,1)
plot(c(t(chol(Cstar))%*%w,t(chol(Cstar))%*%w),type= "l")
plot(c(t(chol(Cstar2))%*%w2),type= "l")

sort_array <- matrix(1:nrow(covariates))
dim(sort_array) <- c(nrow(covariates)/7,7)

min_biweek <- min(covariates$biweek_year)
max_biweek <- max(covariates$biweek_year)

set.subset <- set.subset%>%
  data.frame()%>%
  filter(biweek_year>=min_biweek&biweek_year<=max_biweek)%>%
  mutate(biweek_year=biweek_year-min_biweek+1)

### bayesian model matrix
ModMat <-model.matrix(~region:MEI_runlag2+
                        region:PDOlag+
                        region:NPGOlag-1, data= covariates)

ModMat <- ModMat[,apply(ModMat,2,sd)!=0]
### list of data for model
Data <- with(set.subset,list(
  NO= nrow(set.subset),
  YSP = SP,
  NS = length(unique(SITE)),
  P1= set.agg$P1,
  P2= set.agg$P2,
  NM = dim(sort_array)[1],
  OBS = biweek_year+(as.numeric(SITE)-1)*max(biweek_year),
  NSUB = length((1:nrow(set.subset))[ID>0]),
  NSP= nrow(set.agg),
  ID= ID,
  SUB = (1:nrow(set.subset))[ID>0],
  NYM = 26,
  NP= ncol(ModMat),
  OBS_MONTH= biweek_year,
  PRED_MONTH= covariates$biweek[1:dim(sort_array)[1]],
  OBS_SITE= as.numeric(SITE),
  X = ModMat,
  xind = sort_array,
  D_seas = t(chol(Cstar)),
  D= Duration,
  N= ifelse(ID>0,TOT,SP),
  sigma_scale= 0.25
))
seed <- 12345

mod <- stan_model(file= "2_Analysis_Code/Horseshoe_mod_AR_2.stan") 

## run the Horsehshoe Regression
system.time(p_global <- sampling(mod,
  data=Data,
  iter =1000, warmup=500,
  chains =3,cores = 3,
  pars= c("beta","sigma_S_mu","phi")))
