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

NYM= 26
D_site_star <-as.matrix(dist(1:NYM,upper= T,diag= T))
phi2 = 1

Cstar = array(NA,dim= c(NYM,NYM))

for (i in 1:(NYM-1)) {
  for (j in (i + 1):(NYM)) {
    Cstar[i, j] =  exp(-2*(0.5-0.5*cos(2*pi*abs(D_site_star[i, j]/26)))/((phi)^2));
    Cstar[j, i] = Cstar[i, j];
  }
}
for (k in 1:(NYM)) Cstar[k, k] = 1+0.0001;


w <- rnorm(NYM,0,1)
plot(c(t(chol(Cstar))%*%w,t(chol(Cstar))%*%w),type= "l")

### local covariates 
covariates <- read.csv("1_Data/bayes_covariates.csv")[,-1]%>%
  mutate(region= factor(ifelse(SITE=="FBPC","NC",ifelse(SITE%in%c("OCNBCH","SIO"),"SD","SB"))))%>%
  filter(species=="SP")%>%
  select(SITE,region,biweek,biweek_year,adults,SST_rollmean_30,chla_rollmean_30,
         BK_stand_30,kelp_biomass)%>%
  group_by(SITE)%>%
  dplyr::mutate_at(c("SST_rollmean_30","chla_rollmean_30"),function(x)(ifelse(x==0,0,(x-mean(x,na.rm=T))/(sd(x,na.rm=T)))))%>%
  dplyr::mutate_at(c("adults","kelp_biomass"), function(x)ifelse(is.na(x)|x==0,0,(log(x)-mean(log(x)))/(sd(log(x),na.rm=T))))%>%
  data.frame()%>%
  mutate(region= factor(ifelse(SITE=="FBPC","NC",ifelse(SITE%in%c("OCNBCH","SIO"),"SD","SB"))))%>%
  filter(!is.na(chla_rollmean_30))%>%
  mutate(adults= ifelse(is.na(adults),0,adults),
         kelp_biomass= ifelse(is.na(kelp_biomass),0,kelp_biomass))

sort_array <- matrix(1:nrow(covariates))
dim(sort_array) <- c(nrow(covariates)/7,7)

min_biweek <- min(covariates$biweek_year)
max_biweek <- max(covariates$biweek_year)

sort_array <- matrix(1:nrow(covariates))
dim(sort_array) <- c(nrow(covariates)/7,7)

min_biweek <- min(covariates$biweek_year)
max_biweek <- max(covariates$biweek_year)

set.subset<- set.sum%>%
  data.frame()%>%
  filter(biweek_year>=min_biweek&biweek_year<=max_biweek)%>%
  mutate(biweek_year=biweek_year-min_biweek+1)

set.agg <-  set.subset%>%
  group_by(SITE,biweek,biweek_year)%>%
  dplyr::summarize(P1= sum(SP)+0.5,
                   P2= sum(SF)+0.5,
                   TOT= sum(TOT))%>%
  filter((P1+P2-1)<TOT)%>%
  data.frame()%>%
  mutate(ID = 1:nrow(.))

set.subset <- set.agg%>%
  select(ID,SITE,biweek,biweek_year)%>%
  right_join(set.subset)%>%
  mutate(ID= ifelse(is.na(ID),0,ID))

### bayesian model matrix
ModMat <-model.matrix(~ region:chla_rollmean_30 +
  region:SST_rollmean_30 +
  region:BK_stand_30 +
  region:kelp_biomass +
  region:adults-1, data= covariates)

ModMat <- ModMat[,apply(ModMat,2,sd)!=0]

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
  sigma_scale= 0.25,
  m0=2,
  slab_scale= 1
))

mod <- stan_model(file= "2_Analysis_Code/Horseshoe_mod.stan") 

p_local <- sampling(mod,
    data=Data,
    iter =2000, warmup=1000,
    control= list(adapt_delta= 0.95,
                  max_treedepth = 15),
    chains =3,cores = 3,
    pars= c("beta","phi","mu_S"))