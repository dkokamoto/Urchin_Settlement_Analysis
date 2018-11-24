######################################################
###  Script to run Bayesian model on settlement    ###
###  Author:  D.K. Okamoto                         ###
######################################################

rm(list=ls())
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan",
                "dplyr")

### install packages if you don't have them 
### not run
#lapply(package.list,install.packages)
lapply(package.list,library,character.only=T)

### number of iterations and chains 
n.chains =4
n.iter =2000
n.burnin=1000
set.seed <- 1234
### load urchin data ###
set.sum <- read.csv("1_Data/settlement_data_prepped.csv")

### create the site level data matrix
data_mat <- expand.grid(site=levels(set.sum$SITE),biweek=1:26, YEAR= 1990:2016) %>%
  mutate(site=factor(site,levels= levels(set.sum$SITE)),
         biweek_year  = biweek+(YEAR-1990)*26)

### create the gaussian process covariance matrices (i.e. annual and seasonal scale kernel)
D_site_star <-as.matrix(dist(unique(data_mat$biweek_year),upper= T,diag= T))
yearly <- seq(10,max(unique(data_mat$biweek_year)),by = 26)
NY = length(yearly)
NM = length(unique(unique(data_mat$biweek_year)))
D_star2 <-as.matrix(dist(yearly,upper= T,diag= T))
D_site_star2 <- as.matrix(dist(rbind(data.frame(x=unique(data_mat$biweek_year)),
                                     data.frame(x=yearly))))[1:NY, (NY + 1):(NY +NM)]

### define GP covariance matrices
phi  =1.5
phi2 = 1.5

NYM= 26
Cstar = array(NA,dim= c(NYM,NYM))

NT= NM
Cstar2 = array(NA,dim= c(NT,NT))
Cstar3 = array(NA,dim= c(NT,NT))

### set up the seasonal-scale gaussian process covariance matrix
for (i in 1:(NYM-1)) {
  for (j in (i + 1):(NYM)) {
    Cstar[i, j] =  exp(-2*(0.5-0.5*cos(2*pi*abs(D_site_star[i, j]/26)))/((phi)^2));
    Cstar[j, i] = Cstar[i, j];
  }
}
for (k in 1:(NYM)) Cstar[k, k] = 1;

### plot a realization of the seasonal gaussian process
w <- rnorm(NYM,0,1)
plot(c(t(chol(Cstar))%*%w,t(chol(Cstar))%*%w),type= "l")
abline(v= seq(0,26*2,by= 26))

### set up the annual-scale gaussian process covariance matrix
for (i in 1:(NT-1)) {
  for (j in (i + 1):(NT)) {
    Cstar3[i, j] =   (1+(3^(0.5)*abs(D_site_star[i,j]/26)/(phi2)))*exp(-(3^(0.5)*abs(D_site_star[i, j]/26)/(phi2)));
    Cstar3[j, i] = Cstar3[i, j];
  }
}
for (k in 1:(NT)) Cstar3[k, k] = 1;

### plot a realization of the annual gaussian process
z <- rnorm(NT,0,1)
plot(t(chol(Cstar3))%*%z,type= "l")
lines(t(chol(Cstar3))%*%z,type= "l",col= "red")
abline(v= seq(0,NT,by= 26))

### list of data for model
data <- with(set.sum,list(
  NO= nrow(set.sum),
  P1= SP[(SP+SF)<TOT]+0.5,
  P2= SF[(SP+SF)<TOT]+0.5,
  NB = 1,
  NS = length(unique(SITE)),
  NSUB = length((1:nrow(set.sum))[(SP+SF)<TOT]),
  NM = length(unique(data_mat$biweek+(data_mat$YEAR)*26)),
  NYM = 26,
  OBS = biweek_year+(as.numeric(SITE)-1)*NT,
  PRED_MONTH= rep(unique(data_mat$biweek),length(unique(data_mat$YEAR))),
  D_seas = t(chol(Cstar)),
  D_ann = t(chol(Cstar3)),  ### cholesky factorization of 
  D= Duration,
  N= ifelse((SP+SF)<TOT,TOT,SP),
  SUB = (1:nrow(set.sum))[(SP+SF)<TOT]
))

### compile the model
mod <- stan_model(file= "3_Analysis_Code/Purp_Trends.stan")

### estimate posterior 
post_GP <- sampling(mod, 
                data=data,
                iter = 2000, warmup= 500,
                chains =3,cores =3,
                control= list(adapt_delta= 0.99,
                              max_treedepth= 15),
                pars= c("Omega_S",
                        "Omega_Spatial",
                        "LS",
                        "sigma",
                        "sigma2",
                        "S",
                        "mu_S",
                        "sigma_S_mu",
                        "w"),
                verbose = TRUE,init="random")