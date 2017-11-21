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
### load data and the model 
source("2_Data_Formatting_Code/Data_prep_full.R")

### number of iterations and chains 
n.chains =2
n.iter =1000
n.burnin=500
set.seed <- 1234

# ### compile the models
model_GP3 <- stan_model(file= "3_Analysis_Code/GP_Matern7.stan")
# 
# ### optimize for Max Likelihood (might take a few attemps... )
Fit.test <- optimizing(model_GP3,data= data,iter=100)

plot(Fit.test$par[names(Fit.test$par)%in%paste0("w[1,6,",1:676,"]")],type= "l",col= "blue")

# ### estimate posterior 
 post_GP <- stan(file= "3_Analysis_Code/GP_Matern7.stan", 
                data=data,
                iter = 2000 , warmup= 1000,
                chains =3,cores = 3,
                verbose = TRUE,init="random")
 
pars <- extract(post_GP3)
post_GP3@.MISC <- emptyenv()

### mean values
mus <- adply(pars$mu_SP,2,mean)
names(mus) <- c("site","mu")

### smoothed variance
sigmas <- adply(pars$sigma2[,1,],2,mean)
names(sigmas) <- c("site","sigma")

### seasonal variance
sigmas2 <- adply(pars$sigma[,1,],2,mean)
names(sigmas2) <- c("site","sigma")

fixef <- dcast(melt(apply(pars$w,c(2,3,4),quantile,probs= c(0.05,.5,0.95)),
                       varnames = c("quant","species","site","biweek_year")),species+site+biweek_year~quant)


names(fixef)[4:6] = c("CIL","med","CIU")
fixef$site <- factor(fixef$site,labels= levels(set.sum$SITE))
fixef$site <- factor(fixef$site,levels= levels(fixef$site)[c(2,3,5,6,1,4,7)])
fixef$site2 <- factor(ifelse(fixef$site%in%c("SIO","OCNBCH"),"San Diego",
                                ifelse(fixef$site%in%c("FBPC"),"NorCal","Santa Barbara")))
levels(fixef$site) = c("Fort Bragg","Gaviota","Ellwood","Stearns Wharf",
                       "Anacapa","Ocean Beach","Scripps Pier")
fixef$site2=factor(fixef$site2,levels= c("NorCal","Santa Barbara","San Diego"))

ggplot(aes(biweek_year,med),data= fixef)+
  geom_ribbon(aes(ymin= CIL,ymax= CIU),fill= "grey90")+
  geom_path()+
  facet_grid(site~species)

