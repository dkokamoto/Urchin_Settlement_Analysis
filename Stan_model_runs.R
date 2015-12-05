
setwd("~/Copy/UrchinAnalyses/R Code/Urchin_Settlement_Analysis/")

rm(list=ls())
package.list<-c("abind","AER","lubridate","bitops","car","chron","coda","colorspace","dichromat","digest","Formula",
                "gdata","ggplot2","gpclib","gtable","gtools","Hmisc","labeling","lmtest","lubridate",
                "memoise","munsell","mvtnorm","ncdf","plyr","proto","R.methodsS3","R.oo","R.utils","R2jags",
                "R2WinBUGS","RColorBrewer","RCurl","reshape2","rjags","sandwich","scales","sp",
                "stringr","timeDate","zoo","maps","maptools","classInt","rgeos","fields","rgdal")

### install packages if you don't have them 
### not run
#lapply(package.list,install.packages)
lapply(package.list,library,character.only=T)
### load data and the model 
source("Stan_model.R")
source("Data_prep.R")

### number of iterations and chains 
n.chains =3
n.iter =1000
n.burnin=500
set.seed <- 1234

### parameters to save
params <- c("Omega_SP","Omega_SF","SP","SF","LSP","LSF","beta_SP","beta_SF","phi_SP","phi_SF")

### linear model matrix
MM <- model.matrix(MID~factor(month_ret):factor(SITE)-1,data= data_mat)
Mt <- t(MM)
dim(Mt) <- c(ncol(MM),7,nrow(MM)/7)
MA <- aperm(Mt,c(3,2,1))


### list of data for model 
data <- with(set.sum,list(
  NO= nrow(set.sum),
  YSP = SP,
  NB = NB,
  NS = length(unique(SITE)),
  NM = length(unique(data_mat$month+(data_mat$YEAR)*12)),
  OBS_MONTH= monyr,
  OBS_SITE= as.numeric(SITE),
  MM =MA,
  D= Duration,
  n= SP+SF,
  N= TOT,
  NP = dim(MM)[2]
))

### compile the models
opt_modelAR<-stan_model(model_code=modelAR)

### optimize for Max Likelihood (might take a few attemps... )
Fit.test <- optimizing(opt_modelAR,data= data,iter=10000)

### estimate posterior 
postAR <- stan(model_code = modelAR, 
               data=data,
               pars=params,
               iter = n.iter , warmup= n.burnin,chains =n.chains,
               verbose = FALSE,init="random",seed= set.seed)
