######################################################
###  Script to run Bayesian model on settlement    ###
###  Author:  D.K. Okamoto                         ###
######################################################

rm(list=ls())
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan")

### install packages if you don't have them 
### not run
#lapply(package.list,install.packages)
lapply(package.list,library,character.only=T)
### load data and the model 
source("3_Analysis_Code/Stan_model.R")
source("2_Data_Formatting_Code/Data_prep.R")

### number of iterations and chains 
n.chains =3
n.iter =1000
n.burnin=500
set.seed <- 1234

### parameters to save
params <- c("Omega_SP","Omega_SF","SP","SF","LSP","LSF","beta_SP","beta_SF","phi_SP","phi_SF")



# ### compile the models
 opt_modelAR<-stan_model(model_code=modelAR)
# 
# ### optimize for Max Likelihood (might take a few attemps... )
 Fit.test <- optimizing(opt_modelAR,data= data,iter=10000)
# 
# ### estimate posterior 
 postAR <- stan(model_code = modelAR, 
                data=data,
                pars=params,
                iter = n.iter , warmup= n.burnin,
                chains =n.chains,
                verbose = FALSE,init="random",
                seed= set.seed,cores= 2)
