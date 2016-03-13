######################################################
###  Script to aggregate data and model output     ###
###  Author:  D.K. Okamoto                         ###
######################################################


library(rstan);library(ggplot2);library(MASS)
library(gdata);library(plyr);library(dplyr)
rm(list=ls())
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan","dplyr")

lapply(package.list,library,character.only = T)
### now run "Stan_model_runs.R" ###
#source("3_Analysis_Code/Stan_model_runs.R")

### run data prep script
source("2_Data_Formatting_Code/Data_prep.R")

### load fitted model results
load("5_Model_Output/postAR.RData") 

### generate monthly mean settlement from the data
set.ag <- ddply(set.sum,.(biweek_year,SITE,month_ret,year_ret),summarise, mean_SP = mean(SP_EM, na.rm= T),mean_SF = mean(SF_EM, na.rm= T))

params <- extract(postAR)

### settlement partial autocorrelation function
phi_SP <-  apply(params$phi_SP, c(2), mean)
phi_SF <- apply(params$phi_SF, c(2), mean)

### estimated mean (log-scale) purps and frans from the Stan posterior
E_SP <- apply(exp(params$LSP), c(2,3), mean)

E_SF <- apply(exp(params$LSF), c(2,3), mean)

### get monthly mean (seasonal) estimates
beta_SP <- apply(exp(params$beta_SP),2, mean)
beta_SF <- apply(exp(params$beta_SF),2, mean)

beta_SP_U <- apply(exp(params$beta_SP), c(2), quantile, probs= 0.975)
beta_SF_U <- apply(exp(params$beta_SF), c(2), quantile, probs= 0.975)
beta_SP_L <- apply(exp(params$beta_SP), c(2), quantile, probs= 0.025)
beta_SF_L <- apply(exp(params$beta_SF), c(2), quantile, probs= 0.025)

names(beta_SP) <- dimnames(MM)[2][[1]]
names(beta_SF) <- dimnames(MM)[2][[1]]
seas_mod_SP <- data.frame(matrix(beta_SP,ncol= 7))
names(seas_mod_SP) <- levels(set.sum$SITE)
seas_mod_SF <- data.frame(matrix(beta_SF,ncol= 7))
names(seas_mod_SF) <- levels(set.sum$SITE)
seas_mod_SF$biweek2 <- c(3:26,1:2)
seas_mod_SP$biweek2 <- c(3:26,1:2)

### convert each into a seasonal correlation matrix
seas_cor_SP <- matrix(cor(seas_mod_SP[1:7]), ncol= data$NS,
                   dimnames= list(levels(set.sum$SITE),
                                  levels(set.sum$SITE)))
seas_cor_SF <- matrix(cor(seas_mod_SF[1:7]), ncol= data$NS,
                   dimnames= list(levels(set.sum$SITE),
                                  levels(set.sum$SITE)))
### show correlation matrix 
levelplot(seas_cor_SP)
levelplot(seas_cor_SF)

### get correlation matrix into something we can plot
SP_mod_df <- melt(seas_mod_SP, id.vars= "biweek2",variable.name= "SITE",value.name= "Exp_seas_SP")
SP_mod_df$SP_U <- beta_SP_U 
SP_mod_df$SP_L <- beta_SP_L
SF_mod_df <- melt(seas_mod_SF, id.vars= "biweek2",variable.name= "SITE",value.name= "Exp_seas_SF")
SF_mod_df$SF_U <- beta_SF_U 
SF_mod_df$SF_L <- beta_SF_L
seas_mod <- join(SP_mod_df,SF_mod_df)
levels(seas_mod$SITE) <- site_levels
seas_mod$SITE <- factor(seas_mod$SITE, levels= rev(levels(seas_mod$SITE)[c(7,4,1,6,5,3,2)]))


### now get mean estimates 
SP_mod <- apply(exp(params$SP), c(2,3), mean)
SF_mod <- apply(exp(params$SF), c(2,3), mean)
SP_mod_L <- apply(exp(params$SP), c(2,3),quantile, probs= c(0.025))
SF_mod_L <- apply(exp(params$SF), c(2,3), quantile, probs= c(0.025))
SP_mod_U <- apply(exp(params$SP), c(2,3),quantile, probs= c(0.975))
SF_mod_U <- apply(exp(params$SF), c(2,3), quantile, probs= c(0.975))

format.data <- function(x,name) {
  x[is.na(obs_SF)] <- NA
  x <- data.frame(x)
  names(x) <- levels(set.sum$SITE)
  x$biweek_year <- 1:nrow(x)
  x <- melt(x, id.vars= "biweek_year",value_name= name)
  names(x) <- c("biweek_year","SITE",name)
  return(x)
}

SP_mod <- format.data(SP_mod,"Est_SP")
SF_mod <- format.data(SF_mod,"Est_SF")
SP_mod_U <- format.data(SP_mod_U,"SP_U")
SP_mod_L <- format.data(SP_mod_L, "SP_L")
SF_mod_U <- format.data(SF_mod_U, "SF_U")
SF_mod_L <- format.data(SF_mod_L, "SF_L")

mod_df <- join_all(list(SP_mod,SF_mod,SP_mod_U,SP_mod_L,SF_mod_U,SF_mod_L), by= c("biweek_year","SITE"))
set.ag2 <- join(set.ag,mod_df, by = c("SITE","biweek_year"))
set.ag2$MONTH <- set.ag2$month_ret
set.ag2$YEAR <- set.ag2$year_ret
set.ag2$biweek <- set.ag2$biweek_year-(set.ag2$YEAR-1990)*26
levels(set.ag2$SITE) <- levels(seas_mod$SITE)[c(5,1,2,6,3,4,7)]
set.ag2$SITE <-factor(set.ag2$SITE, levels = c("Fort Bragg[NorCal]","Gaviota[SB]","Ellwood[SB]","Stearns Wharf[SB]" ,"Anacapa[SB]","Ocean Beach[SD]","Scripps[SD]" ))   



set_summary <- with(set.ag2, data.frame(expand.grid(SITE=levels(SITE),biweek= 1:26,YEAR= min(YEAR):max(YEAR))))
set_summary$biweek_year<- set_summary$biweek+(set_summary$YEAR-1990)*26
set_summary$yday<- round((set_summary$biweek-0.5)/26*365)
set_summary$date <- with(set_summary,parse_date_time(paste(yday,YEAR,sep= "-"),"%j-%y"))

set_summary <- join(set_summary, set.ag2)


set_summary$SF_U <- with(set_summary,ifelse(SF_U>0.5,mean_SF*1.2, SF_U))

### generate partial autocorrelation matrix
pacf_seas_SP <- data.frame(list(phi= paste0("phi==",round(phi_SP[c(7,4,1,6,5,3,2)],2)), SITE= levels(set.ag2$SITE)))
pacf_seas_SF <- data.frame(list(phi= paste0("phi==",round(phi_SF[c(7,4,1,6,5,3,2)],2)), SITE= levels(set.ag2$SITE)))

