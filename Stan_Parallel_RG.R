######################################################
###  Script to format data for analysis and plots  ###
###  Author:  D.K. Okamoto                         ###
######################################################

### load necessary packages 
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan",
                "dplyr")

lapply(package.list,library,character.only = T)


 ### load the covariates 
 covariates <- read.csv("1_Data/bayes_covariates.csv")[,-1]

 df <- covariates%>%
   filter(SITE=="GAVIOTA")%>%
   select(-one_of("biweek","biweek_year","species","SITE"))
 
 a <-factanal(df,3)
 
 
 min_biweek <- min(covariates$biweek_year)
 max_biweek <- max(covariates$biweek_year)

  ModMat <-model.matrix(~SITE:species +
                            SITE:species:MEI_runlag2 +
                            SITE:species:chla_rollmean_30 +
                            SITE:species:SST_rollmean_30 +
                            SITE:species:BK_stand_30 +
                            SITE:species:kelp_biomass +
                            SITE:species:adults-1, data= covariates)
  
  vars <- c("MEI_runlag2","chla_rollmean_30","SST_rollmean_30",
    "BK_stand_30","kelp_biomass","adults")
  
  mods <- ldply(0:6,function(x) t(combn(vars,x)))
  
  ModList <- list()
  VarList <- list()
  for(i in 1:nrow(mods)){
    VarList[[i]] <- mods[i,][!is.na(mods[i,])]
    cols_sub <- as.numeric(sapply(VarList[[i]],
      function(x)grep(x,dimnames(ModMat)[[2]])))
    ModList[[i]] <- ModMat[,c(1:14,cols_sub)]
  }
  
  sort_array <- matrix(1:nrow(covariates))
  
  dim(sort_array) <- c(498,7,2)

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
      dplyr::select(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration,SP,SF,TOT)%>%
      mutate(SP_EM= as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SP/(SP+SF)*TOT)/Duration)),
             SF_EM = as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration)),
             julian = julian(DATE_RETRIEVED, origin = "1990-01-01"),
             SITE_NUM = as.numeric(SITE),
             yday = yday(DATE_RETRIEVED),
             biweek  = floor(yday/14)+1) %>%
      mutate(biweek =ifelse(biweek>26,26,biweek),
             biweek_year  = biweek+(year_ret-1990)*26)


  set.sum2 <- settlement %>% group_by(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration) %>%
      mutate(SP= S_PURPURATUS,
             SF= M_FRANCISCANUS,
             TOT= TOTAL_URCHINS) %>%
      dplyr::select(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration,SP,SF,TOT)%>%
      mutate(SP_EM= as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SP/(SP+SF)*TOT)/Duration)),
             SF_EM = as.numeric(ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration)),
             julian = julian(DATE_RETRIEVED, origin = "1990-01-01"),
             SITE_NUM = as.numeric(SITE),
             yday = yday(DATE_RETRIEVED),
             biweek  = floor(yday/14)+1) %>%
      mutate(biweek =ifelse(biweek>26,26,biweek),
             biweek_year  = biweek+(year_ret-1990)*26) %>%
      filter(biweek_year>=min_biweek&biweek_year<=max_biweek)%>%
      mutate(biweek_year=biweek_year-min_biweek+1)

  D_site_star <-as.matrix(dist(unique(covariates$biweek_year),upper= T,diag= T))
  yearly <- seq(10,max(unique(covariates$biweek_year)),by = 26)
  m= length(yearly)
  n= length(unique(unique(covariates$biweek_year)))
  D_star2 <-as.matrix(dist(yearly,upper= T,diag= T))
  D_site_star2 <- as.matrix(dist(rbind(data.frame(x=unique(covariates$biweek_year)),
                                       data.frame(x=yearly))))[1:n, (n + 1):(n +m)]

  ### define GP covariance matrices

  phi  =1.5
  phi2 = 0.75

  NYM= 26
  Cstar = array(NA,dim= c(NYM,NYM))

  NT= dim(sort_array)[1]
  Cstar2 = array(NA,dim= c(NT,NT))
  Cstar3 = array(NA,dim= c(NT,NT))

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

  for (i in 1:(NT-1)) {
      for (j in (i + 1):(NT)) {
          Cstar3[i, j] =   (1+(3^(0.5)*abs(D_site_star[i,j]/26)/(phi2)))*exp(-(3^(0.5)*abs(D_site_star[i, j]/26)/(phi2)));
          Cstar3[j, i] = Cstar3[i, j];
      }
  }
  for (k in 1:(NT)) Cstar3[k, k] = 1+0.001;

  z <- rnorm(NT,0,1)
  plot(t(chol(Cstar2))%*%z,type= "l")
  lines(t(chol(Cstar3))%*%z,type= "l",col= "red")
  abline(v= seq(0,702,by= 26))

  ### list of data for model
  Data_List <- list()
  for(i in 1:length(ModList)){
    Data_List[[i]] <- with(set.sum2,list(
                                NO= nrow(set.sum2),
                                YSP = SP,
                                NB = rep(1,nrow(set.sum2)),
                                NS = length(unique(SITE)),
                                NM = dim(sort_array)[1],
                                NT= m,
                                NYM = 26,
                                NP= length(VarList[[i]]),
                                OBS_MONTH= biweek_year,
                                PRED_MONTH= covariates$biweek[1:dim(sort_array)[1]],
                                OBS_SITE= as.numeric(SITE),
                                X = ModList[[i]],
                                xind = sort_array,
                                D_seas = t(chol(Cstar)),
                                D= Duration,
                                n= SP+SF,
                                N= TOT
                            ))  
  }

  stan_parallel_fun <- function(x,n.iter=1000,n.burnin= 500,n.chains= 1){
    stan(file= "3_Analysis_Code/Settlement_Reg.stan", 
                  data=x,
                  iter =n.iter , warmup= n.burnin,
                  chains =1,cores = 1,
                  verbose = TRUE,init="random",
                  control= list(adapt_delta= 0.95,
                                max_treedepth=15))
  }

  mod_out <- mclapply(Data_List[rep(33:64,each= 3)],
                      stan_parallel_fun,mc.cores= 48)
  
  
iwaic_fun <- 
  
  
  