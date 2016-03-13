### load necessary packages 
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan")

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
         M1 = cos(2*pi*month_ret/26),
         M2 = cos(2*pi*month_ret/26),
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

site <- 1:7
data_mat2 <- expand.grid(site=site_levels,biweek=1:26, YEAR= 1990:2015) %>%
              mutate(site=factor(site,levels= site_levels),
                               biweek_year  = biweek+(YEAR-1990)*26,
                               MID = 1:nrow(data_mat)) %>%
            join(set.monyr)  

obs_SF <- matrix(set.sum2$mean_SF, ncol= 7, byrow= T)
obs_SP <- matrix(set.sum2$mean_SP, ncol= 7, byrow= T)



                                                                                   
