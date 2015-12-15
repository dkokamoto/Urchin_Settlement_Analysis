
### load necessary packages 
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan")

lapply(package.list,library,character.only = T)

### load urchin data ###
settlement <- read.csv("~/Copy/UrchinAnalyses/Data/Invertebrate_Settlement_All_Years.csv",header=T)

### format dates, etc.
settlement$DATE_RETRIEVED <- parse_date_time(settlement$DATE_RETRIEVED, "%Y-%m-%d")
settlement$DATE_DEPLOYED<- parse_date_time(settlement$DATE_DEPLOYED, "%Y-%m-%d")

settlement$DATE_RETRIEVED <- parse_date_time(settlement$DATE_RETRIEVED, "%Y-%m-%d")
settlement$DATE_DEPLOYED<- parse_date_time(settlement$DATE_DEPLOYED, "%Y-%m-%d")
settlement$julian <- as.numeric(julian(settlement$DATE_RETRIEVED))
settlement$ID <- with(settlement,paste(DATE_DEPLOYED,DATE_RETRIEVED,SITE),sep= "_")
settlement$ID2 <- with(settlement,paste(DATE_DEPLOYED,SITE),sep= "_")
settlement$ID3 <- with(settlement,paste(DATE_RETRIEVED,SITE),sep= "_")
settlement$month_dep <- month(settlement$DATE_DEPLOYED)
settlement$month_ret <- month(settlement$DATE_RETRIEVED)
settlement$year_dep <- year(settlement$DATE_DEPLOYED)
settlement$year_ret <- year(settlement$DATE_RETRIEVED)
settlement$day_dep <- day(settlement$DATE_DEPLOYED)
settlement$day_ret <- day(settlement$DATE_RETRIEVED)
settlement$S_FRANCISCANUS[settlement$year_ret==2009&settlement$month_ret==04&settlement$day_ret==07&is.na(settlement$S_FRANCISCANUS)&!is.na(settlement$S_PURPURATUS)] <- 0
settlement <- subset(settlement,!is.na(S_FRANCISCANUS)&!is.na(S_PURPURATUS)&!is.na(TOTAL_URCHINS))

### format data to estimate collection period means per brush per day
settlement$SITE <- factor(settlement$SITE,levels= unique(settlement$SITE))
settlement <- settlement[order(settlement$SITE,settlement$julian),]
settlement$diff <- as.numeric(settlement$Duration)
settlement <- drop.levels(settlement)
settlement<- subset(settlement,Duration>0)
set.sum <- ddply(settlement,.(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration),summarise, SP= sum(S_PURPURATUS),SF= sum(S_FRANCISCANUS),TOT= sum(TOTAL_URCHINS),NB= length(TOTAL_URCHINS))


set.sum$SP_EM <- with(set.sum,ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SP/(SP+SF)*TOT)/Duration/NB))
set.sum$SF_EM <- with(set.sum,ifelse(SP+SF==0&TOT>0,NA, ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration/NB))
set.sum$julian <- julian(set.sum$DATE_RETRIEVED, origin = "1990-01-01")
test2 <- ddply(set.sum,.(SITE), summarize, diff= diff(julian))
set.sum <- subset(set.sum, !is.na(TOT)&Duration<45)
set.sum$SITE_NUM <- as.numeric(set.sum$SITE)
set.sum$M1 <- cos(2*pi*set.sum$month_ret/26)
set.sum$M2 <- cos(2*pi*set.sum$month_ret/26)
set.sum$yday <- yday(set.sum$DATE_RETRIEVED)
set.sum$biweek  <- as.numeric(floor(set.sum$yday/14)+1)
set.sum$biweek  <- ifelse(set.sum$biweek>26,26,set.sum$biweek)
set.sum$biweek_year  <- set.sum$biweek+(set.sum$year_ret-1990)*26
site_levels <- c("Anacapa[SB]","Fort Bragg[NorCal]","Gaviota[SB]","Ocean Beach[SD]","Ellwood[SB]","Stearns Wharf[SB]","Scripps[SD]")

set.monyr <- ddply(set.sum,.(biweek_year,SITE),summarise, mean_SP= mean( SP_EM),mean_SF= mean( SF_EM))

site <- 1:7
data_mat <- expand.grid(site=site,biweek=1:26, YEAR= 1990:2015)
data_mat$SITE <- factor(data_mat$site,labels= levels(set.sum$SITE))
data_mat$biweek_year  <- data_mat$biweek+(data_mat$YEAR-1990)*26
data_mat$MID <- 1:nrow(data_mat)

set.sum2 <- join(data_mat,set.monyr) 


                                                                                   
