
### load necessary packages 
package.list<-c("abind","AER","bitops","car","chron","coda","colorspace","dichromat","digest","Formula",
                "gdata","ggplot2","gpclib","gtable","gtools","Hmisc","labeling","lmtest","lubridate",
                "memoise","munsell","mvtnorm","ncdf","plyr","proto","R.methodsS3","R.oo","R.utils","R2jags",
                "R2WinBUGS","RColorBrewer","RCurl","reshape2","rjags","sandwich","scales","sp",
                "stringr","timeDate","zoo","maps","maptools","classInt","rgeos","fields","rgdal")

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

### get means 
set.ag <- ddply(settlement, .(DATE_RETRIEVED,SITE,Duration),summarize, mean_sum= mean(S_FRANCISCANUS+S_PURPURATUS), mean_tot= mean(TOTAL_URCHINS))

### format data to estimate monthly means per brush per day
settlement$SITE <- factor(settlement$SITE,levels= unique(settlement$SITE))
settlement <- settlement[order(settlement$SITE,settlement$julian),]
settlement$diff <- as.numeric(settlement$Duration)
settlement <- drop.levels(settlement)
settlement<- subset(settlement,Duration>0)
set.sum <- ddply(settlement,.(SITE,DATE_RETRIEVED,DATE_DEPLOYED,month_ret,year_ret,Duration),summarise, SP= sum(S_PURPURATUS),SF= sum(S_FRANCISCANUS),TOT= sum(TOTAL_URCHINS),NB= length(TOTAL_URCHINS))
set.sum$ID <- 1:nrow(set.sum)

set.sum$SP_EM <- with(set.sum,ifelse(SP+SF==0,0,SP/(SP+SF)*TOT))
set.sum$SF_EM <- with(set.sum,ifelse(SP+SF==0,0,SF/(SP+SF)*TOT)/Duration)
set.sum <- subset(set.sum, !is.na(SP)&!is.na(SF)&Duration<45)
set.sum$SITE_NUM <- as.numeric(set.sum$SITE)
set.sum$M1 <- cos(2*pi*set.sum$month_ret/12)
set.sum$M2 <- cos(2*pi*set.sum$month_ret/12)
set.sum$monyr <- set.sum$month_ret+(set.sum$year_ret-1990)*12
head(set.sum)

site_levels <- c("Anacapa[SB]","Fort Bragg[NorCal]","Gaviota[SB]","Ocean Beach[SD]","Ellwood[SB]","Stearns Wharf[SB]","Scripps[SD]")

year <- seq(from= 0,to= diff(range(set.sum$year_ret)),by= 1)
month <- 1:12
site <- 1:7
data_mat <- expand.grid(site=site,month_ret=month,YEAR=year)
data_mat$SITE <- factor(data_mat$site,labels= levels(set.sum$SITE))
data_mat$MID <- 1:nrow(data_mat)
data_mat$year_ret <- data_mat$YEAR+1990
data_mat$monyr <- data_mat$month+(data_mat$YEAR)*12
data_mat$month <- data_mat$month_ret

set.monyr <- ddply(set.sum,.(monyr,SITE),summarise, mean= mean(SP))
set.monyr2 <- join(data_mat,set.monyr)
monyr.mat <- dcast(set.monyr2[,c("monyr","SITE","mean")],monyr~SITE, fun = mean)[,-1]

