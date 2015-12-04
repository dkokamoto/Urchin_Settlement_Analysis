
setwd("/Users/Dan/Copy/UrchinAnalyses/Data/")
rm(list=ls())
package.list<-c("abind","AER","bitops","car","chron","coda","colorspace","dichromat","digest","Formula",
				"gdata","ggplot2","gpclib","gtable","gtools","Hmisc","labeling","lmtest","lubridate",
				"memoise","munsell","mvtnorm","ncdf","plyr","proto","R.methodsS3","R.oo","R.utils","R2jags",
				"R2WinBUGS","RColorBrewer","RCurl","reshape2","rjags","sandwich","scales","sp",
				"stringr","timeDate","zoo","maps","maptools","classInt","rgeos","fields","rgdal")

### install packages if you don't have them 
### not run
lapply(package.list,library,character.only = T)

### get formatted temp, transport and chlorophyll concentration data ###	
SST.TS <- read.csv("/Users/Dan/Data/Urchins/Variables/SST/SST.TS-9.5.14.csv",header=T)
BAKUN.TS <- read.csv("/Users/Dan/Data/Urchins/Variables/BAKUN/BAKUN.TS-9.7.14.csv",header=T)
CHL.TS <- read.csv("/Users/Dan/Data/Urchins/Variables/CHL/CHL.TS-9.6.14.csv",header= T)
SST.TS2 <- join(SST.TS,CHL.TS[,c("chla","julian","site")])
for (i in 1:7){
  x <- na.approx(SST.TS2$chla[SST.TS2$site==levels(SST.TS2$site)[i]&SST.TS2$julian>=10227],na.rm= F)
  SST.TS2$chla2[SST.TS2$site==levels(SST.TS2$site)[i]&SST.TS2$julian>=10227] <- x
}
for (i in 1:7){
  x <- na.spline(SST.TS2$chla[SST.TS2$site==levels(SST.TS2$site)[i]&SST.TS2$julian>=10227],na.rm= F)
  SST.TS2$chla3[SST.TS2$site==levels(SST.TS2$site)[i]&SST.TS2$julian>=10227] <- x
}

  settlement <- read.csv("Invertebrate_Settlement_All_Years.csv",header=T)
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

  settlement$Region <- factor(ifelse(settlement$SITE %in%c("GAVIOTA","SBSTWRF","SBELL","ANACAPA"),
                                     "SB channel",ifelse(settlement$SITE%in%c("SIO","OCNBCH"),"SD","North")))
  
  settlement$Region2 <- ifelse(settlement$SITE=="FBPC","North","South")
  
  settlement$quarter <- ifelse(settlement$month_ret<4,1,ifelse(settlement$month_ret>9,4,ifelse(settlement$month_ret<7,2,3)))
  settlement$DOY <- yday(settlement$DATE_RETRIEVED)
  settlement$meanperday.purps<- settlement$S_PURPURATUS/(settlement$Duration)	
  settlement$meanperday.reds <- settlement$S_FRANCISCANUS/(settlement$Duration)	
  
  settlement$jul_ret <-julian(settlement$DATE_RETRIEVED)	
  
  settlement$jul_dep <- julian(settlement$DATE_DEPLOYED)
  settlement$site2 <- factor(settlement$SITE, levels= levels(settlement$SITE)[c(2,1,3,5,6,4,7)])
  levels(settlement$site2) <- c("Fort Bragg","Anacapa","Gaviota","Ellwood","Stearns","Ocean Beach","Scripps")
  settlement$jul.diff <- settlement$jul_ret-settlement$jul_dep
  
  settlement <- settlement[!(settlement$S_PURPURATUS<0),]
  
  set.ag6 <- with(settlement,aggregate(data.frame(list(meanperday.purps)),
                                       by = list(year_ret,SITE,month_ret,Region,quarter),FUN= mean,na.rm= T))
  
  ### name the 
  names(set.ag6) <- .(year,site,month,Region,quarter,SP.DAY)
  set.ag6$date <- with(set.ag6,parse_date_time(paste(year,month,"15",sep ="-"),"%Y-%m-%d"))
  set.ag6$DOY <- yday(set.ag6$date)
  set.ag6$Region2 <- ifelse(set.ag6$Region=="North","North","SCB")
  set.ag6$julian <- julian(set.ag6$date)
  set.ag6 <- set.ag6[order(set.ag6$site,set.ag6$year,set.ag6$month),]
  
  set.ag6$chla <-as.numeric(NA)
  set.ag6$chla_15 <-as.numeric(NA)
  set.ag6$chla_30 <-as.numeric(NA)
  set.ag6$chla_45 <-as.numeric(NA)
  
  for (i in 1:nrow(set.ag6)){
    a <- subset(CHL.TS,site==set.ag6$site[i]&julian<set.ag6$julian[i]&julian>(set.ag6$julian[i]-30))$chla
    set.ag6$chla[i] <- mean(a,na.rm= T)}
  
  for (i in 1:nrow(set.ag6)){
    a <- subset(CHL.TS,site==set.ag6$site[i]&julian<(set.ag6$julian[i]-15)&julian>(set.ag6$julian[i]-45))$chla
    set.ag6$chla_15[i] <- mean(a,na.rm= T)}
  
  
  for (i in 1:nrow(set.ag6)){
    a <- subset(CHL.TS,site==set.ag6$site[i]&julian<(set.ag6$julian[i]-30)&julian>(set.ag6$julian[i]-60))$chla
    set.ag6$chla_30[i] <- mean(a,na.rm= T)}
  
  for (i in 1:nrow(set.ag6)){
    a <- subset(CHL.TS,site==set.ag6$site[i]&julian<(set.ag6$julian[i]-45)&julian>(set.ag6$julian[i]-75))$chla
    set.ag6$chla_45[i] <- mean(a,na.rm= T)}
  
  
  set.ag6$Region2 <- ifelse(set.ag6$Region=="North","North","SCB")
  set.ag6 <- join(set.ag6,BAKUN.TS,by  =c("Region2","julian"))
  set.ag6 <- join(set.ag6,SST.TS,by  =c("site","julian"))
  
  ag.vars <- names(set.ag6)[grepl("BAKUN|SST|chla|SP.DAY",names(set.ag6))]
  
  set.ag9 <-  with(set.ag6,aggregate(set.ag6[,ag.vars], by = list(site=site,month=month,year=year),FUN= mean,na.rm= T))
  
  ### add overall monthly means for each site (mean seasonality)  ####
  set.ag7 <- with(set.ag6,aggregate(set.ag6[,ag.vars], by = list(site=site,month=month),FUN= mean,na.rm= T))
  names(set.ag7)[-c(1,2)] <- paste("m",ag.vars,sep= "_")
  
  ### add monthly standard deviation (variation about seasonality)  ####
  set.ag8 <- with(set.ag6,aggregate(set.ag6[,ag.vars], by = list(site=site,month=month),FUN= sd,na.rm= T))
  names(set.ag8)[-c(1,2)] <- paste("sd",ag.vars,sep= "_")
  set.joined <- join(set.ag9,set.ag7)
  set.joined <- join(set.joined,set.ag8)
  ### create empty dataframe 
  settlement2 <- data.frame(list(year= rep(1991:2014,each= 12), month= rep(1:12,24), 
                                 site = rep(levels(set.ag6$site), each = 24*12)))
  
  ### put values into the monthly dataframe ###
  set.monthly <- join(settlement2, set.joined)


	### GET PDO, NPGO, MEI values from web ###
 	PDO <- read.table("http://jisao.washington.edu/pdo/PDO.latest", skip=31, nrows=115,fill=T,header= T)[,-1]
names(PDO) <- 1:12
	PDO$YEAR <- 1901:2015
  PDO <- melt(PDO, id.vars= c("YEAR"))
  names(PDO) <- c("year","month","PDO")
	NPGO <- read.table("http://www.o3d.org/npgo/npgo.php", skip=24,nrows=785)
	names(NPGO) <- c("year","month","NPGO")
	
	PDOlag <- PDO
	PDOlag$PDOlag <- PDO$PDO
	PDOlag$month <- ifelse(as.numeric(PDO$month)==1,12,as.numeric(PDO$month)-1)
	PDO <- join(PDO,PDOlag[,names(PDOlag)[!(names(PDOlag)=="PDO")]])

	NPGOlag <- NPGO
	NPGOlag$NPGOlag <- NPGO$NPGO
	NPGOlag$month <- ifelse(NPGO$month==1,12,NPGO$month-1)
	NPGO <- join(NPGO,NPGOlag[,names(NPGOlag)[!(names(NPGOlag)=="NPGO")]])
	
	MEI <- read.table("http://www.esrl.noaa.gov/psd/enso/mei/table.html", skip = 13, fill = T, header= T, nrows= 66)
	MEI.long <- melt(data.frame(MEI), id.vars= "YEAR",variable.name= "month")
		names(MEI.long) <- c("year","month","MEI")
		levels(MEI.long$month) <- c(1:12)
		MEI.long$month <- as.character(as.numeric(MEI.long$month))
	MEIlag <- MEI.long
	MEIlag$MEIlag <- MEIlag$MEI
	MEIlag$month <- ifelse(as.numeric(MEI.long$month)==12,1,as.numeric(MEI.long$month)+1)
	MEI <- join(MEI.long,MEIlag[,names(MEIlag)[!(names(MEIlag)=="MEI")]])
	
  vars <- join(PDO,MEI)
  vars <- join(vars, NPGO)
  set.monthly <- join(set.monthly,vars)
  vars.monthly <- set.monthly[,!names(set.monthly)=="SP.DAY"]
  a <- melt(vars.monthly[,1:77],id = c("month","year","site"))
  b <- ddply(a, .(month, site,variable),summarise, stand =value-mean(value,na.rm= T), year= year)
  stand <- reshape(b,idvar= c("month","site","year"),timevar= "variable", direction= "wide")
  
  vars.monthly <- join(vars.monthly,stand)
