library(rstan);library(ggplot2);library(MASS);library(gdata)
rm(list=ls())
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan")

lapply(package.list,library,character.only = T)

### GET PDO, NPGO, MEI & NINO values from web ###
PDO <- read.table("http://jisao.washington.edu/pdo/PDO.latest", skip=31, nrows=116,fill=T)
PDO$V1 <- 1900:2015
NPGO <- read.table("http://www.o3d.org/npgo/npgo.php", skip=24,nrows=787)
names(PDO) <- c("YEAR",1:12)
PDO <- melt(PDO, id.vars= c("YEAR"))
names(PDO) <- c("YEAR","MONTH","PDO")
names(NPGO) <- c("YEAR","MONTH","NPGO")


### NINO 3.4 or 4 index at NOAA

MEI <- read.table("http://www.esrl.noaa.gov/psd/enso/mei/table.html", skip = 13, fill = T, header= T, nrows= 66)
MEI.long <- melt(data.frame(MEI), id.vars= "YEAR",variable.name= "MONTH")
names(MEI.long) <- c("YEAR","MONTH","MEI")
levels(MEI.long$MONTH) <- c(1:12)
MEI.long$MONTH <- as.numeric(as.character(MEI.long$MONTH))
MEI.long <- subset(MEI.long,YEAR>1989)
MEI.long <- MEI.long[order(MEI.long$YEAR,MEI.long$MONTH),]
MEI.long$MONTH <- as.character(as.numeric(MEI.long$MONTH))
MEI.long$monyr <-(MEI.long$YEAR-1990)*12+as.numeric(MEI.long$MONTH)
MEI.long$MEI_run[2:nrow(MEI.long)] <- rollmean(MEI.long$MEI,k= 2)

add_lag <- function(x,var="MEI") {
  xlag <- x[,c("YEAR","MONTH",var)]
  xlag$MONTH <- ifelse(as.numeric(xlag$MONTH)==12,1,as.numeric(xlag$MONTH)+1)
  xlag$YEAR <- ifelse(as.numeric(xlag$MONTH)==12, xlag$YEAR+1, xlag$YEAR)
  names(xlag)[3] <- paste0(var,"lag")
  xlag
}

add_lag2 <- function(x,var) {
  xlag <- x[,c("YEAR","MONTH","SITE",var)]
  xlag$MONTH <- ifelse(as.numeric(xlag$MONTH)==12,1,as.numeric(xlag$MONTH)+1)
  xlag$YEAR <- ifelse(as.numeric(xlag$MONTH)==12, xlag$YEAR+1, xlag$YEAR)
  names(xlag)[4] <- paste0(var,"lag")
  xlag
}

MEI <- join(MEI.long,add_lag(MEI.long,var= "MEI"))
MEI$YEAR2 <- with(MEI,ifelse(MONTH%in%c(11,12),YEAR+1,YEAR))
PDO <- join(PDO,add_lag(PDO,var= "PDO"))
PDO$YEAR2 <- with(PDO,ifelse(MONTH%in%c(11,12),YEAR+1,YEAR))
NPGO <- join(NPGO,add_lag(NPGO,var= "NPGO"))
NPGO$YEAR2 <- with(NPGO,ifelse(MONTH%in%c(11,12),YEAR+1,YEAR))

glob_vars <- join(MEI,PDO)
glob_vars <- join(glob_vars,NPGO)
head(glob_vars)
glob_vars$YEAR2 <- with(glob_vars,ifelse(MONTH%in%c(11,12),YEAR+1,YEAR))

