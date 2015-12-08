library(rstan);library(ggplot2);library(MASS);library(gdata)
rm(list=ls())
package.list<-c("abind","car","gdata","ggplot2","Hmisc","labeling","lubridate",
                "mvtnorm","plyr","RColorBrewer","reshape2","scales","sp","rstan")

lapply(package.list,library,character.only = T)
### now run "Stan_model_runs.R" ###
source("./3_Analysis_Code/Stan_model_runs.R")
load("~/Copy/UrchinAnalyses/Data/postAR.RData") 

### generate monthly mean settlement from the data
set.ag <- ddply(set.sum,.(monyr,SITE,month_ret,year_ret),summarise, mean_SP = mean(SP_EM, na.rm= T),mean_SF = mean(SF_EM, na.rm= T))



params <- extract(postAR)

### settlement partial autocorrelation function
phi_SP <-  apply(params$phi_SP, c(2), mean)
phi_SF <- apply(params$phi_SF, c(2), mean)

### estimated mean (log-scale) purps and frans from the Stan posterior
E_SP <- apply(exp(params$LSP), c(2,3), mean)

E_SF <- apply(exp(params$LSF), c(2,3), mean)

### get monthly mean (seasonal) estimates
beta_SP <- apply(exp(params$beta_SP), c(2), mean)
beta_SF <- apply(exp(params$beta_SF), c(2), mean)

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
seas_mod_SF$month <- c(3:12,1:2)
seas_mod_SP$month <- c(3:12,1:2)

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
SP_mod_df <- melt(seas_mod_SP, id.vars= "month",variable.name= "SITE",value.name= "Exp_seas_SP")
SP_mod_df$SP_U <- beta_SP_U 
SP_mod_df$SP_L <- beta_SP_L
SF_mod_df <- melt(seas_mod_SF, id.vars= "month",variable.name= "SITE",value.name= "Exp_seas_SF")
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
  x[is.na(monyr.mat)] <- NA
  x <- data.frame(x)
  names(x) <- levels(set.sum$SITE)
  x$monyr <- 1:nrow(x)
  x <- melt(x, id.vars= "monyr",value_name= name)
  names(x) <- c("monyr","SITE",name)
  return(x)
}

SP_mod <- format.data(SP_mod,"Est_SP")
SF_mod <- format.data(SF_mod,"Est_SF")
SP_mod_U <- format.data(SP_mod_U,"SP_U")
SP_mod_L <- format.data(SP_mod_L, "SP_L")
SF_mod_U <- format.data(SF_mod_U, "SF_U")
SF_mod_L <- format.data(SF_mod_L, "SF_L")

mod_df <- join_all(list(SP_mod,SF_mod,SP_mod_U,SP_mod_L,SF_mod_U,SF_mod_L), by= c("monyr","SITE"))
set.ag2 <- join(set.ag,mod_df, by = c("SITE","monyr"))
set.ag2$MONTH <- set.ag2$month_ret
set.ag2$YEAR <- set.ag2$year_ret

set_summary <- with(set.ag2, data.frame(expand.grid(SITE=levels(SITE),MONTH= min(MONTH):max(MONTH),YEAR= min(YEAR):max(YEAR))))


set.ag2$month_ret2 <- with(set.ag2,ifelse(month_ret==11,1,ifelse(month_ret==12,2,month_ret+2)))
set_summary <- join(set_summary, set.ag2)
set_summary$date <- with(set_summary,parse_date_time(paste(15,MONTH,YEAR,sep= "-"),"%d-%m-%y"))
levels(set.ag2$SITE) <- levels(seas_mod$SITE)[c(5,1,2,6,3,4,7)]
set.ag2$SITE <-factor(set.ag2$SITE, levels = rev(levels(set.ag2$SITE))[c(7,4,1,6,5,3,2)])   
set_summary$monyr <- set_summary$MONTH+(set_summary$YEAR-1990)*12
set_summary$SF_U <- with(set_summary,ifelse(SF_U>0.5,mean_SF*1.2, SF_U))



pacf_seas_SP <- data.frame(list(phi= paste0("phi==",round(phi_SP[c(7,4,1,6,5,3,2)],2)), SITE= levels(set.ag2$SITE)))
pacf_seas_SF <- data.frame(list(phi= paste0("phi==",round(phi_SF[c(7,4,1,6,5,3,2)],2)), SITE= levels(set.ag2$SITE)))


options <- theme(strip.text.y = element_text(size =12,angle= 0),
                 axis.text.y = element_text( colour= "black",size =12),
                 axis.text.x = element_text(colour= "black",size =12),
                 axis.title.y = element_text( colour= "black",size =12),
                 axis.title.x = element_text(colour= "black",size =12),
                 panel.grid.minor = element_line(colour = NA),
                 panel.grid.major = element_line(colour = NA),
                 panel.margin= unit(.5,"lines"),
                 strip.background = element_rect(fill= NA,colour= NA,size= 0),
                 legend.key = element_rect(colour=NA),
                 panel.grid = element_line(colour = NA),
                 legend.position = "right",
                 panel.background=element_rect(fill= NA,colour=NA),
                 plot.background=element_rect(fill= NA,colour=NA),
                 legend.direction="vertical",
                 legend.background=element_rect(fill= "NA"),
                 legend.key.height= unit(1,"lines"),
                 legend.title = element_text(size = 0),
                 legend.text = element_text(size = 14),
                 legend.title.align = 0,
                 panel.border= element_rect(fill= NA),
                 legend.key.width = unit(3,"lines"),
                 plot.margin = unit(rep(0.2, 4), "inches"))

ggplot(aes(date, mean_SF),data=set_summary)+
  geom_ribbon(aes(ymin= SF_L, ymax= SF_U),fill= "grey70")+
  geom_point(size= 1)+
  geom_path(aes(y= Est_SF),colour= "red")+
  theme_bw()+
  facet_wrap(~SITE,scales= 'free_y')+
  options
  
ggplot(aes(date, mean_SP),data=subset(set_summary,YEAR%in%c(2010:2012)))+
  geom_ribbon(aes(ymin= SP_L, ymax= SP_U),fill= "grey50")+
  geom_point(size= 1)+
  geom_path(aes(y= Est_SP),colour= "red",size= 0.25)+
  theme_bw()+
  facet_wrap(~SITE,scales= 'free_y')

Season_patterns <- ggplot(aes(month,Exp_seas_SP),data= seas_mod)+
  geom_ribbon(aes(ymax= SP_U, ymin= SP_L), fill= "grey80")+
  geom_line()+
  facet_wrap(~SITE,scales= "free_y",ncol=1)+
  geom_point(aes(month_ret2,mean_SP),data= set.ag2,size= 1, shape= 21, fill= "grey",alpha= 0.5)+
  ylab(expression(paste("monthly settlement (# ",brush^-1," ",day^-1,")")))+
  xlab("")+
 # geom_text(aes(x= 10,y= 5,label = phi),parse= TRUE,data= pacf_seas,size =4)+
  scale_y_continuous(trans= "log10")+
  annotation_logticks(side= "l")+
  theme_bw()+
  options+
  scale_x_continuous(breaks= c(2,4,6,8,10,12),labels= c("Dec","Feb","Apr","Jun",
                                                        "Aug","Oct"),expand= c(0,0))

pdf(width=3.5,height= 8.5,file= "~/Copy/UrchinAnalyses/Figures/Seas_patterns.pdf",family = "serif",pointsize = 16)
Season_patterns
dev.off()

sp_res <- SP_mod[,1:7]-seas_mod[rep(1:12,26),1:7]
sp_res$MONTH <- rep(1:12,26)
sp_res$YEAR <- rep(1990:2015,each= 12)
sp_res_df <- melt(sp_res, id.vars= c("MONTH","YEAR"),variable.name= "SITE")
sp_res_df$monyr <- with(sp_res_df,MONTH+(YEAR-1990)*12)
sp_res_df <- join(sp_res_df,SP_mod_df)


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

### local measurements 
### get formatted temp, transport and chlorophyll concentration data ###	
SST.TS <- read.csv("~/Copy/UrchinAnalyses/Data/SST.TS_10.12.15.csv",header=T)
BAKUN.TS <- read.csv("~/Copy/UrchinAnalyses/Data/BAKUN.TS_10.12.15.csv",header=T)
CHL.TS <- read.csv("/Users/Dan/Data/Urchins/Variables/CHL/CHL.TS-9.6.14.csv",header= T)
LOC_VARS <- join(SST.TS,CHL.TS[,c("chla","julian","site")])
LOC_VARS$Region2 <- ifelse(LOC_VARS$site%in%c("FBPC"), "NORTH","SCB")
LOC_VARS <- join(LOC_VARS,BAKUN.TS)

means <- aggregate(LOC_VARS[,-c(1:7,17)], 
                  by= list(MONTH= LOC_VARS$month,
                           YEAR=LOC_VARS$year,SITE=LOC_VARS$site), 
                  FUN= mean, na.rm= TRUE)

LOC_VARS_std <- ddply(means, .(MONTH,SITE), summarize,YEAR= YEAR,
                  SST_ANOM= SST-mean(SST,na.rm= T),
                  ST_ANOM_30= SST_rollmean_30-mean(SST_rollmean_30,na.rm= T),
                  chla_anom= chla-mean(chla,na.rm= T), 
                  BK2_anom = BAKUN2- mean(BAKUN2,na.rm=T),
                  BK_anom = BAKUN- mean(BAKUN,na.rm=T))

LOC_VARS_std$BK_anom[is.nan(LOC_VARS_std$BK_anom)] <- NA

LOC_VARS_std <- LOC_VARS_std[order(LOC_VARS_std$SITE,LOC_VARS_std$YEAR,LOC_VARS_std$MONTH),]

LOC_VARS_std <- join(LOC_VARS_std,add_lag2(LOC_VARS_std,var= "chla_anom"))
LOC_VARS_std <- join(LOC_VARS_std,add_lag2(LOC_VARS_std,var= "SST_ANOM"))
LOC_VARS_std <- join(LOC_VARS_std,add_lag2(LOC_VARS_std,var= "BK_anom"))
LOC_VARS_std$YEAR2 <- with(LOC_VARS_std,ifelse(MONTH%in%c(11,12),YEAR+1,YEAR))

sp_monthly <- join(sp_res_df,glob_vars, by= c("YEAR","MONTH"))
sp_monthly <- join(sp_monthly,LOC_VARS_std, by= c("YEAR","MONTH","SITE"))
varlist <- c("MEIlag","PDOlag","NPGOlag","SST_ANOMlag","chla_anomlag","BK_anomlag")

b <- subset(sp_monthly,SITE%in%c("GAVIOTA")&!is.na(value)&!is.na(MEIlag)&!is.na(SST_ANOMlag)&!is.na(chla_anomlag)&!is.na(BK_anomlag)&MONTH%in%c(1:6))
b <- b[,c("value",varlist)]
names(b) <- c("set","ENSO","PDO","NPGO","SST","chla","BK")
b <- data.frame(apply(b,2,function(x)x/sd(x)))

### create the blacklist 
df1 <- data.frame(list(from = "set",to = rep(names(b)[-1])))
df2 <- data.frame(list(from= "chla",to= c("ENSO","PDO","NPGO","BK","SST")))
df3 <- data.frame(list(from= "SST",to= c("ENSO","PDO","NPGO","BK")))
df4 <- data.frame(list(from= "BK",to= c("ENSO","PDO","NPGO")))
df5 <- data.frame(list(from= "PDO",to= c("ENSO","NPGO")))
df6 <- data.frame(list(from= "NPGO",to= c("ENSO","PDO")))
df7 <- data.frame(list(from= "ENSO",to= c("NPGO","PDO")))
blk <- rbind(df1,df2,df3,df4,df5,df6,df7)

### create the whitelist 
bnSB <-hc(b,blacklist= blk)
bnSB3 <- arc.strength(bnSB,b,criterion= "cor")
strength.plot(bnSB,bnSB3)
test5 <- boot.strength(b,algorithm= "hc",algorithm.args= list(blacklist= blk),R= 1000,cpdag= FALSE)
bnSB2 <- bn.fit(bnSB,b)

plot(bnSB)

ci.test("set","ENSO",c("NPGO","SST"),data= b)
ci.test("set","SST",c("NPGO","ENSO"),data= b)
ci.test("set","NPGO",c("SST","ENSO"),data= b)
ci.test("BK","PDO",data= b)
ci.test("SST","ENSO",c("BK"),data= b)
ci.test("SST","BK",c("ENSO"),data= b)
ci.test("chla","SST",data= b)










bnSB <-tabu(b,blacklist= df,score= "bde")
test6 <- cpdist(bnSB2,nodes= names(b))

test5 <- boot.strength(b,algorithm= "tabu",algorithm.args= list(blacklist= df1),R= 200,cpdag= FALSE)


pdf(width = 4, height =6, file = "~/Copy/UrchinAnalyses/Figures/GAVIOTA_Bayesian_network.pdf", family = "Times",pointsize = 14)
strength.plot(bnSB,bnSB3)
dev.off()

bnSB2 <- boot.strength(b[,c("value",varlist)], algorithm = "hc")
strength <- custom.strength(list(bnSB), nodes = names(b[,c("value",varlist)]),
                cpdag = FALSE)
bn.SB <- set.arc(bnSB, "value","chla_anom")
arc.strength(bn.SB, b[,c("value",varlist)])

a <- glm(exp(MEAN)~MEI,data= subset(sp_monthly,SITE== "GAVIOTA"),family= gaussian(link= "log"))
b2 <- bn.fit(b,a)
b$learning


sp_monthly2 <- subset(sp_monthly,MONTH%in%c(1:6))
sp_ann4 <- aggregate(sp_monthly2[,c("value","MEIlag")], by= list(YEAR= sp_monthly2$YEAR2,SITE= sp_monthly2$SITE), FUN= mean, na.rm= T)

sp_monthly3 <- subset(sp_monthly,MONTH%in%c(1:3))
sp_ann5 <- aggregate(sp_monthly3[,c("value","MEIlag")], by= list(YEAR= sp_monthly3$YEAR2,SITE= sp_monthly3$SITE),FUN= mean)

names(sp_ann5)[4] <- "MEI_win" 

sp_ann4 <- join(sp_ann4,sp_ann5[,c("YEAR","SITE","MEI_win")])


sp_ann4 <- subset(sp_ann4,!is.na(value)&!is.na(MEI_win))
sp_ann4$st_val <- ddply(sp_ann4,.(SITE),summarize,st_val= exp(value)/max(exp(value),na.rm=T))$st_val

fit1 <- rq(log(st_val)~MEI_win:SITE+SITE,data=sp_ann4,tau=c(0.14,0.86))
sp_ann4$pred_l <- exp(predict(fit1)[,1])
sp_ann4$pred_u <- exp(predict(fit1)[,2])
sp_ann4$SITE <-factor(sp_ann4$SITE, levels= levels(sp_ann4$SITE)[c(2,1,3,5,6,4,7)]) 
levels(sp_ann4$SITE) <- c("Ft Bragg [NorCal]","Anacapa [SB]","Gaviota [SB]","Ellwood [SB]","Stearrns Wharf [SB]","Ocean Beach [SD]","Scripps [SD]")

pdf(width = 13, height =2.8, file = "~/Copy/UrchinAnalyses/Figures/MEI_annual.pdf", family = "Times",pointsize = 20)                  
ggplot(aes(MEIlag,st_val),dat= sp_ann4)+
  geom_point(shape= 21,fill= "grey")+
  facet_wrap(~SITE,ncol= 7,scales= "free_y")+
  stat_smooth(method= "glm",colour= "black",family= gaussian(link= "log"),se= FALSE)+
  geom_line(aes(MEI_win,pred_l), linetype= "dotted")+
  geom_line(aes(MEI_win,pred_u), linetype= "dotted")+
  theme_bw()+
  scale_y_continuous(breaks= c(0,0.5,1))+
  options+theme(axis.text.y= element_text(angle =90,hjust= 0.5))+
  ylab("Scaled Annual\nSettlement Index")+
  xlab("Multivariate ENSO Index")
dev.off()


Omega2 <- cor(sp_res[,1:7],use= "pairwise.complete.obs")

sigma <- Fit.test$par[names(Fit.test$par)%in%paste("L_sigma_SP","[",rep(1:data$NS,each= data$NS),rep(1:data$NS,data$NS),"]",sep= "")]

Omega <- matrix(Fit.test$par[names(Fit.test$par)%in%
                               paste("Omega","[",rep(1:data$NS,data$NS),
                                     ",",rep(1:data$NS,each= data$NS),"]",sep= "")],
                ncol= data$NS,dimnames= list(levels(set.sum$SITE),levels(set.sum$SITE)))

rgb.palette <- colorRampPalette(c(brewer.pal(11,"RdYlBu")[1:5],"white"),space= "rgb",bias=1.3)
rgb.palette2 <- colorRampPalette(rev(c("white",brewer.pal(11,"RdYlBu")[7:11])),space= "rgb",bias=1.3)


Omega_df2 <- melt(Omega2)
Omega_df2[,1] <- factor(Omega_df2[,1])
Omega_df2[,2] <- factor(Omega_df2[,2])
levels(Omega_df2[,1]) <- site_levels
levels(Omega_df2[,2]) <- site_levels
Omega_df2[,2] <- factor(Omega_df2[,2],levels= levels(Omega_df2[,2])[c(7,4,1,6,5,3,2)])
Omega_df2[,1] <- factor(Omega_df2[,1],levels= levels(Omega_df2[2,1])[c(7,4,1,6,5,3,2)])

options <- theme(
  strip.text = element_text(size= 10),
  panel.grid.minor = element_line(colour = NA),
  panel.grid.major = element_line(colour = NA),
  strip.background = element_rect(fill= NA,colour= NA),
  legend.key = element_blank(),
  legend.key.height = unit(0.2, "cm"),
  legend.background = element_rect(fill = NA),
  legend.position = c(0.3,1.11),
  legend.title = element_text(size = 12,face= "plain"),
  legend.text = element_text(size = 8,face= "plain"),
  legend.title.align = 0.5,
  axis.ticks.length = unit(.25,"lines"),
  axis.title = element_blank(),
  axis.text.y= element_text(vjust=0.5,size=10,colour= "black"),
  axis.text.x= element_text(angle= 35,hjust=1,size=10,colour= "black"),
  axis.title.x = element_blank(),
  axis.ticks.margin = unit(.25,"lines"),
  panel.background=element_rect(fill= NA,colour=NA),
  plot.background=element_rect(fill= NA,colour=NA),
  legend.direction="horizontal",
  panel.border = element_rect(colour="black",fill=NA),
  legend.key.width = unit(3.5,"lines"),
  panel.margin= unit(5,"lines"),
  plot.margin = unit(c(0.5,0.1,0.1,0.1), "inches"))

corr_plot <- ggplot(aes(x=Var1,y=Var2),data= Omega_df2)+
  geom_tile(aes(fill= value))+  
  scale_x_discrete(expand= c(0,0))+
  scale_y_discrete(expand= c(0,0))+
  options+
  scale_fill_gradientn(colours= c(rgb.palette2(50),
                                  rev(rgb.palette(50))),
                       name= "correlation in monthly settlement anomalies",
                       guide= guide_colourbar(title.position= "top"),
                       limits= c(-1,1),
                       breaks =c(-10:10)/5)+
  geom_text(aes(label= round(value,2)), size= 2)
corr_plot

pdf("../Figures/Resid_Corrs.pdf", width= 4.1,height =4)
corr_plot
dev.off()

seas_cor_df <- melt(seas_cor)
levels(seas_cor_df[,1]) <- site_levels
levels(seas_cor_df[,2]) <- site_levels
seas_cor_df[,2] <- factor(seas_cor_df[,2],levels= levels(seas_cor_df[,2])[c(7,4,1,6,5,3,2)])
seas_cor_df[,1] <- factor(seas_cor_df[,1],levels= levels(seas_cor_df[,1])[c(7,4,1,6,5,3,2)])

SSC <-  ggplot(aes(x=Var1,y=Var2),data= seas_cor_df)+
  geom_tile(aes(fill= value))+  
  scale_x_discrete(expand= c(0,0))+
  scale_y_discrete(expand= c(0,0))+
  options+
  scale_fill_gradientn(colours= c(rgb.palette2(50),
                                  rev(rgb.palette(50))),
                       name= "correlation in seasonal settlement trend",
                       guide= guide_colourbar(title.position= "top"),
                       limits= c(-1,1),
                       breaks =c(-10:10)/5)+
  geom_text(aes(label= round(value,2)), size= 2)

pdf("../Figures/Seas_Corrs2.pdf", width= 4.1,height =4)
SSC
dev.off()

year_anom[is.na(monyr.mat.year[-1,2:8])] <- NA

ann_cor_df <- melt(cor(year_anom,use="pairwise.complete.obs"))
ann_cor_df[,1] <- factor(ann_cor_df[,1])
ann_cor_df[,2] <- factor(ann_cor_df[,2])
levels(ann_cor_df[,1]) <- site_levels
levels(ann_cor_df[,2]) <- site_levels
ann_cor_df[,2] <- factor(ann_cor_df[,2],levels= levels(ann_cor_df[,2])[c(7,4,1,6,5,3,2)])
ann_cor_df[,1] <- factor(ann_cor_df[,1],levels= levels(ann_cor_df[,1])[c(7,4,1,6,5,3,2)])

ASC <- ggplot(aes(x=Var1,y=Var2),data= ann_cor_df)+
  geom_tile(aes(fill= value))+  
  scale_x_discrete(expand= c(0,0))+
  scale_y_discrete(expand= c(0,0))+
  options+
  scale_fill_gradientn(colours= c(rgb.palette2(50),rev(rgb.palette(50))),name= "correlation in seasonal settlement patterns",guide= guide_colourbar(title.position= "top"),limits= c(-1,1),breaks =c(-5:5)/5)+
  geom_text(aes(label= round(value,2)));ASC

ggplot(aes(month,exp(pred)*14),data= a)+geom_point()+
  facet_wrap(~site)+theme_bw()

ggplot(aes(month_ret,exp(beta_SP)),data= a)+geom_point()+
  facet_wrap(~site)+theme_bw()+
  scale_y_continuous(trans= "log10")+annotation_logticks(side= "l")
  geom_point(aes(month_ret,exp(SP_mod))

plot(exp(SP_mod),data$YSP)             
             
n.chains =3
n.iter =1500
n.burnin=500
set.seed <- 1234
params <-c("LSP","LSF")
Fit.Bayes = stan(model_code = model,
                data=data,
                pars=params,
                iter = n.iter, warmup= n.burnin,chains = n.chains,
                verbose = FALSE,init="random",seed= set.seed)
a <- extract(Fit.Bayes)
### plot to see date and monthly mean estimates

df1 <- data.frame(list(SP_E= exp(Fit.test$par[names(Fit.test$par)%in%paste("LSP","[",1:data$NO,"]",sep= "")]),SF_E=exp(Fit.test$par[names(Fit.test$par)%in%paste("LSF","[",1:data$NO,"]",sep= "")]), ID=1:data$NO))

set.sum2 <- join(set.sum,df1,by = "ID")
plot.set <- join(plot.set,set.ag)
plot(I(SP_E/diff/NB)~I(SP_EM/diff/NB),data= plot.set)


df1 <- expand.grid(list(RET_MONTH=1:12,RET_YEAR= 1991:2014,SITE= levels(set.sum$SITE)))
a <- ddply(plot.set,.(RET_MONTH,RET_YEAR,SITE),summarise, mean= mean(SP_E,na.rm=T),E_mean= mean(SP_EM,na.rm=T))
df1 <- join(df1,a)
b <- reshape(df1[,1:4], timevar = "SITE", idvar = c("RET_YEAR","RET_MONTH"), direction = "wide")
names(b)[3:9]<-gsub("mean.","",names(b)[3:9])
corrs <- rcorr(log(as.matrix(b[,3:9],ncol= 7)))
N <- adply(corrs$n,c(1,2))
names(N) <- c("Site","Ref_Site","N")
plot(log(b[,3:9]))

x<-acf(log(as.matrix(b[,3:9],ncol= 7)),na.action=na.pass)
a <- ddply(plot.set,.(RET_MONTH,RET_YEAR,SITE),summarise, mean= mean(SP))

x2<-adply(x$acf,c(1,2,3))
names(x2) <- c("Lag","Site","Ref_Site","r")
x2$Lag <- as.numeric(as.character(x2$Lag))
x2$Ref_Site <- factor(x2$Ref_Site,labels= x$snames)
x2$Site <- factor(x2$Site,labels= x$snames)
x2<- join(x2,N)
x2$sig <- with(x2,ifelse(abs(r)>=2/(sqrt(N)),1,0))
x2$SITE_LAB <- ifelse(x2$Ref_Site==x2$Site,x2$Site,"")
ggplot(aes(x=Lag,xend= Lag,y=r,yend=0),data= x2)+
  #geom_point(aes(colour= factor(sig),fill= factor(sig)),size= 2)+
  geom_segment(aes(colour= factor(sig),fill= factor(sig)),size= 1)+
  facet_grid(Ref_Site~Site)+plot.options+
  geom_hline(yintercept=0)+
  geom_text(aes(label= SITE_LAB,y= 0.8,x= 0.5),hjust= 0)+
  scale_colour_manual(values= c("grey","black"))+
  scale_fill_manual(values= c("grey","black"))

a <- a[order(a$RET_YEAR,a$RET_MONTH),]
b <- reshape(a[,c("SITE","RET_YEAR","RET_MONTH","mean")], timevar = "SITE", idvar = c("RET_YEAR","RET_MONTH"), direction = "wide")


acf(log(b[,3:9]),na.action= na.pass)

rcorr(as.matrix(b[,3:9],ncol=7))
)s.set3  <-s.set[!duplicated(s.set$ID2),c("SITE","RET_MONTH","RET_YEAR","ID2")]
cor.set <- join(s.set3,df2)

sp.cor.error <-  as.matrix(reshape(cor.set[,c("SITE","RET_MONTH","RET_YEAR","sp_error")], timevar = "SITE", idvar = c("RET_MONTH", "RET_YEAR"), direction = "wide")[,-c(1,2)],ncol=7))
colnames(sp.cor.error) <-levels(set$SITE) 
rcorr(sp.cor.set)

sp.cor.set <-  as.matrix(reshape(cor.set[,c("SITE","RET_MONTH","RET_YEAR","MMSP")], timevar = "SITE", idvar = c("RET_MONTH", "RET_YEAR"), direction = "wide")[,-c(1,2)],ncol=7)
colnames(sp.cor.set) <-levels(set$SITE) 
rcorr(sp.cor.set)


sf.cor.set <-  reshape(cor.set[,c("SITE","RET_MONTH","RET_YEAR","sf_error")], timevar = "SITE", idvar = c("RET_MONTH", "RET_YEAR"), direction = "wide")

rcorr(sp.cor.set)
sp.corr <- rcorr(as.matrix(sp.cor.set[,-c(1,2)],ncol= 7))


head(set)
ggplot(aes(RET_MONTH,SP_Empirical_Mean/diff),data= plot.set)+
  geom_line(aes(RET_MONTH,SP),col= "black")+
  geom_point(colour= "grey",size= 0.5)+
  #geom_point(aes(y=SF_Empirical_Mean/diff), colour= "red")+
  #geom_point(aes(julian,SF),col= "red")+
  facet_wrap(~SITE,scales= "free_y",ncol=1)+plot.options

head(set)
