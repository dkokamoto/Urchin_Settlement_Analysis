#### requires you run the "Stan_model_runs.R" and "Analyses" First ###
setwd("/Users/Dan/Data/Urchins/SettlementProject/SettlementData")
package.list<-c("abind","AER","bitops","biwavelet","car","chron","coda","colorspace","dichromat","digest","Formula",
				"gdata","ggplot2","gpclib","nlme","gtable","gtools","Hmisc","labeling","lmtest","lubridate",
				"memoise","munsell","mvtnorm","ncdf","plyr","proto","R.methodsS3","R.oo","R.utils","R2jags",
				"R2WinBUGS","RColorBrewer","gridExtra","RCurl","reshape2","rjags","sandwich","scales","sp",
				"stringr","timeDate","zoo","maps","maptools","classInt","rgeos","fields","rgdal")

### install packages if you don't have them 
### not run
#lapply(package.list,library,character.only = T)
lapply(package.list,library,character.only=T)

### load data ###
set.gav <- subset(sp_monthly,SITE=="GAVIOTA")[,c("MEI","NPGO","chla_anom","PDO","value","BK_anom","BK2_anom","SST_ANOM")]
date <- with(subset(sp_monthly,SITE=="GAVIOTA"),parse_date_time(paste(YEAR,MONTH,"15",sep= "-"),"%Y-%m-%d"))
set.gav$value2 <- na.approx(set.gav$value, 1:nrow(set.gav),na.rm= FALSE)

### name rows ###
row.names(set.gav) <- 1:nrow(set.gav)
set.gav <- set.gav[14:307,]
date <- date[14:307]
### create the necessary variables ###
MEI <- cbind(1:nrow(set.gav),set.gav[,"MEI"])
urch <- cbind(1:nrow(set.gav),set.gav[,"value2"])
chl <- cbind(1:nrow(set.gav),set.gav[,"chla_anom"])
NPGO <- cbind(1:nrow(set.gav),set.gav[,"NPGO"])
PDO  <- cbind(1:nrow(set.gav),set.gav[,"PDO"])
BK <- cbind(1:nrow(set.gav),set.gav[,"BK_anom"])
BK2 <- cbind(1:nrow(set.gav),set.gav[,"BK2_anom"])
SST <- cbind(1:nrow(set.gav),set.gav[,"SST_ANOM"])

### wavelet coherence analysis for each potentially influential variable on the monthly scale ###
 pwtc.MEI <- wtc(urch,MEI,nrands=1000,dj= 1/36)
 pwtc.chl <- wtc(urch[84:269,],chl[84:269,],nrands=1000,dj= 1/36)
 pwtc.NPGO <- wtc(urch,NPGO,nrands=1000,dj= 1/36)
 pwtc.PDO<- wtc(urch,NPGO,nrands=1000,dj= 1/36)
 pwtc.BK <- wtc(urch,BK,nrands=0,dj= 1/36)
 pwtc.BK2 <- wtc(urch,BK2,nrands=0,dj= 1/36) 
 pwtc.SST <- wtc(urch,SST,nrands=0,dj= 1/36)

 plot(pwtc.SST)
 plot(pwtc.NPGO)
 plot(pwtc.MEI)

  options <- theme(strip.text.x = element_text(size =16),
			axis.text.y = element_text( colour= "black",size =12),
			axis.text.x = element_text(colour= "black",size =12),
			axis.text.y = element_text( colour= "black",size =12),
			axis.title.x = element_text(colour= NA,size= 0),
			panel.grid.minor = element_line(colour = NA),
			panel.grid.major = element_line(colour = NA),
			panel.margin= unit(.5,"lines"),
			strip.background = element_rect(fill= NA,colour= NA,size= 0),
			legend.key = element_rect(colour=NA),
			panel.grid = element_line(colour = NA),
			legend.position = "none",
			panel.background=element_rect(colour= "black"),
			plot.background=element_blank(),
			plot.background= element_blank(),
			legend.direction="horizontal",
			legend.background=element_rect(fill= "NA"),
			legend.key.height= unit(5,"lines"),
			legend.title.align = 0.5,
			panel.border= element_rect(fill= NA),
			legend.key.width = unit(2.5,"lines"),
			plot.margin = unit(rep(0.01, 4), "inches"))


pdf(width = 4, height =3, file = "~/Copy/UrchinAnalyses/Figures/Gaviota_Wavelet_Coherence_MEI.pdf", family = "Times",pointsize = 14)
### plot the partial wavelet coherence		
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
pushViewport(viewport(layout=grid.layout(100,100)))


X <- pwtc.MEI
df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
             s= matrix(X$signif,ncol= 1),
             sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
             phase= matrix(X$phase,ncol= 1))
df1$dx <- cos(df1$phase)
df1$dy <- sin(df1$phase)
df.arrowsb <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.8,]

a1<- ggplot(df1)
a<- a1+geom_tile(aes(x=x,y=log(y,base=2),fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
  geom_contour(aes(x=x,y=log(y,base=2),z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
  #geom_segment(aes(x=x,y=log(y,base=2),xend=(x+dx),yend=log((y+dy),base=2)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsb)+
  scale_y_continuous(trans= "reverse",limits=c(log(108,base=2),1),expand= c(0,0),breaks= log(c(3,6,12,24,36,60,108),base=2),labels= c(3,6,12,24,36,60,108))+
  scale_fill_gradientn(colours= tim.colors(100),limits= c(0,1), name= "time series wavelet coherency")+
  geom_ribbon(aes(x=X$t,ymin=log(X$coi,base=2),ymax= log(108,base=2)),fill= "white",alpha=0.4)+
  scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,296,24), 
     labels= paste0("'",substr(as.character(year(date)[seq(12,296,24)]),start=3,stop=4)))+
  options+ylab("")+theme(legend.position=c(0.5,1.23),legend.key.height= unit(.25,"cm"),plot.margin = unit(c(0.7,0.01,0.01,0.01), "inches"))+guides(fill = guide_colourbar(title.position= "top"))

print(a, vp=vplayout(1:100,1:99), more = T)
# X <- pwtc.NPGO
# 	df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
# 				s= matrix(X$signif,ncol= 1),
# 				sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
# 				phase= matrix(X$phase,ncol= 1))
# 	df1$dx <- cos(df1$phase)
#  	df1$dy <- sin(df1$phase)
#  	df.arrowsb <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.95,]
# 
# b1<- ggplot(df1)
# b<- b1+geom_tile(aes(x=x,y=log(y,base=2),fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
# 		geom_contour(aes(x=x,y=log(y,base=2),z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
# 	#	geom_segment(aes(x=x,y=log(y,base=2),xend=(x+dx),yend=log((y+dy),base=2)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsb)+
# 		scale_y_continuous(trans= "reverse",limits=c(log(108,base=2),1),expand= c(0,0),breaks= log(c(3,6,12,24,36,60,108),base=2),labels= c(3,6,12,24,36,60,108))+
# 		scale_fill_gradientn(colours= tim.colors(100),limits= c(0,1), name= expression(paste(r^2)))+
# 		geom_ribbon(aes(x=X$t,ymin=log(X$coi,base=2),ymax= log(108,base=2)),fill= "white",alpha=0.4)+
# 		scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,296,24), 
# 					labels= paste0("'",substr(as.character(year(date)[seq(12,296,24)]),start=3,stop=4)))+
# 				options+ylab("")
# 
# print(b, vp=vplayout(62:100,3:99), more = T)

grid.text("MEI vs Settlement", y=unit(1, "npc") - unit(.21, "npc"),x=unit(1, "npc") - unit(.85, "npc"),hjust=0, gp=gpar(fontsize=12))
#grid.text("b) NPGO vs Settlement", y=unit(1, "npc") - unit(.58, "npc"),x=unit(1, "npc") - unit(.85, "npc"),hjust=0, gp=gpar(fontsize=12))
grid.text("period (months)",rot= 90, y=unit(1, "npc") - unit(.7, "npc"),x=unit(1, "npc") - unit(.97, "npc"),hjust=0)
dev.off() 



pdf(width = 5, height =7, file = "/Users/Dan/Copy/DissertationBackup/ExitSeminar/Gaviota_Wavelet_Coherence8.20.14.pdf", family = "Times",pointsize = 14)
### plot the partial wavelet coherence		
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
pushViewport(viewport(layout=grid.layout(100,100)))

X <- pwtc.MEI
	df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
				s= matrix(X$signif,ncol= 1),
				sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
				phase= matrix(X$phase,ncol= 1))
	df1$dx <- cos(df1$phase)
 	df1$dy <- sin(df1$phase)
 	df.arrowsb <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.8,]

a1<- ggplot(df1)
a<- a1+geom_tile(aes(x=x,y=y,fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
		geom_contour(aes(x=x,y=y,z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
		geom_segment(aes(x=x,y=y,xend=(x+dx),yend=(y+dy)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsb)+
		scale_y_continuous(trans= "log2", limits = c(2,78),expand= c(0,0),breaks= c(4,8,16,32,64))+
		scale_fill_gradient2(low= "blue",mid= "pink",high="red",midpoint= 0.5,limits= c(0,1), name= expression(paste(r^2)))+
		geom_ribbon(aes(x=X$t,ymin=X$coi,ymax= 78),fill= "white",alpha=0.75)+
		scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), 
					labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+ylab("")+theme(axis.text.x = element_text(colour= NA))
print(a, vp=vplayout(1:25,1:99), more = T)
X <- pwtc.SST
	df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
				s= matrix(X$signif,ncol= 1),
				sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
				phase= matrix(X$phase,ncol= 1))
	df1$dx <- cos(df1$phase)
 	df1$dy <- sin(df1$phase)
 	df.arrowsb <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.8,]


b1<- ggplot(df1)
b<- b1+geom_tile(aes(x=x,y=y,fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
		geom_contour(aes(x=x,y=y,z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
		geom_segment(aes(x=x,y=y,xend=(x+dx),yend=(y+dy)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsb)+
		scale_y_continuous(trans= "log2", limits = c(2,78),expand= c(0,0),breaks= c(4,8,16,32,64))+
		scale_fill_gradient2(low= "blue",mid= "pink",high="red",midpoint= 0.5,limits= c(0,1), name= expression(paste(r^2)))+
		geom_ribbon(aes(x=X$t,ymin=X$coi,ymax= 78),fill= "white",alpha=0.75)+
		scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), 
					labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+ylab("Period (months)                                        ")+theme(axis.text.x = element_text(colour= NA))
print(b, vp=vplayout(25:50,1:99), more = T)
X <- pwtc.UP
	df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
				s= matrix(X$signif,ncol= 1),
				sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
				phase= matrix(X$phase,ncol= 1))
	df1$dx <- cos(df1$phase)
 	df1$dy <- sin(df1$phase)
 	df.arrowsc <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.8,]
c1 <- ggplot(df1)
c <- c1+geom_tile(aes(x=x,y=y,fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
		geom_contour(aes(x=x,y=y,z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
		geom_segment(aes(x=x,y=y,xend=(x+dx),yend=(y+dy)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsc)+
		scale_y_continuous(trans= "log2", limits = c(2,78),expand= c(0,0),breaks= c(4,8,16,32,64))+
		scale_fill_gradient2(low= "blue",mid= "pink",high="red",midpoint= 0.5,limits= c(0,1), name= expression(paste(r^2)))+
		geom_ribbon(aes(x=X$t,ymin=X$coi,ymax= 78),fill= "white",alpha=0.75)+
			scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), 
					labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+ylab("")+theme(axis.text.x = element_text(colour= NA))
print(c, vp=vplayout(50:75,1:99), more = T)
X <- pwtc.chl
	df1 <- cbind(expand.grid(y=X$period,x=X$xaxis),z=matrix(X$rsq,ncol= 1),
				s= matrix(X$signif,ncol= 1),
				sig = as.numeric(matrix(X$rsq,ncol= 1)>matrix(X$signif,ncol= 1)),
				phase= matrix(X$phase,ncol= 1))
	df1$dx <- cos(df1$phase)
 	df1$dy <- sin(df1$phase)
 	df.arrowsc <- df1[df1$y%in%X$period[c(1,3,5,7,seq.int(9,length(X$period),length.out=32))]&df1$x%in%X$xaxis[seq.int(1,length(X$xaxis),length.out=70)]&df1$sig>0.8,]
d1 <- ggplot(df1)
d <- d1+geom_tile(aes(x=x,y=y,fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
		geom_contour(aes(x=x,y=y,z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
		geom_segment(aes(x=x,y=y,xend=(x+dx),yend=(y+dy)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsc)+
		scale_y_continuous(trans= "log2", limits = c(2,78),expand= c(0,0),breaks= c(4,8,16,32,64))+
		scale_fill_gradient2(low= "blue",mid= "pink",high="red",midpoint= 0.5,limits= c(0,1), name= expression(paste(r^2)))+
		geom_ribbon(aes(x=X$t,ymin=X$coi,ymax= 78),fill= "white",alpha=0.75)+
			scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), 
					labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+ylab("")+theme(legend.position=c(0.1,0.5),legend.key.height= unit(0.5,"cm"))
print(d, vp=vplayout(75:100,1:99), more = T)
grid.text("a) MEI", y=unit(1, "npc") - unit(.03, "npc"),x=unit(1, "npc") - unit(.85, "npc"),hjust=0)
grid.text("b) SST", y=unit(1, "npc") - unit(.28, "npc"),x=unit(1, "npc") - unit(.85, "npc"),hjust=0)
grid.text("c) Upwelling ", y=unit(1, "npc") - unit(.53, "npc"),x=unit(1, "npc") - unit(.85, "npc"),hjust=0)
grid.text("d) SS chl", y=unit(1, "npc") - unit(.78, "npc"),x=unit(1, "npc") - unit(.65, "npc"),hjust=0)
dev.off() 


options <- theme(strip.text.x = element_text(size =16),
			axis.text.y = element_text( colour= "black",size =12),
			axis.text.x = element_text(colour= "black",size =12),
			axis.text.y = element_text( colour= "black",size =12),
			axis.title.x = element_text(colour= NA,size= 0),
			panel.grid.minor = element_line(colour = NA),
			panel.grid.major = element_line(colour = NA),
			panel.margin= unit(.5,"lines"),
			strip.background = element_rect(fill= NA,colour= NA,size= 0),
			legend.key = element_rect(colour=NA),
			panel.grid = element_line(colour = NA),
			legend.position = "none",
			panel.background=element_blank(),
			plot.background=element_blank(),
			plot.background= element_blank(),
			legend.direction="vertical",
			legend.background=element_rect(fill= "NA"),
			legend.key.height= unit(5,"lines"),
			legend.title.align = 0.25,
			panel.border= element_rect(fill= NA),
			legend.key.width = unit(0.7,"lines"),
			plot.margin = unit(rep(0.01, 4), "inches"))



pdf(width = 7, height =2.5, file = "/Users/Dan/Copy/DissertationBackup/ExitSeminar/SSTGaviota_Wavelet.pdf", family = "Times",pointsize = 14)

### plot the partial wavelet coherence		


b1<- ggplot(df1)
b1+geom_tile(aes(x=x,y=y,fill= z),show_guide = FALSE,hjust=0,vjust=0,interpolate= T)+
		geom_contour(aes(x=x,y=y,z= sig),alpha= 0.25,colour= "black",linend= "round",weight=0.1,size=0.1)+
		geom_segment(aes(x=x,y=y,xend=(x+dx),yend=(y+dy)),size=0.2,arrow = arrow(length = unit(0.08,"cm")),data= df.arrowsb)+
		scale_y_continuous(trans= "log2", limits = c(2,78),expand= c(0,0),breaks= c(4,8,16,32,64))+
		scale_fill_gradient2(low= "blue",mid= "pink",high="red",midpoint= 0.5,limits= c(0,1), name= expression(paste(r^2)))+
		geom_ribbon(aes(x=X$t,ymin=X$coi,ymax= 78),fill= "white",alpha=0.75)+
		scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), 
					labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+ylab("Period (months)")+scale_x_continuous(limits=c(4,max(X$xaxis)-3),expand= c(0,0),breaks= seq(12,256,24), labels= paste0("'",substr(as.character(year(date)[seq(12,256,24)]),start=3,stop=4)))+
				options+theme(legend.position="right",legend.key.height= unit(0.5,"cm"))
dev.off()



a <- ts(set.gav[-c(1:12),c("month","year","SP.gav")]$SP.gav,frequency=12)
b <- ts(set.gav[-c(1:12),c("month","year","MEIlag")]$MEIlag,frequency=12)
c <- ts(set.gav[-c(1:12),c("month","year","BAKUN")]$BAKUN,frequency=12)

astl <- stl(log(a+1/224),s.window="per")
bstl <- stl(b,s.window="per")
cstl <- stl(c,s.window="per")