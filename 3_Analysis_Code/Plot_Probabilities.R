
options <- theme(strip.text.y = element_text(size =12,hjust=0.5),
                 axis.text.y = element_text( colour= "black",size =12,hjust=0.5),
                 axis.text.x = element_text(colour= "black",size =12),
                 axis.title.y = element_text( colour= "black",size =12,vjust=1),
                 axis.title.x = element_blank(),
                 panel.grid.minor = element_line(colour = NA),
                 panel.grid.major = element_line(colour = NA),
                 strip.background = element_rect(fill= NA,colour= NA,size= 0),
                 legend.key = element_rect(colour=NA,fill=NA),
                 panel.grid = element_line(colour = NA),
                 legend.position = "none",
                 panel.background=element_rect(fill=NA),
                 plot.background=element_blank(),
                 panel.border = element_rect(colour = "black",fill= NA),
                 legend.direction="vertical",
                 legend.background=element_rect(fill= NA, colour= NA),
                 legend.key.height= unit(1,"lines"),
                 legend.title = element_text(size = 0),
                 legend.text = element_text(size = 12),
                 legend.title.align = 0,
                 legend.key.width = unit(2,"lines"),
                 plot.margin = unit(c(-0.1,0.01,0.01,0.2), "lines"),
                 panel.margin = unit(c(-0.1,0.01,0.01,0.2), "lines"))

year_labs <- paste0("'",substr(as.character(1991:2015),start=3,stop=4))
year_labs[seq(from = 2, to= 24, by= 2)] <- ""

sp_monthly$Region <- factor(with(sp_monthly,ifelse(SITE=="FBPC","Fort Bragg",ifelse(SITE%in%c("OCNBCH","SIO"),"San Diego","SB Channel"))))
levels(sp_monthly$Region)
sp_monthly$Region <- factor(sp_monthly$Region,levels= c("Fort Bragg", "SB Channel","San Diego"))
sp_monthly$date <- parse_date_time(with(sp_monthly,paste(15,MONTH,YEAR, sep= "-")),"%d-%m-%y")
sp_monthly$site <- factor(sp_monthly$SITE)
levels(sp_monthly$site) <- c("Anacapa","Fort Bragg","Gaviota","Ocean Beach","Ellwood","Stearns Wharf","Scripps Pier")
sp_monthly$YEAR_FRAC <- sp_monthly$YEAR+sp_monthly$MONTH/12
a<- ggplot(data=sp_monthly[sp_monthly$Region=="San Diego",])+
  scale_x_continuous(breaks= 1991:2015, labels= year_labs,limits= c(1990,2015.8),expand= c(0,0))+
  scale_y_continuous(trans="log10", limits=c(0.00085,30),expand=c(0,0), breaks= c(0.001,0.01,0.1,1,10,100,500), labels= c("<0.001",0.01,0.1,1,10,100,500))+options+
  theme_bw()+geom_path(aes(YEAR_FRAC,ifelse(exp(Est_SP)<0.001,0.001,exp(Est_SP)),colour= site),alpha= 0.6,size= .7,alpha= 0.8)+
  facet_grid(Region~.)+options+xlab("")+ylab("")+
  scale_colour_manual(values= c("red","black"))+annotation_logticks(base=10,sides= "l")+
  theme(legend.position=c(0.5,0.93),
        legend.direction="horizontal")

b<- ggplot(data=sp_monthly[sp_monthly$Region=="SB Channel",])+
  scale_x_continuous(breaks= 1991:2015, labels= year_labs,limits= c(1990,2015.8),expand= c(0,0))+
  scale_y_continuous(trans="log10", limits=c(0.00085,30),expand=c(0,0), breaks= c(0.001,0.01,0.1,1,10,100,500), labels= c("<0.001",0.01,0.1,1,10,100,500))+options+
  theme_bw()+geom_path(aes(YEAR_FRAC,ifelse(exp(Est_SP)<0.001,0.001,exp(Est_SP)),colour= site),alpha= 0.6,size= .7,alpha= 0.8)+
  facet_grid(Region~.)+options+xlab("")+ylab("")+
  scale_colour_manual(values= c("red","black","blue","yellow"))+annotation_logticks(base=10,sides= "l")+
  theme(axis.text.x = element_text(colour= "black",size =0),
        axis.title.x = element_text(colour= "black",size =0),
        legend.position=c(0.5,0.93),legend.direction="horizontal")

c<- ggplot(data=sp_monthly[sp_monthly$Region=="Fort Bragg",])+
  scale_x_continuous(breaks= 1991:2015, labels= year_labs,limits= c(1990,2015.8),expand= c(0,0))+
  scale_y_continuous(trans="log10", limits=c(0.00085,30),expand=c(0,0), breaks= c(0.001,0.01,0.1,1,10,100,500), labels= c("<0.001",0.01,0.1,1,10,100,500))+options+
  theme_bw()+geom_path(aes(YEAR_FRAC,ifelse(exp(Est_SP)<0.001,0.001,exp(Est_SP)),colour= site),alpha= 0.6,size= .7,alpha= 0.8)+
  facet_grid(Region~.)+options+xlab("")+ylab("")+
  scale_colour_manual(values= c(1,1,"red","grey","purple",1,"grey"))+annotation_logticks(base=10,sides= "l")+
  theme(axis.text.x = element_text(colour= "black",size =0),
        axis.title.x = element_text(colour= "black",size =0),
        plot.margin= unit(c(0.4,0.01,0,0.2), "lines"))

pdf(width=8,height= 4.5,file= "~/Copy/UrchinAnalyses/Figures/monthly_probabilities.pdf",family = "sans",pointsize = 16)
vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
grid.arrange(c,b,a,heights= c(1.1,1,1.1))
grid.text(expression(paste("monthly settlement (# ",brush^-1," ",day^-1,")")),rot= 90, y=unit(1, "npc") - unit(.5, "npc"),x=unit(1, "npc") - unit(.98, "npc"),hjust=0.5)

dev.off()