######################################################
###  Script to generate settlement plots           ###
###  Author:  D.K. Okamoto                         ###
######################################################


### run analytical summary script
source("3_Analysis_Code/Analysis_Summary.R")

### load plotting options 
source("3_Analysis_Code/Ancillary_Funs.R")

### generate plots
fran_plot <- ggplot(aes(as.Date(date), mean_SF),data=subset(set_summary,!(SITE=="Anacapa[SB]"&YEAR==2009)))+
  geom_linerange(aes(ymin= SF_L, ymax= SF_U),colour= "pink",size=0.25)+
  geom_path(aes(y= Est_SF),colour= "red")+
  geom_point(size= 0.1,alpha= 0.5,shape= 21,fill= alpha("grey20",0.5))+
  theme_bw()+
  scale_x_date(date_breaks = "5 years", date_labels = "%y",limits= c(as.Date("1990-01-01"),as.Date(parse_date_time("2016-12-31","%Y-%m-%d"))),expand= c(0,0))+
  facet_grid(SITE~.,scales= 'free_y')+
  ylab(expression(paste(italic("M. franciscanus"),  " ",brush^-1," ",day^-1)))+
  options;fran_plot

# pdf(width = 7, height =7, file = "4_Figures/SFran_settlers.pdf", family = "Times",pointsize = 14)
# fran_plot
# dev.off()  

purp_plot<- ggplot(aes(as.Date(date), mean_SP),data=set_summary)+
  geom_linerange(aes(ymin= SP_L, ymax= SP_U),colour= "pink",size=0.25)+
  geom_path(aes(y= Est_SP),colour= "red")+
  geom_point(size= 0.1,alpha= 0.5,shape= 21,fill= alpha("grey20",0.5))+
  theme_bw()+
  scale_x_date(date_breaks = "5 years", date_labels = "%y",limits= c(as.Date("1990-01-01"),as.Date(parse_date_time("2016-12-31","%Y-%m-%d"))),expand= c(0,0))+
  facet_grid(SITE~.,scales= 'free_y')+
  ylab(expression(paste(italic("S.purpuratus"),  " ",brush^-1," ",day^-1)))+
  options;purp_plot
# 
# pdf(width = 7, height =7, file = "4_Figures/Spurp_settlers.pdf", family = "Times",pointsize = 14)
# purp_plot
# dev.off()

log_frans_plot <- ggplot(aes(as.Date(date), mean_SF+1/14),data=subset(set_summary,!(SITE=="Anacapa[SB]"&YEAR==2009)))+
  geom_linerange(aes(ymin= SF_L+1/14, ymax= SF_U+1/14),colour= "pink",size=0.25)+
  geom_path(aes(y= Est_SF+1/14),colour= "red")+
  geom_point(size= 0.1,alpha= 0.5,shape= 21,fill= alpha("grey20",0.5))+
  theme_bw()+
  scale_x_date(date_breaks = "5 years", date_labels = "%y",limits= c(as.Date("1990-01-01"),as.Date(parse_date_time("2016-12-31","%Y-%m-%d"))),expand= c(0,0))+
  facet_grid(SITE~.,scales= 'free_y')+
  scale_y_continuous(trans= "log10")+
  annotation_logticks(side= "l",base= 10)+
  ylab(expression(paste(italic("M. franciscanus"),  " ",brush^-1," ",day^-1,"+1/14")))+
  options;log_frans_plot

# pdf(width = 7, height =7, file = "4_Figures/SFran_settlers_log.pdf", family = "Times",pointsize = 14)
# log_frans_plot
# dev.off()  

log_purps_plot<- ggplot(aes(as.Date(date), mean_SP+1/14),data=set_summary)+
  geom_linerange(aes(ymin= SP_L+1/14, ymax= SP_U+1/14),colour= "pink",size=0.25)+
  geom_path(aes(y= Est_SP+1/14),colour= "red")+
  geom_point(size= 0.1,alpha= 0.5,shape= 21,fill= alpha("grey20",0.5))+
  theme_bw()+
  scale_y_continuous(trans= "log10")+
  annotation_logticks(side= "l",base= 10)+
  scale_x_date(date_breaks = "5 years", date_labels = "%y",limits= c(as.Date("1990-01-01"),as.Date(parse_date_time("2016-12-31","%Y-%m-%d"))),expand= c(0,0))+
  facet_grid(SITE~.,scales= 'free_y')+
  ylab(expression(paste(italic("S. purpuratus"),  " ",brush^-1," ",day^-1,"+1/14")))+
  options;log_purps_plot

# pdf(width = 7, height =7, file = "4_Figures/Spurp_settlers_log.pdf", family = "Times",pointsize = 14)
# log_purps_plot
# dev.off()

Season_patterns <- ggplot(aes(biweek2,Exp_seas_SP+1/14),data= seas_mod)+
  geom_ribbon(aes(ymax= SP_U+1/14, ymin= SP_L+1/14), fill= "grey80")+
  geom_line()+
  facet_grid(SITE~.,scales= "free_y")+
  geom_point(aes(biweek,mean_SP+1/14),data= set.ag2,size= 1, shape= 21, fill= "grey",alpha= 0.5)+
  ylab(expression(paste("monthly settlement (# ",brush^-1," ",day^-1,")+1/14")))+
  xlab("")+
  # geom_text(aes(x= 10,y= 5,label = phi),parse= TRUE,data= pacf_seas,size =4)+
  scale_y_continuous(trans= "log10")+
  annotation_logticks(side= "l")+
  theme_bw()+
  options+theme(axis.text.x = element_text(size=12))+
  scale_x_continuous(breaks= c(2,4,6,8,10,12)*2,
                     labels= c("Dec","Feb","Apr","Jun","Aug","Oct"),expand= c(0,0));Season_patterns

# pdf(width=3.5,height= 8.5,file= "~/Copy/UrchinAnalyses/Figures/Seas_patterns.pdf",family = "serif",pointsize = 16)
# Season_patterns
# dev.off()


