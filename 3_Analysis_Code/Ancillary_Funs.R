### plotting options
options <- theme(strip.text.y = element_text(size =20,angle= 0),
                 axis.text.y = element_text( colour= "black",size =12),
                 axis.text.x = element_text(colour= "black",size =20),
                 axis.title.y = element_text( colour= "black",size =20),
                 axis.title.x = element_blank(),
                 panel.grid.minor = element_line(colour = NA),
                 panel.grid.major.y = element_line(colour = NA),
                 panel.grid.major.x = element_line(colour = "grey40"),
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


every_nth <- function(x, nth, empty = TRUE, inverse = FALSE) 
{
  if (!inverse) {
    if(empty) {
      x[1:nth == 1] <- ""
      x
    } else {
      x[1:nth != 1]
    }
  } else {
    if(empty) {
      x[1:nth != 1] <- ""
      x
    } else {
      x[1:nth == 1]
    }
  }
}

custom_breaks <- seq(from = 0, to =2,by = 0.05)

