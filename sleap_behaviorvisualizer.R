## Chase Analysis: Individual Variation

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)
library(ggpubr)

#based on bfd generated in 'behaviorextraction'

#tmp check with a single chase, to make sure everything works:

behaviorlist <- split(bdf, bdf$BehaviorID)

#Big ol figure generator:

behaviorlist <- split(bdf, bdf$BehaviorID)

arrangeplots <- function (n){
  p1 = ggplot(n, aes(x=second, y=velocity_mmpersec, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p2 = ggplot(n, aes(x=second, y=acceleration_mmpsps, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p3 = ggplot(n, aes(x=second, y=theta, group=BehaviorID))+
    geom_line()+
    theme_classic()+
    facet_wrap(~track, ncol=1)
  
  p4 = ggplot(n, aes(x=second, y=truenn, fill=truenndist))+
    geom_tile()+
    theme_classic()+
    ylab("Fish Combination")+
    xlab("Time (Seconds)")+
    labs(fill="Dist Btwn\nFish (mm)")+
    theme(text=element_text(size=15),
          legend.title = element_text(size=10),
          legend.text = element_text(size=8),
          legend.position = c(0.95, 0.9),
          legend.key.height = unit(0.25, 'cm'),
          legend.background = element_rect(fill="transparent"))
  
  p5 = ggplot(n, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
    geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25)+
    theme_classic()+
    ylab("Vector Y Component")+
    xlab("Vector X Component")+
    labs(color="Individual")+
    theme(text=element_text(size=10),
          legend.position = c(0.95, 0.9),
          legend.key.height = unit(0.25, 'cm'),
          legend.background = element_rect(fill="transparent"))
  
  p6 = ggplot(n)+
    geom_histogram(mapping=aes(x=angletofish0, fill="Angle to fish 0"), binwidth=5)+
    geom_histogram(mapping=aes(x=angletofish1, fill="Angle to fish 1"), binwidth=5)+
    geom_histogram(mapping=aes(x=angletofish2, fill="Angle to fish 2"), binwidth=5)+
    geom_histogram(mapping=aes(x=angletofish3, fill="Angle to fish 3"), binwidth=5)+
    geom_histogram(mapping=aes(x=angletofish4, fill="Angle to fish 4"), binwidth=5)+
    coord_polar()+
    theme_bw()+
    theme(legend.position = c(0.8, 0.28),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x = element_text(size=8))+
    facet_wrap(~track)
  
  finalplot = ggarrange(
    ggarrange(p1, p2, p3, p4, ncol=4, labels = c("Velocity (MMPS)", "Acceleration (MMPSPS)", "Tail Angle (Deg)", "Nearest Neighbor Dist (MM)")),
    ggarrange(p5, p6, ncol=2, labels =c("Vector Field", "Relative Heading (Deg)")),
    nrow=2)
}

arrangedplotlist <- lapply(behaviorlist, arrangeplots)
plotnames <- paste("chase", unique(bdf$BehaviorID), sep='-')
names(arrangedplotlist) <- plotnames
lapply(names(arrangedplotlist), 
       function(x) ggsave(filename=paste(x,".png",sep=""),
                          plot=arrangedplotlist[[x]], 
                          width = 50,
                          height = 60,
                          unit = "cm",
                          path=here("Figures")))


###################################################################################

#tmp invidual checks

#idk if I'll use it but this is a really nice palette, for future:
# 
# palette <- c("#729ECE", "#6FB899", "#31A1B3", "#CCB22B","#FF9E4A")
#   
# #still need to get rid of the 0 bin but other than that this is workable.
# ggplot(chase10)+
#   geom_histogram(mapping=aes(x=angletofish0, fill="Angle to fish 0"), binwidth=5)+
#   geom_histogram(mapping=aes(x=angletofish1, fill="Angle to fish 1"), binwidth=5)+
#   geom_histogram(mapping=aes(x=angletofish2, fill="Angle to fish 2"), binwidth=5)+
#   geom_histogram(mapping=aes(x=angletofish3, fill="Angle to fish 3"), binwidth=5)+
#   geom_histogram(mapping=aes(x=angletofish4, fill="Angle to fish 4"), binwidth=5)+
#   coord_polar()+
#   theme_bw()+
#   theme(legend.position = c(0.8, 0.28),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         axis.text.x = element_text(size=8))+
#   facet_wrap(~track)
# 
# #vector field:
# ggplot(chase10, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
#   geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), linewidth=0.25)+
#   theme_classic()+
#   ylab("Vector Y Component")+
#   xlab("Vector X Component")+
#   labs(color="Individual")+
#   theme(text=element_text(size=10),
#         legend.position = c(0.95, 0.9),
#         legend.key.height = unit(0.25, 'cm'),
#         legend.background = element_rect(fill="transparent"))
# 
# #nearest neighbor:
# ggplot(chase10, aes(x=second, y=truenn, fill=truenndist))+
#   geom_tile()+
#   theme_classic()+
#   ylab("Fish Combination")+
#   xlab("Time (Seconds)")+
#   labs(fill="Dist Btwn\nFish (mm)")+
#   theme(text=element_text(size=15),
#         legend.title = element_text(size=10),
#         legend.text = element_text(size=8),
#         legend.position = c(0.95, 0.9),
#         legend.key.height = unit(0.25, 'cm'),
#         legend.background = element_rect(fill="transparent"))