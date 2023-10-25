#SLEAP Data Basic visualizations + calculations

#library import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#this script is designed for use with dataframes formatted by sleaporganization.r
#this dataframe has information about behaviors annotated by boris data

#import dataset:
trialdata <- read.csv(here("WT031623_26dpf_readyforanalysis.csv"))

#############################################################################################################

#Script To-dos

#Easy stuff
  #add a heat map, both for all animals and for each individual's performance
  #make it so that for the 5 individual figures (position, velocity, acceleration, heatmaps) get saved in one big file together for easy visualization

#Nearest Neighbor
  #calculate Nearest Neighbor among all individuals
    #Molly has a script that does this, but I don't like it and want to work on it
  #filter for moments when animals are close (<2 body lengths, possibly <1 bodylength) away from each other

#"Spike Filter"
  #at Juan's suggestion, the other filter we can use is when velocity and acceleration are high
  #I want to do a couple of statistics on chases before I go in this direction
  #as that will help us identify what the 'average chase' looks like

#Incorporate group cohesion
  #this is an old script I wrote that I'd like to use a bit more of

#Tail beat calc
  #this is mostly just checking over the calculation I wrote

#############################################################################################################
###Figures for each individual fish

##individual position across entire assay:
individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")

triallist <- split(dataframe_full, dataframe_full$track)

plotposition <- function(n) {
  ggplot(data=n, aes(x=x, y=y, color=as.factor(nodes)))+
    geom_point(alpha=0.4, size=0.5)+
    geom_path(alpha=0.4, linewidth=0.3)+
    #title(paste('Individual', n$track, sep="_"))+
    xlab("X Direction (mm)")+
    ylab("Y Direction (mm)")+
    labs(color="Skeletal Nodes")+
    theme_classic()+
    coord_cartesian(x=c(0, 1250), y=c(0, 1000))+
    theme(text = element_text(size = 15))  
}

#save to designated folder
positionslist <- lapply(triallist, plotposition)
names(positionslist) <- figuretitle
lapply(names(positionslist), 
       function(x) ggsave(filename=paste(x,"_positions",".png",sep=""),
                          plot=positionslist[[x]], width = 13, height = 7, units = "in",
                          path=here("HotFigures")))
#Ideally I would like to stitch all five figures together into a panel as well

##Velocity Over Time

#velocity calculation, using the head node as our centroid
head_df <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>% #this assumes the video is filmed at 25fps
  mutate(second = frame/25) %>% #stop here if no boris data
  mutate(Behavior = na_if(Behavior, "S")) %>% #these filter out extra behaviors noted in boris
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#individual velocity over time:
individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")
headlist <- split(head_df, head_df$track)

plotheadvelocity <- function(n) {
  ggplot(n, aes(y=velocity_mmpersec, x=second))+
    geom_point(size=0.5)+
    geom_path(linewidth=0.4)+
    xlab("Time (Seconds")+
    ylab("Velocity (mm/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

#save to designated folder
headvelocitylist <- lapply(headlist, plotheadvelocity)
names(headvelocitylist) <- figuretitle
lapply(names(headvelocitylist), 
       function(x) ggsave(filename=paste(x,"_velocity",".png",sep=""),
                          plot=headvelocitylist[[x]], width = 13, height = 7, units = "in",
                          path=here("HotFigures")))

#and the version where all figures displayed in the same chart
facetvelocity <- ggplot(head_df, aes(y=velocity_mmpersec, x=second))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.4)+
  xlab("Time (Seconds")+
  ylab("Velocity (mm/sec)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)

ggsave(filename="facetvelocity.png" ,plot=facetvelocity, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

#color coded for behavior, individual plots:
plotheadvelocity_colorcoded <- function(n) {
  ggplot(n, aes(y=velocity_mmpersec, x=second))+
    geom_point(mapping=aes(color=Behavior, size=Behavior))+
    geom_path(alpha=0.5, linewidth=0.4)+
    scale_size_manual(values= c(5,5,5,5))+
    xlab("Time (Seconds)")+
    ylab("Velocity (mm/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

headvelocitylist_colorcoded <- lapply(headlist, plotheadvelocity_colorcoded)
names(headvelocitylist_colorcoded) <- figuretitle
lapply(names(headvelocitylist_colorcoded), 
       function(x) ggsave(filename=paste(x,"_velocity_colorcoded",".png",sep=""),
                          plot=headvelocitylist_colorcoded[[x]], width = 13, 
                          height = 7, units = "in",
                          path=here("HotFigures")))

#color coded for behavior, facet
facetvelocity_colors <-   
  ggplot(head_df, aes(y=velocity_mmpersec, x=second))+
  geom_point(mapping=aes(color=Behavior, size=Behavior))+
  geom_path(alpha=0.5, linewidth=0.4)+
  scale_size_manual(values= c(2,2,2,2))+
  xlab("Time (Seconds)")+
  ylab("Velocity (mm/sec)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)

ggsave(filename="facetvelocity_colors.png" ,plot=facetvelocity_colors, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

##Acceleration

#calculate from velocity data
head_df <- head_df %>%
  group_by(track) %>%
  mutate(acceleration = velocity_mmpersec - lag(velocity_mmpersec))

individ <- unique(trialdata$track)
figuretitle <- paste('Individual', individ, sep="_")
headlist <- split(head_df, head_df$track)

plotheadacceleration_colorcoded <- function(n) {
  ggplot(n, aes(y=acceleration, x=second))+
    geom_point(mapping=aes(color=Behavior, size=Behavior))+ #color coded for behavior
    geom_path(alpha=0.5, linewidth=0.4)+
    scale_size_manual(values= c(5,5,5,5))+
    xlab("Time (Seconds)")+
    ylab("Acceleration (mm/sec/sec)")+
    theme_classic()+
    theme(text = element_text(size = 15))  
}

headaccelerationlist_colorcoded <- lapply(headlist, plotheadacceleration_colorcoded)
headaccelerationlist_colorcoded

#save individual figures to designated folder
names(headaccelerationlist_colorcoded) <- figuretitle
lapply(names(headaccelerationlist_colorcoded), 
       function(x) ggsave(filename=paste(x,"_acceleration_colorcoded",".png",sep=""),
                          plot=headaccelerationlist_colorcoded[[x]], width = 13, 
                          height = 7, units = "in",
                          path=here("HotFigures")))

#############################################################################################################
### Average Performance

##average velocity

#without outliers
averagevelocity_nooutliers <- ggplot(head_df, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Fish Identity")+
  ylab("Median Velocity")+
  theme_classic()+
  coord_cartesian(y=c(0, 500))+
  theme(text=element_text(size=15))
averagevelocity_nooutliers
#save to designated folder:
ggsave(filename="averagevelocity_nooutliers.png", plot=averagevelocity_nooutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

#average velocity with outliers
averagevelocity_withoutliers <- ggplot(head_df, aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_boxplot(outlier.size=0.5)+
  xlab("Fish Identity")+
  ylab("Median Velocity")+
  theme_classic()+
  theme(text=element_text(size=15))
averagevelocity_withoutliers
#save to designated folder:
ggsave(filename="averagevelocity_withoutliers.png", plot=averagevelocity_withoutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

#average velocity chases only (needs work)
head_df_chasesonly <- head_df %>%
  filter(Behavior == "C")

averagevelocity_chasesonly <- ggplot(head_df_chasesonly, 
                                     aes(x=track, y=velocity_mmpersec, group=as.factor(track)))+
  geom_boxplot(outlier.size=0.5)+
  xlab("Fish Identity")+
  ylab("Median Velocity")+
  theme_classic()+
  theme(text=element_text(size=15))
averagevelocity_chasesonly
#save to designated folder:
ggsave(filename="averagevelocity_chasesonly.png", plot=averagevelocity_chasesonly, 
       width = 7, height = 5, units = "in", path=here("HotFigures"))

## average acceleration

#acceleration with outliers
averageacceleration <- ggplot(head_df, aes(x=track, y=acceleration, group=as.factor(track)))+
  geom_boxplot(outlier.size=0.5)+
  xlab("Fish Identity")+
  ylab("Median Acceleration (mm/sec/sec)")+
  theme_classic()+
  #coord_cartesian(y=c(0, 500))+
  theme(text=element_text(size=15))
averageacceleration
ggsave(filename="averageacceleration.png", plot=averageacceleration, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="averageacceleration_nooutliers.png", plot=averageacceleration_nooutliers, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

#no outliers
averageacceleration_nooutliers <- ggplot(head_df, aes(x=track, y=acceleration, group=as.factor(track)))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Fish Identity")+
  ylab("Median Acceleration (mm/sec/sec)")+
  theme_classic()+
  coord_cartesian(y=c(-10, 10))+
  theme(text=element_text(size=15))
averageacceleration_nooutliers

#chases only
head_df_chasesonly <- head_df %>%
  filter(Behavior == "C")

averageacceleration_chasesonly <- ggplot(head_df_chasesonly, 
                                         aes(x=track, y=acceleration, group=as.factor(track)))+
  geom_boxplot(outlier.shape=NA)+
  xlab("Fish Identity")+
  ylab("Median Acceleration (mm/sec/sec)")+
  theme_classic()+
  coord_cartesian(y=c(-700, 500))+
  theme(text=element_text(size=15))
averageacceleration_chasesonly
ggsave(filename="averageacceleration_chasesonly.png", plot=averageacceleration_chasesonly, 
       width = 7, height = 5, units = "in", path=here("HotFigures"))

##############################################################################################################
### Further Calculations

#this calculates instantaneous vector components from the head node
head_df_vectors <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(vx = ((x - lag(x))*25)) %>% #again, assuming we film at 25fps
  mutate(vy = ((y - lag(y))*25)) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>% #filtering for behaviors other than chases
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#I also have a calculation working on tail beats, but I think it needs some work


#############################################################################################################

## Nearest Neighbor

#to calculate nearest neighbor, we'll need a separate dataframe:

#all we really need is frame, track and xy position coordinates, thus:
nearestneighbordf <- trialdata[c(1,4,5,6)]
#this does some data manipulation to organize the dataframe:
nearestneighbordf$xy <- paste(nearestneighbordf$x, nearestneighbordf$y, sep="-")
nearestneighbordf <- nearestneighbordf[-c(3,4)]
nearestneighbordf <- reshape2::dcast(nearestneighbordf, c(frame)~track)
names(nearestneighbordf) <- c("frame", "track0xy", "track1xy", "track2xy", "track3xy", "track4xy")
nearestneighbordf[c('track0x', 'track0y')] <- str_split_fixed(nearestneighbordf$track0xy, '-', 2)
nearestneighbordf[c('track1x', 'track1y')] <- str_split_fixed(nearestneighbordf$track1xy, '-', 2)
nearestneighbordf[c('track2x', 'track2y')] <- str_split_fixed(nearestneighbordf$track2xy, '-', 2)
nearestneighbordf[c('track3x', 'track3y')] <- str_split_fixed(nearestneighbordf$track3xy, '-', 2)
nearestneighbordf[c('track4x', 'track4y')] <- str_split_fixed(nearestneighbordf$track4xy, '-', 2)
nearestneighbordf <- nearestneighbordf[-c(2:6)]
names(nearestneighbordf) <- c("frame", "X0", "Y0", "X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4")
nearestneighbordf$X0 <- as.numeric(nearestneighbordf$X0)
nearestneighbordf$Y0 <- as.numeric(nearestneighbordf$Y0)
nearestneighbordf$X1 <- as.numeric(nearestneighbordf$X1)
nearestneighbordf$Y1 <- as.numeric(nearestneighbordf$Y1)
nearestneighbordf$X2 <- as.numeric(nearestneighbordf$X2)
nearestneighbordf$Y2 <- as.numeric(nearestneighbordf$Y2)
nearestneighbordf$X3 <- as.numeric(nearestneighbordf$X3)
nearestneighbordf$Y3 <- as.numeric(nearestneighbordf$Y3)
nearestneighbordf$X4 <- as.numeric(nearestneighbordf$X4)
nearestneighbordf$Y4 <- as.numeric(nearestneighbordf$Y4)

#now our distance function:
distance_function <- function(x, y, n, p) {
  sqrt((n-x)^2+(p-y)^2)
}

#grab the frame vector, as we'll need it for each individual matrix:
frame <- nearestneighbordf$frame

#we'll briefly grab the closest two fish at any given time here
#to bind that in at the very end:
F0_F1 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X1, nearestneighbordf$Y1)
F0_F2 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X2, nearestneighbordf$Y2)
F0_F3 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X3, nearestneighbordf$Y3)
F0_F4 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X4, nearestneighbordf$Y4)
F1_F2 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X2, nearestneighbordf$Y2)
F1_F3 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X3, nearestneighbordf$Y3)
F1_F4 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X4, nearestneighbordf$Y4)
F2_F3 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X3, nearestneighbordf$Y3)
F2_F4 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X4, nearestneighbordf$Y4)
F3_F4 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X4, nearestneighbordf$Y4)
nnearestneighbor_fulldist<- data.frame(frame, F0_F1, F0_F2, F0_F3, F0_F4, F1_F2, F1_F3, F1_F4, F2_F3, F2_F4, F3_F4)
names(nnearestneighbor_fulldist) <- c("frame", "F0_F1", "F0_F2", "F0_F3", "F0_F4", "F1_F2", "F1_F3", "F1_F4", "F2_F3", "F2_F4", "F3_F4")
rm(F0_F1, F0_F2, F0_F3, F0_F4, F1_F2, F1_F3, F1_F4, F2_F3, F2_F4, F3_F4)
nnearestneighbor_fulldist <- reshape2::melt(nnearestneighbor_fulldist, id="frame")
nnearestneighbor_fulldist <- nnearestneighbor_fulldist %>%
  group_by(frame) %>%
  slice(which.min(value))
names(nnearestneighbor_fulldist) <- c("frame", "truenn", "truenndist")


#in order to get actual nearest neighbor, we'll have to be a little funky with our distance matrix.
#We'll calculate a distance matrix for each animal
#then select the animal that is the closest to that particular individual at any particular frame:

#fish one:
F1 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X0, nearestneighbordf$Y0, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish0_NN <- data.frame(frame, F1, F2, F3, F4)
rm(F1, F2, F3, F4)
Fish0_NN <- reshape2::melt(Fish0_NN, id="frame")
Fish0_NN <- Fish0_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish0_NN) <- c("frame", "NNtoF0", "disttoF0")

#fish two
F0 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X0, nearestneighbordf$Y0)
F2 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X1, nearestneighbordf$Y1, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish1_NN <- data.frame(frame, F0, F2, F3, F4)
rm(F0, F2, F3, F4)
Fish1_NN <- reshape2::melt(Fish1_NN, id="frame")
Fish1_NN <- Fish1_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish1_NN) <- c("frame", "NNtoF1", "disttoF1")

#fish three
F0 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X1, nearestneighbordf$Y1)
F3 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X3, nearestneighbordf$Y3)
F4 <- distance_function(nearestneighbordf$X2, nearestneighbordf$Y2, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish2_NN <- data.frame(frame, F0, F1, F3, F4)
rm(F0, F1, F3, F4)
Fish2_NN <- reshape2::melt(Fish2_NN, id="frame")
Fish2_NN <- Fish2_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish2_NN) <- c("frame", "NNtoF2", "disttoF2")

#fish four
F0 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X2, nearestneighbordf$Y2)
F4 <- distance_function(nearestneighbordf$X3, nearestneighbordf$Y3, nearestneighbordf$X4, nearestneighbordf$Y4)
Fish3_NN <- data.frame(frame, F0, F1, F2, F4)
rm(F0, F1, F2, F4)
Fish3_NN <- reshape2::melt(Fish3_NN, id="frame")
Fish3_NN <- Fish3_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish3_NN) <- c("frame", "NNtoF3", "disttoF3")

#fish five
F0 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X0, nearestneighbordf$Y0)
F1 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X1, nearestneighbordf$Y1)
F2 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X2, nearestneighbordf$Y2)
F3 <- distance_function(nearestneighbordf$X4, nearestneighbordf$Y4, nearestneighbordf$X3, nearestneighbordf$Y3)
Fish4_NN <- data.frame(frame, F0, F1, F2, F3)
rm(F0, F1, F2, F3)
Fish4_NN <- reshape2::melt(Fish4_NN, id="frame")
Fish4_NN <- Fish4_NN %>%
  group_by(frame) %>%
  slice(which.min(value))
names(Fish4_NN) <- c("frame", "NNtoF4", "disttoF4")

#final df:
NNDF <- list(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN, nnearestneighbor_fulldist)
NNDF <- NNDF %>% reduce(full_join, by="frame")
#cleanup:
rm(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN, nnearestneighbor_fulldist, nearestneighbordf, frame)
