#Characterize a chase

#library import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#data import, if necessary:
trialdata <- read.csv(here("WT031623_26dpf_readyforanalysis.csv"))

#filter for chase moments, 2 seconds before and after:

#filter for head:
head_df <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#seq along for each instance of chase, for grouping later:
tmp <- head_df %>%
  group_by(track) %>%
  filter(Behavior == "C") %>%
  mutate(ChaseID = seq_along(Behavior))
#then, merge the labels back in for analysis
annotatedchasedf <- merge(x = head_df, 
                          y = tmp, 
                          by = c("frame", "nodetrackid", "nodes", "track", "x", "y", "Behavior", 
                                 "distance", "velocity_mmpersec", "second"),
                          all.x = TRUE)
annotatedchasedf <- annotatedchasedf %>%
  arrange(track, frame)

#chases are actually happening really close together to each other, so we'll do +/- 1 second (for now)
leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}
annotatedchasedf_chases <- annotatedchasedf[leadlag(annotatedchasedf$Behavior == "C", 25, 25),]
#and then from here, we can fill out the rest

roll <- function(x, n) {
  idx <- which(!is.na(x))
  for (i in idx) x[pmax(i - n, 0):pmin(i + n, length(x))] <- x[i]
  return(x)
}
annotatedchasedf_chases$ChaseID <- roll(annotatedchasedf_chases$ChaseID, 25)

#and mutate a time column that resets for each bout:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track,ChaseID) %>%
  mutate(timereset = seq_along(cumsum(ChaseID)))

write.csv(annotatedchasedf_chases, file="isolatedchases.csv")

#really quick, a df of just blanket categories of whether or not a chase is happening
annotatedchasedf_nonchases <- annotatedchasedf[!leadlag(annotatedchasedf$Behavior == "C", 100, 100),]
annotatedchasedf_nonchases <- annotatedchasedf_nonchases %>%
  mutate(chaseorna = "NoChase")
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  mutate(chaseorna = "Chases")
annotatedchasedf <- rbind(annotatedchasedf_nonchases, annotatedchasedf_chases)
#Right, so. Now we can actually do some visualizations!

annotatedchasedf_chases <- read.csv(here("isolatedchases.csv"))

#better time reset bind-on:
annotatedchasedf_chases$timereset <- annotatedchasedf_chases$timereset-26
annotatedchasedf_chases$secondreset <- annotatedchasedf_chases$timereset/25

annotatedchasedf_chases_select <- annotatedchasedf_chases %>%
  group_by(ChaseID) %>%
  filter(ChaseID %in% c(5, 13, 27, 39, 45, 51, 59, 63, 72, 79))

ggplot(annotatedchasedf_chases_select, aes(x=secondreset, 
                                    y=velocity_mmpersec, 
                                    group=ChaseID))+
  geom_point(size=0.2)+
  geom_path(linewidth=0.3)+
#  coord_cartesian(xlim=c(0, 50))+
  theme_classic()+
  theme(text = element_text(size = 15))+ 
  xlab("Time (Seconds)")+
  ylab("Velocity (mm/sec)")+
  facet_wrap(~track, ncol=1)

#acceleration too, for the poster:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track, ChaseID) %>%
  mutate(acceleration = velocity_mmpersec - lag(velocity_mmpersec))

ggplot(annotatedchasedf_chases_select, aes(x=secondreset, 
                                    y=acceleration, 
                                    group=ChaseID))+
  geom_point(size=0.2)+
  geom_path(linewidth=0.3)+
#  coord_cartesian(xlim=c(0, 50))+
  theme_classic()+
  theme(text = element_text(size = 15))+  
  xlab("Time (Seconds)")+
  ylab("Acceleration (mm/sec/sec)")+
  facet_wrap(~track, ncol=1)

#nearest neighbor stuff >.<
nearestneighbordf <- annotatedchasedf_chases[c(1,4,5,6)]
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

names(nearestneighbordf) <- c("frame", "X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4", "X5", "Y5")
trialdata <- nearestneighbordf
trialdata$X1 <- as.numeric(trialdata$X1)
trialdata$Y1 <- as.numeric(trialdata$Y1)
trialdata$X2 <- as.numeric(trialdata$X2)
trialdata$Y2 <- as.numeric(trialdata$Y2)
trialdata$X3 <- as.numeric(trialdata$X3)
trialdata$Y3 <- as.numeric(trialdata$Y3)
trialdata$X4 <- as.numeric(trialdata$X4)
trialdata$Y4 <- as.numeric(trialdata$Y4)
trialdata$X5 <- as.numeric(trialdata$X5)
trialdata$Y5 <- as.numeric(trialdata$Y5)

distance_function <- function(x, y, n, p) {
  sqrt((n-x)^2+(p-y)^2)
}

#ok, now we build distance df:
F1_F2 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X2, trialdata$Y2)
F1_F3 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X3, trialdata$Y3)
F1_F4 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X4, trialdata$Y4)
F1_F5 <- distance_function(trialdata$X1, trialdata$Y1, trialdata$X5, trialdata$Y5)
F2_F3 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X3, trialdata$Y3)
F2_F4 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X4, trialdata$Y4)
F2_F5 <- distance_function(trialdata$X2, trialdata$Y2, trialdata$X5, trialdata$Y5)
F3_F4 <- distance_function(trialdata$X3, trialdata$Y3, trialdata$X4, trialdata$Y4)
F3_F5 <- distance_function(trialdata$X3, trialdata$Y3, trialdata$X5, trialdata$Y5)
F4_F5 <- distance_function(trialdata$X4, trialdata$Y4, trialdata$X5, trialdata$Y5)

trialdata_coordinatedistance <- bind_cols(F1_F2, F1_F3, F1_F4, F1_F5, F2_F3, F2_F4, F2_F5, F3_F4, F3_F5, F4_F5)
names(trialdata_coordinatedistance) <- c("F0_F1", "F0_F2", "F0_F3", "F0_F4", "F1_F2", "F1_F3", "F1_F4", "F2_F3", "F2_F4", "F3_F4")
rm(F1_F2, F1_F3, F1_F4, F1_F5, F2_F3, F2_F4, F2_F5, F3_F4, F3_F5, F4_F5)

tmp <- annotatedchasedf_chases %>%
  filter(track==0)
trialdata_coordinatedistance <- cbind(tmp$frame, trialdata_coordinatedistance)
names(trialdata_coordinatedistance)[1] <- "frame"
specificnndforfigure <- trialdata_coordinatedistance %>%
  filter(frame %in% c(1131:1181))

#quick melt so I can visualize:
specificnndfforfigure <- reshape2::melt(specificnndforfigure, id="frame")
specificnndfforfigure$second <- specificnndfforfigure$frame/25

ggplot(specificnndfforfigure, aes(x=second, y=variable, fill=value))+
  geom_tile()+
  theme_classic()+
  theme(text=element_text(size=15))+
  ylab("Fish Combination")+
  xlab("Time (Seconds)")+
  labs(fill="Distance\nBetween\nFish (mm)")

#ok. we're gonna do a boxplot of each chase, and facet by individual. Just to give a sense of speed for each one
ggplot(annotatedchasedf_chases, aes(x=as.factor(ChaseID), y=velocity_mmpersec))+
  geom_boxplot()+
  facet_wrap(~track, ncol=1)

#let's try doing all chases compared to everything that is not a chase, per individual...
ggplot(annotatedchasedf, aes(x=chaseorna, y=velocity_mmpersec))+
  geom_boxplot(outlier.shape=NA)+
  coord_cartesian(ylim=c(0, 1000))+
  facet_wrap(~track, ncol=1)

#calculate some stats for this...
chasestats <- annotatedchasedf_chases %>%
  group_by(track, ChaseID) %>%
  summarize(minimum_vel = min(velocity_mmpersec),
            quantile25_vel = quantile(velocity_mmpersec, probs=c(.25)),
            quantile50_vel = quantile(velocity_mmpersec, probs=c(.50)),
            quantile75_vel = quantile(velocity_mmpersec, probs=c(.75)),
            maximum_vel = max(velocity_mmpersec),
            mean_vel = mean(velocity_mmpersec)
  )

chasestats$minimum_vel <- as.numeric(chasestats$minimum_vel)
chasestats$quantile25_vel <- as.numeric(chasestats$minimum_vel)
chasestats$quantile50_vel <- as.numeric(chasestats$quantile50_vel)
chasestats$quantile75_vel <- as.numeric(chasestats$quantile75_vel)
chasestats$maximum_vel <- as.numeric(chasestats$maximum_vel)
chasestats$mean_vel <- as.numeric(chasestats$mean_vel)

#real quick, will try and visualize this:
chasestats_forgraph <- reshape2::melt(chasestats, id=c("track",
                                                       "ChaseID"))

ggplot(chasestats_forgraph, aes(x=as.factor(variable), y=value))+
  geom_point(size=1)+
  theme_classic()+
  scale_y_log10()+
  facet_wrap(~track, ncol=1)

write.csv(chasestats, file="chasestats.csv")