#Feature Extraction

#library Import
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(zoo)
library(REdaS)

#import dataset:
trialdata <- read.csv(here("WT031623_26dpf_readyforanalysis.csv"))

########velocity and acceleration##########
featuredf <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>% #frame-by-frame instantaneous vel calculation
  mutate(velocity_mmpersec = distance*25) %>% #instantaneous vel conversion to mmpersec
  mutate(acceleration_mmpsps = velocity_mmpersec - lag(velocity_mmpersec)) %>% #instantaneous acceleration calculation
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>% #set behaviors we annotated for sleap but don't need for analysis to NA
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

#######tail angle#######

tailbeats <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "head", "spine1", "caudal"))
tailbeats <- tailbeats[-c(2)]

tailbeats <- unite(tailbeats, col='x-y', c('x', 'y'), sep='-') #group x and y together for manipulation
tailbeats <- spread(tailbeats, key=nodes, value="x-y") #spread out so that nodes are in their own column
tailbeats <- tailbeats %>% #put it in order
  arrange(track,frame)
#string separations
tailbeats <- separate(tailbeats, col=caudal, into=c('caudal.x', 'caudal.y'), sep='-')
tailbeats <- separate(tailbeats, col=head, into=c('head.x', 'head.y'), sep='-')
tailbeats <- separate(tailbeats, col=nose, into=c('nose.x', 'nose.y'), sep='-')
tailbeats <- separate(tailbeats, col=spine1, into=c('spine1.x', 'spine1.y'), sep='-')
#and for some reason there was NA nonsense and str problems, so:
tailbeats[tailbeats == "NA"] <- NA
tailbeats$caudal.x <- as.numeric(tailbeats$caudal.x)
tailbeats$caudal.y <- as.numeric(tailbeats$caudal.y)
tailbeats$head.x <- as.numeric(tailbeats$head.x)
tailbeats$head.y <- as.numeric(tailbeats$head.y)
tailbeats$nose.x <- as.numeric(tailbeats$nose.x)
tailbeats$nose.y <- as.numeric(tailbeats$nose.y)
tailbeats$spine1.x <- as.numeric(tailbeats$spine1.x)
tailbeats$spine1.y <- as.numeric(tailbeats$spine1.y)

calculatetailangle <- function(ny, hy, nx, hx, cy, sy, cx, sx) {
  slope1 = (ny-hy)/(nx-hx)
  slope2 = (cy-sy)/(cx-sx)
  angle1rad = atan((slope1-slope2)/(1+slope1*slope2))
  angle1 = rad2deg(angle1rad)
  angle1
}

tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(theta = calculatetailangle(nose.y, head.y, nose.x, head.x, caudal.y, spine1.y, caudal.x, spine1.x))

featuredf$theta <- tailbeats$theta

######### Nearest Neighbor ######################

#to calculate nearest neighbor, we'll need the head dataframe once more:
head_df <- trialdata %>%
  filter(nodes == "head")

#all we really need is frame, track and xy position coordinates, thus:
nearestneighbordf <- head_df[c(1,4,5,6)]
#this does some data manipulation to organize the dataframe:
nearestneighbordf$xy <- paste(nearestneighbordf$x, nearestneighbordf$y, sep="-")
nearestneighbordf <- nearestneighbordf[-c(3,4)]
nearestneighbordf <- reshape2::dcast(nearestneighbordf, frame~track)
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

#at some point we can decide if we only need to keep the true NN column, but for now I think we can keep it all
featuredf <- merge(featuredf,
                       NNDF,
                       by=c("frame"),
                       all.x=TRUE)
featuredf <- featuredf %>%
  arrange(track, frame)

write_csv(featuredf_tmp, "featureextractiondataframe.csv")


