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

########velocity, acceleration, vector x, vector y, cleanup ####################
featuredf <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(second = frame/25) %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>% #frame-by-frame instantaneous vel calculation
  mutate(velocity_mmpersec = distance*25) %>% #instantaneous vel conversion to mmpersec
  mutate(acceleration_mmpsps = velocity_mmpersec - lag(velocity_mmpersec)) %>% #instantaneous acceleration calculation
  mutate(vx = ((x - lag(x))*25)) %>% #velocity vector x component
  mutate(vy = ((y - lag(y))*25)) %>% #velocity vector y component
  mutate(Behavior = na_if(Behavior, "S")) %>% #set behaviors we annotated for sleap but don't need for analysis to NA
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

####### tail angle #####################################################

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
tailbeats[] <- lapply(tailbeats, as.numeric)

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

######### Nearest Neighbor ############################################

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
nearestneighbordf[] <- lapply(nearestneighbordf, as.numeric)

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

frame <- data.frame(frame) #for frame length preservation

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
names(Fish0_NN) <- c("frame", "relativeNN", "relativedist")
Fish0_NN <- merge(frame, Fish0_NN, by=c("frame"), all.x=TRUE) #preserve frame length
Fish0_NN$track <- 0

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
names(Fish1_NN) <- c("frame", "relativeNN", "relativedist")
Fish1_NN <- merge(frame, Fish1_NN, by=c("frame"), all.x=TRUE)
Fish1_NN$track <- 1

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
names(Fish2_NN) <- c("frame", "relativeNN", "relativedist")
Fish2_NN <- merge(frame, Fish2_NN, by=c("frame"), all.x=TRUE)
Fish2_NN$track <- 2

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
names(Fish3_NN) <- c("frame", "relativeNN", "relativedist")
Fish3_NN <- merge(frame, Fish3_NN, by=c("frame"), all.x=TRUE)
Fish3_NN$track <- 3

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
names(Fish4_NN) <- c("frame", "relativeNN", "relativedist")
Fish4_NN <- merge(frame, Fish4_NN, by=c("frame"), all.x=TRUE)
Fish4_NN$track <- 4

#final df:
NNDF <- rbind(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN)
NNDF <- merge(NNDF, nnearestneighbor_fulldist, by=c("frame"), all.x=T)
NNDF <- NNDF %>%
  arrange(track, frame)

featuredf <- merge(featuredf, NNDF, by=c("frame", "track"), all.x=T) #this line has not been checked

#cleanup:
rm(Fish0_NN, Fish1_NN, Fish2_NN, Fish3_NN, Fish4_NN, nnearestneighbor_fulldist, frame)

########## relative heading angle #############################################

#using the same dataframe as nearest neighbor above, just do a quick rename:
bearingangledf <- nearestneighbordf
frame <- bearingangledf$frame

calculatebearingangle <- function(X1, Y1, X2, Y2) {
  bearing_rad = atan2(X2 - X1, Y1 - Y2)
  bearing = rad2deg(bearing_rad)
  bearing
}

#the calculation for each fish pair with repeats
F0_F0 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X0, bearingangledf$Y0)
F0_F1 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X1, bearingangledf$Y1)
F0_F2 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X2, bearingangledf$Y2)
F0_F3 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X3, bearingangledf$Y3)
F0_F4 <- calculatebearingangle(bearingangledf$X0, bearingangledf$Y0, bearingangledf$X4, bearingangledf$Y4)
F1_F0 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X0, bearingangledf$Y0)
F1_F1 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X1, bearingangledf$Y1)
F1_F2 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X2, bearingangledf$Y2)
F1_F3 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X3, bearingangledf$Y3)
F1_F4 <- calculatebearingangle(bearingangledf$X1, bearingangledf$Y1, bearingangledf$X4, bearingangledf$Y4)
F2_F0 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X0, bearingangledf$Y0)
F2_F1 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X1, bearingangledf$Y1)
F2_F2 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X2, bearingangledf$Y2)
F2_F3 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X3, bearingangledf$Y3)
F2_F4 <- calculatebearingangle(bearingangledf$X2, bearingangledf$Y2, bearingangledf$X4, bearingangledf$Y4)
F3_F0 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X0, bearingangledf$Y0)
F3_F1 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X1, bearingangledf$Y1)
F3_F2 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X2, bearingangledf$Y2)
F3_F3 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X3, bearingangledf$Y3)
F3_F4 <- calculatebearingangle(bearingangledf$X3, bearingangledf$Y3, bearingangledf$X4, bearingangledf$Y4)
F4_F0 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X0, bearingangledf$Y0)
F4_F1 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X1, bearingangledf$Y1)
F4_F2 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X2, bearingangledf$Y2)
F4_F3 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X3, bearingangledf$Y3)
F4_F4 <- calculatebearingangle(bearingangledf$X4, bearingangledf$Y4, bearingangledf$X4, bearingangledf$Y4)

#build the df relative to each fish to line up with our feature extration
f0df <- data.frame(frame,F0_F0, F0_F1, F0_F2, F0_F3, F0_F4)
f1df <- data.frame(frame,F1_F0, F1_F1, F1_F2, F1_F3, F1_F4)
f2df <- data.frame(frame,F2_F0, F2_F1, F2_F2, F2_F3, F2_F4)
f3df <- data.frame(frame,F3_F0, F3_F1, F3_F2, F3_F3, F3_F4)
f4df <- data.frame(frame,F4_F0, F4_F1, F4_F2, F4_F3, F4_F4)
names(f0df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f1df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f2df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f3df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
names(f4df) <- c("frame", "angletofish0", "angletofish1", "angletofish2", "angletofish3", "angletofish4")
bearingangle_featuredf <- rbind(f0df, f1df, f2df, f3df, f4df) #bind dataframe together
bearingangle_featuredf <- bearingangle_featuredf[-c(1)]
rm(F0_F0, F0_F1, F0_F2, F0_F3, F0_F4,
   F1_F0, F1_F1, F1_F2, F1_F3, F1_F4,
   F2_F0, F2_F1, F2_F2, F2_F3, F2_F4,
   F3_F0, F3_F1, F3_F2, F3_F3, F3_F4,
   F4_F0, F4_F1, F4_F2, F4_F3, F4_F4,
   f0df, f1df, f2df, f3df, f4df)

featuredf <- cbind(featuredf, bearingangle_featuredf)

######## body length ############################################

#this may actually be useful for a potential frame alerter
bodylength <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "caudal")) %>%
  select(frame, nodes, track, x, y) %>%
  unite(col='x-y', c('x', 'y'), sep='-') %>%
  spread(key=nodes, value="x-y") %>%
  separate(col=nose, into=c('nose.x', 'nose.y'), sep='-') %>%
  separate(col=caudal, into=c('caudal.x', 'caudal.y'), sep='-') %>%
  arrange(track, frame) %>%
  mutate_at(c("caudal.x", "caudal.y", "nose.x", "nose.y"), as.numeric) %>%
  mutate(bodylength = distance_function(nose.x, nose.y, caudal.x, caudal.y))


write_csv(featuredf_tmp, "featureextractiondataframe.csv")


