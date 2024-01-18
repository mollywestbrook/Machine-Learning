### SLEAP Organization

library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)

#set file name to call in via rdhdf5 package
filename <- "predictions_MW_05-29-2023_group5_R.mp4.000_05-29-2023_group5_R.analysis.h5"

trialdata_bits <- h5ls(here(filename))

trialdata_tracks <- as.data.frame(h5read(here(filename), "/tracks"))
# nrow(trialdata_tracks)
# ncol(trialdata_tracks)
# 
# #calls for accessing each of the elements of the h5
# trialdata_edge_inds <- as.data.frame(h5read(here(filename), "/edge_inds"))
# trialdata_edge_names <- as.data.frame(h5read(here(filename), "/edge_names"))
# trialdata_instance_scores <- as.data.frame(h5read(here(filename), "/instance_scores"))
# trialdata_labels_path <- as.data.frame(h5read(here(filename), "/labels_path"))
# trialdata_node_names <- as.data.frame(h5read(here(filename), "/node_names"))
# 
# #calculate a mean for point scores
# trialdata_point_scores <- as.data.frame(h5read(here(filename), "/point_scores"))
# trialdata_point_scores <- trialdata_point_scores %>% mutate_all(~ifelse(is.nan(.), NA, .))
# pointscore_mean <- as.data.frame(sapply(trialdata_point_scores, mean, na.rm=TRUE))
# 
# trialdata_track_names <- as.data.frame(h5read(here(filename), "/track_names"))
# 
# trialdata_track_occupancy <- as.data.frame(h5read(here(filename), "/track_occupancy"))
# trialdata_track_occupancy <- as.data.frame(t(trialdata_track_occupancy)) #transpose the track occupancy df
# unique(trialdata_track_occupancy)
# 
# #calculate mean of tracking scores
# trialdata_tracking_scores <- as.data.frame(h5read(here(filename), "/tracking_scores"))
# trialdata_tracking_scores <- trialdata_tracking_scores %>% mutate_all(~ifelse(is.nan(.), NA, .))
# trackscore_mean <- as.data.frame(sapply(trialdata_tracking_scores, mean, na.rm=TRUE))

#########################################################################################

#the file to be working with is trialdata_tracks

#step one: interpolate missing values + smooth over 10-frame window

trialdata_tracks <- trialdata_tracks %>% mutate_all(~ifelse(is.nan(.), NA, .)) #replace nans with na
tracks_interp <- sapply(trialdata_tracks, na.approx, na.rm=FALSE) #interpolate
tracks_interp <- as.data.frame(tracks_interp)
tracks_smooth <- sapply(tracks_interp, function(x) rollmean(x, 10, fill=NA)) #smooth
tracks_smooth <- as.data.frame(tracks_smooth)

rm(trialdata_tracks, tracks_interp) #drop tmp files we don't need

#step 2: grab track names from other elements of the h5 file, to keep this universal to any h5 file input

#organize track names
trialdata_track_names <- as.data.frame(h5read(here(filename), "/track_names"))
names(trialdata_track_names) <- "tracks"
trialdata_track_names$tracks <- as.numeric(str_replace(trialdata_track_names$tracks, "track_", ""))

#organize node names
trialdata_node_names <- as.data.frame(h5read(here(filename), "/node_names"))
names(trialdata_node_names) <- "nodes"

#this chunk assembles an info dataframe of our labels
track_names <- rep(trialdata_node_names$nodes, count(trialdata_track_names)*2)
track_names <- reshape2::melt(track_names)
names(track_names) <- "nodes"
tracks <- rep(trialdata_track_names, count(trialdata_node_names)*2)
tracks <- reshape2::melt(tracks)
tracks$value <- sort(tracks$value)
names(tracks) <- c("tracks, label")
track_labels <- as.data.frame(cbind(track_names$nodes, tracks$tracks))
names(track_labels) <- c("nodes", "track")

#this chunk gets our xy coordinate labels
coordinate_mini <- c(rep("x", count(trialdata_node_names)), rep("y", count(trialdata_node_names)))
coordinate <- rep(coordinate_mini, count(trialdata_track_names))
track_labels <- cbind(track_labels, coordinate)

#this chunk removes tmp files we don't need anymore
rm(track_names, tracks, trialdata_node_names, trialdata_track_names)

#this takes our info dataframe and our tracks dataframe and smooshes them together
rep_track_labels <- as.data.frame(lapply(track_labels, rep, each=nrow(tracks_smooth)))
tracks_smooth_tmp <- reshape2::melt(tracks_smooth, value.name="data")

#this sets up our new dataframe to get ready to be cast into the x/y coordinates
tracks_final <- cbind(rep_track_labels, tracks_smooth_tmp)
tracks_final$track <- as.numeric(tracks_final$track)
tracks_final <- tracks_final %>%
  arrange(track, nodes)
tracks_final$nodetrackid <- cumsum(!duplicated(tracks_final[c(1,2)]))
tracks_final <- tracks_final[-c(4)]

tracks_final <- tracks_final %>%
  group_by(nodetrackid, coordinate) %>%
  mutate(frame = row_number())

#the cast into the new coordinates
tracks_final <- reshape2::dcast(tracks_final, nodetrackid + nodes + track + frame  ~ coordinate, value.var = "data")
rm(tracks_smooth, tracks_smooth_tmp, rep_track_labels, track_labels)

#these two chunks are for molly's own set up, but can adjust as needed for yours
#and we need to make sure to convert these to mm:
# tracks_final$x <- tracks_final$x/5.79
# tracks_final$y <- tracks_final$y/5.79
#add a second:
# tracks_final$second <- tracks_final$frame/100

fwrite(tracks_final, file = "WT031623_26dpf_filetofocuson.csv", row.names=FALSE)

#and tracks_final will be our dataset with everything, all the things, so beautiful

######################################################################################

#skip here if already organized

tracks_final <- read_csv(here("pilotdataset_organized_tracksonly_raw.csv"))

#basic xyplots as a visual to ensure everything is looking ok:

xypositionplot <- function(n) {
  ggplot()+
    geom_point(n, mapping=aes(x=x, y=y, group=nodes, color=nodes), size=0.7, alpha=0.5)+
    theme_classic()+
    labs(title=paste("fish", n$track, sep="_"))+
    xlab("Fish X Position (mm)")+
    ylab("Fish Y Position (mm)")+
    coord_cartesian(xlim=c(0, 224.5), ylim=c(0, 168.39))+
    theme_classic()+
    theme(text = element_text(size = 15))   
}

tracklist <- split(tracks_final, tracks_final$track)
plotlist <- lapply(tracklist, xypositionplot)
plotnames <- paste("fish", unique(tracks_final$track), sep="_")
names(plotlist) <- plotnames

lapply(names(plotlist), 
       function(x) ggsave(filename=paste(x,"_xyposition",".png",sep=""), plot=plotlist[[x]], width = 7, height = 4, units = "in", path=here("Figures")))