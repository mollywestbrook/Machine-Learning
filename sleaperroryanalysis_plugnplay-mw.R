## Plug n Play Error Analysis: Single Model Work

#libraries
library(rhdf5)
library(here)
library(tidyverse)
library(zoo)
library(reshape2)
library(data.table)
library(reticulate)
library(patchwork)

#file import
#change the file name in quotes to your current .h5 file
modeldata <- "model_comparison_6-13-2023.003_HiQualCamTestVid2_25fps_107.analysis.h5"
figuretitle <- "input-everything-before-'analysis.h5'"

modeldata_bits <- h5ls(here(modeldata))

#overall model histograms

#instance scores
modeldata_instancescores <- as.data.frame(h5read(here(modeldata), "/instance_scores"))
modeldata_instancescores <- modeldata_instancescores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_instancescores <- as.matrix(modeldata_instancescores)
hist(modeldata_instancescores, breaks=50, xlab="Instance Scores", ylab="Frequency", 
     main=paste("Instance Scores For", figuretitle, sep=" "))

#point scores
modeldata_pointscores <- as.data.frame(h5read(here(modeldata), "/point_scores"))
modeldata_pointscores <- modeldata_pointscores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_pointscores <- as.matrix(modeldata_pointscores)
hist(modeldata_pointscores, breaks=50, xlab="Point Scores", ylab="Frequency", 
     main=paste("Point Scores For", figuretitle, sep=" "))

#tracking scores
modeldata_trackingscores <- as.data.frame(h5read(here(modeldata), "/tracking_scores"))
modeldata_trackingscores <- modeldata_trackingscores %>% mutate_all(~ifelse(is.nan(.), NA, .))
modeldata_trackingscores <- as.matrix(modeldata_trackingscores)
hist(modeldata_trackingscores, breaks=50, xlab="Tracking Scores", ylab="Frequency", 
     main=paste("Tracking Scores For", figuretitle, sep=" "))

####

#node-by-node point scores
#organize the data first
modeldata_nodes <- t(as.data.frame(h5read(here(modeldata), "/node_names"))) #transpose node names frame
modeldata_nodes <- rep(modeldata_nodes, modeldata_numtracks) #repeat the number of nodes for each number of tracks
modeldata_pointscores <- as.data.frame(modeldata_pointscores) #turn the pointscores back into a data frame
names(modeldata_pointscores) <- c(modeldata_nodes) #change the name of the pointscores to the node names
modeldata_pointscores$frame <- seq_along(modeldata_pointscores[,1]) #create a frame column
modeldata_pointscores_melted <- reshape2::melt(modeldata_pointscores, id=c("frame")) #melt to the frame column so data is in wide format and can be filtered

#nose
modeldata_pointscores_noses <- modeldata_pointscores_melted %>%
  filter(variable == "nose")
hist(as.numeric(modeldata_pointscores_noses$value), breaks=50, xlab="Nose Scores", ylab="Frequency", 
     main=paste("Nose Scores For", figuretitle, sep=" "))

#eyeL
modeldata_pointscores_eyeL <- modeldata_pointscores_melted %>%
  filter(variable == "eyeL")
hist(as.numeric(modeldata_pointscores_eyeL$value), breaks=50, xlab="EyeL Scores", ylab="Frequency", 
     main=paste("EyeL Scores For", figuretitle, sep=" "))

#eyeR
modeldata_pointscores_eyeR <- modeldata_pointscores_melted %>%
  filter(variable == "eyeR")
hist(as.numeric(modeldata_pointscores_eyeR$value), breaks=50, xlab="EyeR Scores", ylab="Frequency", 
     main=paste("EyeR Scores For", figuretitle, sep=" "))

#head
modeldata_pointscores_head <- modeldata_pointscores_melted %>%
  filter(variable == "head")
hist(as.numeric(modeldata_pointscores_head$value), breaks=50, xlab="Head Scores", ylab="Frequency", 
     main=paste("Head Scores For", figuretitle, sep=" "))

#spine1
modeldata_pointscores_spine1 <- modeldata_pointscores_melted %>%
  filter(variable == "spine1")
hist(as.numeric(modeldata_pointscores_spine1$value), breaks=50, xlab="Spine1 Scores", ylab="Frequency", 
     main=paste("Spine1 Scores For", figuretitle, sep=" "))

#spine2
modeldata_pointscores_spine2 <- modeldata_pointscores_melted %>%
  filter(variable == "spine2")
hist(as.numeric(modeldata_pointscores_spine2$value), breaks=50, xlab="Spine2 Scores", ylab="Frequency", 
     main=paste("Spine2 Scores For", figuretitle, sep=" "))

#caudal

#tail

     
