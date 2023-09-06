# Draft for sleap data pipeline

#####################################################################################
#Characterize a chase

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

#this will select +/- two seconds before and after each chase
#sec increased to +/- 4 seconds, as we don't really know how long chases tend to last on average
leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}
annotatedchasedf_chases <- annotatedchasedf[leadlag(annotatedchasedf$Behavior == "C", 100, 100),]
#and then from here, we can fill out the rest
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  fill(ChaseID, .direction=c("downup"))
#and mutate a time column that resets for each bout:
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  group_by(track,ChaseID) %>%
  mutate(timereset = seq_along(cumsum(ChaseID)))

#really quick, a df of just blanket categories of whether or not a chase is happening
annotatedchasedf_nonchases <- annotatedchasedf[!leadlag(annotatedchasedf$Behavior == "C", 100, 100),]
annotatedchasedf_nonchases <- annotatedchasedf_nonchases %>%
  mutate(chaseorna = "NoChase")
annotatedchasedf_chases <- annotatedchasedf_chases %>%
  mutate(chaseorna = "Chases")
annotatedchasedf <- rbind(annotatedchasedf_nonchases, annotatedchasedf_chases)
#Right, so. Now we can actually do some visualizations!

ggplot(annotatedchasedf_chases, aes(x=timereset, y=velocity_mmpersec))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.3)+
  facet_wrap(~track, ncol=1)

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
write.csv(chasestats, file="chasestats.csv")

####################################################################################

#Vector field of chases

#let's use a clean head_df so we can calculate the vector components
#and we'll also clean up the behaviors as well
head_df_vectors <- trialdata %>%
  filter(nodes == "head") %>%
  mutate(vx = ((x - lag(x))*25)) %>%
  mutate(vy = ((y - lag(y))*25)) %>%
  mutate(second = frame/25) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

leadlag <- function(lgl, bef = 1, aft = 1) {
  n <- length(lgl)
  bef <- min(n, max(0, bef))
  aft <- min(n, max(0, aft))
  befx <- if (bef > 0) sapply(seq_len(bef), function(b) c(tail(lgl, n = -b), rep(FALSE, b)))
  aftx <- if (aft > 0) sapply(seq_len(aft), function(a) c(rep(FALSE, a), head(lgl, n = -a)))
  rowSums(cbind(befx, lgl, aftx), na.rm = TRUE) > 0
}

chase_vectors <- head_df_vectors[leadlag(head_df_vectors$Behavior == "C", 10, 10),]

#chase 1:
chase1 <- chase_vectors %>%
  group_by(track) %>%
  slice(1:21)
chase1plot <- ggplot(chase1, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase2 <- chase_vectors %>%
  group_by(track) %>%
  slice(22:42)
chase2plot <-  ggplot(chase2, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase3 <- chase_vectors %>%
  group_by(track) %>%
  slice(43:63) %>%
  mutate(distance = sqrt((x - lag(x))^2) + ((y - lag(y))^2)) %>%
  mutate(velocity_mmpersec = distance*25) %>%
  mutate(acceleration = velocity_mmpersec - lag(velocity_mmpersec))

chase3plot <- ggplot(chase3, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=20))
#to go with chase three I'd like to do a facet of acceleration + velocity with our fish of interest:
chase3_velocity <- ggplot(chase3, aes(x=second, y=velocity_mmpersec))+
  geom_point()+
  geom_path()+
  theme_classic()+
  ylab("Velocity (mm/sec)")+
  xlab("Time (Seconds)")+
  theme(text=element_text(size=18))+
  facet_wrap(~track, ncol=1)
ggsave("chase3_velocity.png", plot=chase3_velocity, width=4, height=8, unit="in")
chase3_acceleration <- ggplot(chase3, aes(x=second, y=acceleration))+
  geom_point()+
  geom_path()+
  theme_classic()+
  ylab("Acceleration (mm/sec/sec)")+
  xlab("Time (Seconds)")+
  theme(text=element_text(size=18))+
  facet_wrap(~track, ncol=1)
ggsave("chase3_acceleration.png", plot=chase3_acceleration, width=4, height=8, unit="in")


chase4 <- chase_vectors %>%
  group_by(track) %>%
  slice(64:84)
chase4plot <- ggplot(chase4, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

chase5 <- chase_vectors %>%
  group_by(track) %>%
  slice(85:105)
chase5plot <- ggplot(chase5, aes(x=x, y=y, group=as.factor(track), color=as.factor(track)))+
  geom_segment(aes(xend=x+vx, yend=y+vy), arrow=arrow(length=unit(0.1, "cm")), size=0.25)+
  theme_classic()+
  ylab("Vector Y Component")+
  xlab("Vector X Component")+
  labs(color="Individual")+
  theme(text=element_text(size=15))

ggsave(filename="chase1.png", plot=chase1plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase2.png", plot=chase2plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase3.png", plot=chase3plot, 
       width = 7, height = 5, units = "in", path=here("HotFigures"))
ggsave(filename="chase4.png", plot=chase4plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))
ggsave(filename="chase5.png", plot=chase5plot, 
       width = 13, height = 7, units = "in", path=here("HotFigures"))

####################################################################################
#Make the nearest neighbor calculation more flexible!

#once again we want to focus on the head df. 
#let's assemble just the essential components of this df, as most of this we don't need

nndf <- head_df[-c(2,3,7:10)]
nndf$xy <- paste(nndf$x, nndf$y, sep='-')
nndf$track <- paste("track", nndf$track, sep='')
nndf <- nndf[-c(3,4)]

#################################################################################

#Tail Beat Dynamics:
#in order to model the tail beat, we'll need these:
#vector between nose and head (this will be our 'base' to measure against)
#and the vector between spine and caudal (this will be the thing that rotates)

#therefore, we'll arrange the dataset with the things we need first.
#this will require some casting.
tailbeats <- trialdata %>%
  group_by(track) %>%
  filter(nodes %in% c("nose", "head", "spine1", "caudal"))
tailbeats <- tailbeats[-c(2)]

tailbeats <- unite(tailbeats, col='x-y', c('x', 'y'), sep='-') #group x and y together for manipulation
tailbeats <- spread(tailbeats, key=nodes, value="x-y") #spread out so that nodes are in their own column
tailbeats <- tailbeats %>% #put it in order
  arrange(track,frame)

#spread it back out so columns have their proper labels
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

################################################

#now for the actual calculus :(

#we'll need the vector between nose and head, and the vector between spine 1 and caudal:
tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(theta = acos( (((head.x - nose.x)*(caudal.x - spine1.x))+((head.y - nose.y)*(caudal.y - spine1.y)))/
                         (sqrt((head.x - nose.x)^2+(head.y - nose.y)^2)*sqrt((caudal.x - spine1.x)^2+(caudal.y - spine1.y)^2)) ) ) %>%
  mutate(thetadeg = rad2deg(theta))

#and now the plots:

#just for visualization's sake, I will filter this a little bit heh.

tailbeats <- tailbeats %>%
  group_by(track) %>%
  mutate(Behavior = na_if(Behavior, "S")) %>%
  mutate(Behavior = na_if(Behavior, "X")) %>%
  mutate(Behavior = na_if(Behavior, ""))

individ <- unique(tailbeats$track)
figuretitle <- paste('Individual', individ, sep="_")
thetalist <- split(tailbeats, tailbeats$track)

thetafacet <-  
  ggplot(tailbeats, aes(y=thetadeg, x=frame))+
  geom_point(mapping=aes(size=Behavior, color=Behavior, group=Behavior))+
  geom_path(linewidth=0.3)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle Relative to Head (Radians)")+
  scale_size_manual(values= c(2,2,2,2,0.5))+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)
thetafacet

ggsave(filename="thetafacet.png",
       plot=thetafacet, width = 13, 
       height = 7, units = "in",
       path=here("HotFigures"))

plottheta <- function(n) {
  ggplot(n, aes(y=thetadeg, x=frame))+
    geom_point(mapping=aes(color=Behavior, size=Behavior, group=Behavior))+
    geom_path(linewidth=0.3)+
    xlab("Time (Seconds)")+
    ylab("Tail Angle Relative to Head (Radians)")+
    scale_size_manual(values= c(5,5,5,5))+
    theme_classic()+
    theme(text = element_text(size = 18))  
}

thetaplotlist <- lapply(thetalist, plottheta)
thetaplotlist

names(thetaplotlist) <- figuretitle
lapply(names(thetaplotlist), 
       function(x) ggsave(filename=paste(x,"_tailangle_colorcoded",".png",sep=""),
                          plot=thetaplotlist[[x]], width = 10, 
                          height = 5, units = "in",
                          path=here("HotFigures")))

#of just the lateral displays:
lateraldisplays <- tailbeats[leadlag(tailbeats$Behavior == "L", 20, 20),]
lateraldisplay1 <- lateraldisplays %>%
  group_by(track) %>%
  slice(1:41)

lateraldisplayfacet <-  
  ggplot(lateraldisplay1, aes(y=thetadeg, x=frame))+
  geom_point(size=0.5)+
  geom_path(linewidth=0.3)+
  xlab("Time (Seconds)")+
  ylab("Tail Angle Relative to Head (Degrees)")+
  theme_classic()+
  theme(text = element_text(size = 12))+
  facet_wrap(~track, ncol=1)
lateraldisplayfacet

ggsave("lateraldisplayfacet.png", plot=lateraldisplayfacet, width=7, height=5,
       unit="in", path=here("HotFigures"))