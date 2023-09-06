# Draft for sleap data pipeline

#actually, more low hanging fruit: chase stuff

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

#Make the nearest neighbor calculation more flexible!

#once again we want to focus on the head df. 
#let's assemble just the essential components of this df, as most of this we don't need

nndf <- head_df[-c(2,3,7:10)]
nndf$xy <- paste(nndf$x, nndf$y, sep='-')
nndf$track <- paste("track", nndf$track, sep='')
nndf <- nndf[-c(3,4)]