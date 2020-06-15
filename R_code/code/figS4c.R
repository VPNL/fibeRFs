## Creates the basic plots for supplemental figure S4c 
# (Upper or lower field bias of ROIs)
# Note that .mat files are organized such that 1st col is ROI, 2nd col is subject
# All of EVC is included here so ROI indices are: 1 = V1, 2 = V2, 3 = V3, 
# 4 = IOG, 5 = pFus, 6 = mFus, 7 = pSTS, 8 = mSTS, 9 = CoS

rm(list=ls())
library(ggthemes)
library(R.matlab)
library(tidyverse)
library(plyr)

sem <- function(x) {sd(x, na.rm=TRUE) / sqrt(sum(!is.na((x))))}

#set results path
path <- "../results/study1/pRFs/"

files <- dir(paste0(path), 
             pattern = "*UpperVSLower_contraOnly_10_40.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  lh_lat_index = data$left.lat.index
  rh_lat_index = data$right.lat.index
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

##Create data frame
# Right
V1 <- rh_lat_index[1,]
IOG <- rh_lat_index[4,]
pFus <- rh_lat_index[5,]
mFus <- rh_lat_index[6,]
pSTS <- rh_lat_index[7,]
mSTS <- rh_lat_index[8,]
CoS <- rh_lat_index[9,]
hemi <- c("right","right","right","right","right",
          "right","right","right","right","right",
          "right","right","right","right","right",
          "right","right","right","right","right","right")
rh <- data.frame(subject,hemi,V1,IOG,pFus,mFus,pSTS,mSTS,CoS)

# Left
V1 <- lh_lat_index[1,]
IOG <- lh_lat_index[4,]
pFus <- lh_lat_index[5,]
mFus <- lh_lat_index[6,]
pSTS <- lh_lat_index[7,]
mSTS <- lh_lat_index[8,]
CoS <- lh_lat_index[9,]
hemi <- c("left","left","left","left","left",
          "left","left","left","left","left",
          "left","left","left","left","left",
          "left","left","left","left","left","left")
lh <- data.frame(subject,hemi,V1,IOG,pFus,mFus,pSTS,mSTS,CoS)

full <- full_join(rh, lh)

tidy <- gather(full, ROI, value, 3:9)
levels(tidy$ROI) <- c("V1","IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")
positions <- c("V1","IOG", "pFus", "mFus", "pSTS", "mSTS","CoS")

p <- ggplot(tidy, aes(x=ROI, y=value, col=ROI)) + 
  scale_x_discrete(limits = positions) +
  geom_bar(aes(fill=hemi),stat = "summary", fun.y = "mean", alpha = 0.9, position = "dodge",width=.75,show.legend=FALSE) + #avoid plotting outliers twice
  geom_point(aes(fill=hemi),position=position_dodge(width=.75), alpha = 0.75,size = 1.5,show.legend=FALSE) +
  theme_classic() +
  labs(x = "Region of interest", y = "Upper visual field - Lower visual field") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-1,1) +
  scale_colour_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1","#000000"), guide = 'none')  +
  scale_fill_manual(values=c("#C0C0C0","#696969")) 
p
ggsave("figs/S4c_ve10.pdf", p, width=7, height=5)

##Stats

# Some summary stats 
stable <- ddply(tidy, c("ROI", "hemi"), summarise,
                N    = length(value),
                mean = mean(value, na.rm=TRUE),
                se   = sem(value)
)
stable
## t-tests against 0
all_ps <- c(t.test(rh_lat_index[1,])$p.value,
            t.test(rh_lat_index[4,])$p.value,
            t.test(rh_lat_index[5,])$p.value,
            t.test(rh_lat_index[6,])$p.value,
            t.test(rh_lat_index[7,])$p.value,
            t.test(rh_lat_index[8,])$p.value,
            t.test(rh_lat_index[9,])$p.value,
            t.test(lh_lat_index[1,])$p.value,
            t.test(lh_lat_index[4,])$p.value,
            t.test(lh_lat_index[5,])$p.value,
            t.test(lh_lat_index[6,])$p.value,
            t.test(lh_lat_index[7,])$p.value,
            #t.test(lh_lat_index[8,])$p.value, #no mSTS
            t.test(lh_lat_index[9,])$p.value)
all_ps_corrected <- p.adjust(all_ps, method = "bonferroni")
all_ps_corrected
