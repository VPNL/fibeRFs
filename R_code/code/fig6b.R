## Creates figure 6B (per ROI fiber comparison)

rm(list=ls())
library(ggthemes)
library(R.matlab)
library(tidyverse)
library(plyr)

sem <- function(x) {sd(x, na.rm=TRUE) / sqrt(sum(!is.na((x))))}

#set results path
path <- "../results/study1/fibers/"

files <- dir(paste0(path), 
             pattern = "*_fibers_overall_per_ROI.mat")

for (f in files) {
  mf <- paste0(path,f)
  
  data <- readMat(mf)
  
  if (grepl("lh",mf)) {
    lh_fibers = data$overall
  } else if (grepl("rh",mf)) {
    rh_fibers = data$overall
  }
}

subject <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)

ROI <- rep(c("IOG"), times = 16)
right <- rh_fibers[,1]
left <- lh_fibers[,1]
IOG <- data.frame(subject,ROI,left,right)
IOG_gather <- gather(IOG,"hemi","proportion", 3:4)
ROI <- rep(c("pFus"), times = 16)
right <- rh_fibers[,2]
left <- lh_fibers[,2]
pFus <- data.frame(subject,ROI,left,right)
pFus_gather <- gather(pFus,"hemi","proportion", 3:4)
ROI <- rep(c("mFus"), times = 16)
right <- rh_fibers[,3]
left <- lh_fibers[,3]
mFus <- data.frame(subject,ROI,left,right)
mFus_gather <- gather(mFus,"hemi","proportion", 3:4)
ROI <- rep(c("pSTS"), times = 16)
right <- rh_fibers[,4]
left <- lh_fibers[,4]
pSTS <- data.frame(subject,ROI,left,right)
pSTS_gather <- gather(pSTS,"hemi","proportion", 3:4)
ROI <- rep(c("mSTS"), times = 16)
right <- rh_fibers[,5]
left <- lh_fibers[,5]
mSTS <- data.frame(subject,ROI,left,right)
mSTS_gather <- gather(mSTS,"hemi","proportion", 3:4)
ROI <- rep(c("CoS"), times = 16)
right <- rh_fibers[,6]
left <- lh_fibers[,6]
CoS <- data.frame(subject,ROI,left,right)
CoS_gather <- gather(CoS,"hemi","proportion", 3:4)


test <- bind_rows(IOG_gather, pFus_gather, mFus_gather, pSTS_gather, mSTS_gather, CoS_gather)
levels(test$ROI) <- c("IOG", "pFus", "mFus", "pSTS", "mSTS", "CoS")
positions <- c("IOG", "pFus", "mFus", "pSTS", "mSTS", "CoS")

## Key plotting
b <- ggplot(test, aes(x=ROI, y=proportion, col=ROI)) + 
  scale_x_discrete(limits = positions) +
  scale_y_continuous(limits = c(0,0.9), breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) +
  geom_bar(aes(fill=hemi),stat = "summary", fun.y = "mean", alpha = 0.9, position = "dodge",width=.75) + #avoid plotting outliers twice
  geom_point(aes(fill=hemi),position=position_dodge(width=.75), alpha = 0.75,size = 1.5,show.legend=FALSE) +
  theme_classic() +
  labs(title = "Proportion of all fROI fibers connecting to EVC",
       y = "Proportion", x = "ROI") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values=c("#009900","#d73027", "#fdae61",  "#4575b4", "#f46d43",  "#74add1"), guide = 'none')  +
  scale_fill_manual(values=c("#C0C0C0","#696969"), guide = 'none') 
b
ggsave("figs/6b.pdf", b, width=5, height=5)

# Some summary stats 
cdata <- ddply(test, c("ROI", "hemi"), summarise,
               N    = length(proportion),
               mean = mean(proportion, na.rm=TRUE),
               se   = sem(proportion)
)
cdata

