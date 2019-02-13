rm(list=ls())
setwd("~/Desktop/MSResearch/BSAL_EDNA_X2")

library("readxl")
# Read in raw data from excel file 16
Plate_16_raw <- read_excel(path="Bsal_Plate_16_12152018_Resave.xls",
                           sheet=1,
                           skip=19)
# Select columns of interest and rows corresponding to Bsal
P16 <- subset(Plate_16_raw, select=c(1,4,5,6,7), subset=Target=="Bsal")
colnames(P16) <- c("Well", "Target", "Sample", "Starting_Quantity", "Cq")

# Remove wells corresponding to gBlock standards or extraction/PCR negatives
P16 <- P16[!grepl(pattern="gBlock", x=P16$Sample), ]
P16 <- P16[!grepl(pattern="NEG", x=P16$Sample, ignore.case=T), ]

# Import data from excel file 17
Plate_17_raw <- read_excel(path="Bsal_Plate_17_01092019_Resave.xls",
                           sheet=1,
                           skip=19)
P17 <- subset(Plate_17_raw, select=c(1,4,5,6,7), subset=Target=="Bsal")
colnames(P17) <- c("Well", "Target", "Sample", "Starting_Quantity", "Cq")

# In addition to standards/negatives, remove observations from X1C
P17 <- P17[!grepl(pattern="gBlock", x=P17$Sample), ]
P17 <- P17[!grepl(pattern="NEG", x=P17$Sample, ignore.case=T), ]
P17 <- P17[!grepl(pattern="X1C", x=P17$Sample, ignore.case=T), ]
P17$Sample <- gsub(pattern="X2.", replacement="", x=P17$Sample)

# Combine datasets
Bsal <- rbind(P16, P17)

# Aggregate by sample
Bsal.agg <- aggregate(x=Bsal$Starting_Quantity,
                      by=list(Sample=Bsal$Sample),
                      FUN=mean,
                      na.rm=T)
colnames(Bsal.agg)[2] <- "EstimatedCopies"

# Manually remove controls and samples that were re-processed
BsalPos <- Bsal.agg[-c(1:7, 10, 13, 15, 17, 19, 22, 24, 29, 33, 37, 42), ]

# Adjust estimated copies for dilutions
diluted.index <- grepl(pattern="D", x=BsalPos$Sample)
BsalPos$EstimatedCopiesAdj <-
  ifelse(grepl(pattern="D", x=BsalPos$Sample),
         BsalPos$EstimatedCopies*10,
         BsalPos$EstimatedCopies)
BsalPos$EstimatedCopies <- NULL

# Extract info from sample name and prune for filter data only
sample.expanded <- strsplit(BsalPos$Sample, ".", fixed=T)
BsalPos$FilterVolume <- sapply(sample.expanded,function(x) x[2])
BsalPos$ReplicateContainer <- sapply(sample.expanded,function(x) x[1])
BsalPos$ReplicateContainer <-
  gsub(x=BsalPos$ReplicateContainer, pattern="F", replacement="")

plot(x=BsalPos$FilterVolume,
     y=BsalPos$EstimatedCopiesAdj,
     #xlim=c(0, 1000),
     ylim=c(0, 4.3e6),
     xlab="Volume of Water Filtered",
     ylab="Estimated Copies")

plot(x=BsalPos$ReplicateContainer,
     y=BsalPos$EstimatedCopiesAdj,
     #xlim=c(0, 1000),
     ylim=c(0, 4.3e6),
     xlab="Replicate Container",
     ylab="Estimated Copies")

# Some alternate plots
library(tidyverse)
library(broom)

qplot(as.numeric(FilterVolume), EstimatedCopiesAdj, color=ReplicateContainer, data=BsalPos, geom="point") +
  geom_smooth(method = "lm", se=F) +
  facet_grid(ifelse(ReplicateContainer==1,"Y","N") ~ ., space="free_y", scales = "free_y")

qplot(as.numeric(FilterVolume), EstimatedCopiesAdj,
      color=ReplicateContainer,
      data=subset(BsalPos, ReplicateContainer != 1),
      geom=c("point", "line")) +
  geom_smooth(method = "lm", se=F)

# linear regression by container
lms <- BsalPos %>%
  group_by(ReplicateContainer) %>%
  do(mod = lm(EstimatedCopiesAdj ~ as.numeric(FilterVolume), data = .)) %>%
  tidy(mod, conf.int=TRUE)

ggplot(lms, aes(ReplicateContainer, y = estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() +
  facet_grid(term ~ ., scale="free_y")

