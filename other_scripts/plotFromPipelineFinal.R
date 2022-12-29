library(tidyverse)
library(pheatmap)
library(pvclust)

#setwd("/home/carmoma/projects/pollen/downloaded_bam/tmp/merge")

bad_samples <- c("FAR74611-1-NB02", "FAR74611-1-NB03", "FAR76967-1-NB01")
setwd("~/projects/pollen/results/mock1/11_countReadsPerReference/")
f <- list.files()
