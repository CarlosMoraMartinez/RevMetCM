library(tidyverse)
library(pheatmap)
library(pvclust)
#library(dendextend)
library(pvclust)
library(ggtext)
library(khroma)

mytheme <-  theme_bw()+
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5,
                                  colour = "black", face = "bold")) +
  theme(legend.title = element_text(face = "bold")) +
  theme(axis.text.y = element_text(size = 12, 
                                   colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.y = element_text(size = 12, 
                                    colour = "black", angle = 0, face = "bold")) +
  theme(strip.text.x = element_text(size = 12, 
                                    colour = "black", angle = 0, face = "bold")) +
  theme(axis.text.x = element_text( size = 14, 
                                   colour = "black", angle = 0, 
                                   face = "bold"))+
  theme(axis.title.x = element_text(vjust = 1, hjust = 0.5, 
                                    size = 16, colour = "black", 
                                    angle = 0, face = "bold")) +
  theme(axis.title.y= element_text(vjust = 1, hjust = 0.5, 
                                   size = 16, colour = "black", 
                                   angle = 90, face = "bold"))+
  theme(strip.background = element_blank())

read_all_bin_csvs <- function(recycle=TRUE){
  if(file.exists("ALL_READS_BINNED.csv") & recycle){
    cat("Reading file...")
    res <- read.table("ALL_READS_BINNED.csv", head=T, sep="\t", stringsAsFactors=F)
    return(res)
  }
  f <- list.files(pattern="_all.csv")
  res <- data.frame()
  for(ont in f){
    cat(ont, "\n")
    a <- read.table(ont, sep="\t", head=T, stringsAsFactors = F)
    a <- a %>% group_by(rname) %>% 
      top_n(n = 1, wt = coverage) %>% 
      ungroup
    b <- table(a$rname)
    b <- b[b > 1]
    a <- a[!a$rname %in% names(b), ]
    res <- rbind(res, a)
  }
  
  write.table(res, "ALL_READS_BINNED.csv", sep="\t", row.names = F, quote=F)
  return(res)
}

getSpeciesSummary<- function(res){
  res3 <- res %>% group_by(ONT) %>% 
    mutate(n_ont_including_unassigned = n()) %>% 
    filter(binned_species2 != "Unassigned") %>% 
    mutate(n_ont = n()) %>% 
    group_by(ONT, binned_species2) %>% 
    summarise(n_species = n(),
              n_ont = unique(n_ont),
              n_ont_including_unassigned = unique(n_ont_including_unassigned)) %>% 
    mutate(species = binned_species2,
           pct_species = n_species / n_ont) %>% 
    mutate(
      pct_species_trimmed = ifelse(pct_species < MIN_PERC_SP, 0, pct_species),
      pct_species_trimmed = pct_species_trimmed/sum(pct_species_trimmed),
      species_ref = species,
      species_num = as.numeric(gsub("AS0", "", species)),
      species = spnames$Species[match(species_num, spnames$EASI_ID)],
      mock_sample = mocksamplenum$Mock_commiunity[match(ONT, mocksamplenum$proj2)],
      mock_mix = mocksamplenum$mix[match(ONT, mocksamplenum$proj2)],
      dilution = mocksamplenum$rep[match(ONT, mocksamplenum$proj2)]
    )
  write.table(res3, "percents_norm.csv", sep="\t", row.names = F, quote = F)
  return(res3)
}

spread2matrixshape <- function(res3){
  res4 <- res3 %>% 
    select(mock_sample, species, pct_species_trimmed) %>% 
    spread(species, pct_species_trimmed)  %>% 
    data.frame %>% select(-ONT)
  res4$mock_sample <- as.character(res4$mock_sample)
  return(res4)
}

buildMatrix <- function(res4, mock_species2){
  mat <- res4 %>% select_if(is.numeric) %>% as.matrix
  mat[is.na(mat)] <- 0
  rownames(mat) <- res4$mock_sample
  colnames(mat) <- gsub("\\.", " ", colnames(mat)) %>% gsub(" $", ".", .)
  ord <- order(rownames(mat))
  mat <- mat[ord, ]
 return(mat)
}

makePvclust <- function(mat){
  zerocols <- apply(mat, MAR=2, any)
  mat2 <- mat[, zerocols]
  pvc <- pvclust(t(mat2), method.dist = "canberra", method.hclust = "ward.D2")
  pdf("pvclust.pdf", width = 6, height = 5)
  plot(pvc)
  dev.off()
  #tryCatch(dev.off(),cat(" device already closed. "))
}

makeAllHeatmaps<-function(mat){
  spsitalic <-sapply(
    colnames(mat),
    function(x) bquote(italic(.(x))))
  
  p <- pheatmap(t(mat), 
                #annotation_col = ann, 
                #annotation_row = ann_species,  
                filename = "heatmap_cluster.png", 
                labels_row = as.expression(spsitalic),
                angle_col=0,
                width = 9, height = 8)
  #dev.off()
  p <- pheatmap(t(mat), 
                #annotation_col = ann, 
                #annotation_row = ann_species,
                filename = "heatmap_bygroup.png", 
                labels_row = as.expression(spsitalic),
                cluster_cols = F,
                #gaps_col = seq(3, nrow(mat)-1, by=3),
                angle_col=0,
                width = 9, height = 8)
  
  mcc <- cor(mat)
  mcc[is.na(mcc)] <- 0
  spsitalic2 <-sapply(
    colnames(mcc),
    function(x) bquote(italic(.(x))))
  pheatmap(mcc,
           labels_row = as.expression(spsitalic2),
           filename = "heatmap_plant_correlation.png", 
           labels_col = as.expression(spsitalic2),
           width = 9, height = 8)
  #dev.off()
  
}

getBarPlot<- function(res3, mat){
  sps <- mat %>% apply(MAR=2, FUN=function(x)any(x==0)) %>% which %>% names
  res3$proportion <- res3$pct_species_trimmed
  res3$sample <- res3$mock_sample %>% as.factor
  g3 <- ggplot(res3, aes(y = proportion, x=sample, fill=species))+
    geom_bar(stat = "identity") +
    #scale_fill_tokyo(discrete=T)+
    mytheme
  ggsave("barplot_together.png", g3, width = 12, height = 7)
  
  ordersp <- apply(mat, MAR=2, mean) %>% order(decreasing = F)
  res3$species <- factor(res3$species, levels = colnames(mat)[ordersp])
  (g4 <- ggplot(res3, aes(y = proportion, x=species, fill=sample))+
    geom_bar(stat = "identity") +
    facet_grid(. ~sample)+
    scale_fill_muted()+
    mytheme +
    theme(axis.text.y = element_text(size = 12, 
                                      colour = "black", angle = 0, 
                                      face = "italic")) +
      theme(axis.text.x = element_text(size = 8, 
                                       colour = "black", angle = 0, 
                                       face = "bold"))+
      coord_flip() 
  )
  ggsave("barplot_bysample.png", g4, width = 14, height = 10)
  ggsave("barplot_bysample.pdf", g4, width = 12, height = 9)
}

calc_diversity<-function(res4){
  library(vegan)
  mat <- res4 %>% select(-c(mock_sample, condition)) %>% as.matrix()
  rownames(mat) <- res4$mock_sample
  res4$`Shannon diversity index` <- apply(mat, MAR=1, diversity, index="shannon")
  res4$sample <- res4$mock_sample %>% factor( levels=res4$sample[order(res4$sample)])
  #apply(mat, MAR=1, FUN=function(p)-sum(p[p!=0]*log(p[p!=0]))) # --> gives the same result
  gdiv <- ggplot(res4, aes(x=sample, y=`Shannon diversity index`,
                           fill=sample, col=sample, label=round(`Shannon diversity index`, 2)))+
    geom_segment( aes(x=sample, xend=sample, y=0, yend=`Shannon diversity index`)) +
    geom_point( size=4, alpha=0.7, shape=21, stroke=2) +
    ylim(0,4)+
    scale_color_muted()+
    scale_fill_muted()+
    geom_text(size=6, vjust= -1.5)+
    #geom_line(aes(col="black", fill="black"))+
    mytheme +
    theme(axis.title.x = element_text(vjust = 1, hjust = 0.5, 
                                      size = 12, colour = "black", 
                                      angle = 0, face = "bold")) +
    theme(axis.title.y= element_text(vjust = 1, hjust = 0.5, 
                                     size = 12, colour = "black", 
                                     angle = 90, face = "bold"))
  ggsave("shannon_index.pdf", gdiv, width = 9, height = 3)
  write.table(res4, "results_withShannonDiv.csv", sep="\t", quote=F, row.names = F)
}
getPCA<- function(mat){
  library(stats)
  library(factoextra)
  library(scatterplot3d)
  zerocols <- apply(mat, MAR=2, any)
  mat2 <- mat[, zerocols]
  pcacor <- prcomp(mat2, center=T, scale=T)
  
  grupo1 <-c( "Schedonorus pratensis", "Holcus lanatus", "Agrostis stolonifera", 
              "Phleum pratense", "Deschampsia cespitosa L.", 
              "Glyceria declinata", "Arrhenatherum elatius", 
              "Agrostis gigantea", "Agrostis capillaris", "Bromus hordeaceus", 
              "Agrostis canina", "Cynosurus cristatus L.")
  
  grupo2 <-c( "Poa pratensis", "Avellena flexuosa", "Glyceria maxima", "Molinia caerulea", "Donthania decumbuns", 
              "Ammophila arenaria L.", "Elymus repens")
  colorsvars <- ifelse(rownames(pcacor$rotation) %in% grupo1, "group 1", 
                   ifelse(rownames(pcacor$rotation) %in% grupo2, "group 2", "other"))
  pdf("PCA_scale.pdf")
  biplot(pcacor)
  
  fviz_pca_ind(pcacor, repel=T)
  fviz_eig(pcacor, choice="eigenvalue")
  fviz_eig(pcacor, choice="variance")
  fviz_pca_biplot(pcacor, axes = c(1,2), col.ind = "firebrick2",
                  repel = TRUE, col.var = "contrib") +
    scale_color_gradient(low="grey", high="dodgerblue2")  #gradient2 with low, mid and high
  fviz_pca_biplot(pcacor, axes = c(1,3), col.ind = "firebrick2",
                  repel = TRUE, col.var = "contrib") +
    scale_color_gradient(low="grey", high="dodgerblue2") 

  
  fviz_pca_biplot(pcacor, axes = c(1,2), col.ind = "black",
                  repel = TRUE, col.var = colorsvars)
    #scale_color_discrete(c("green","grey", "dodgerblue2"))  #gradient2 with low, mid and high
  fviz_pca_biplot(pcacor, axes = c(1,3), col.ind = "black",
                  repel = TRUE, col.var = colorsvars) 
    #scale_color_gradient(low="grey", high="dodgerblue2") 
  
   scatterplot3d(pcacor$x[,1], pcacor$x[,2], pcacor$x[,3],
                xlab="PC1", ylab="PC2", zlab="PC3",angle=40,pch=16)
   
   dev.off()
}
#setwd("/home/carmoma/projects/pollen/downloaded_bam/tmp/merge")

################################ READ DATA #############################

mock_species <- read_tsv("~/projects/pollen/METADATA/mock_composition.tsv")
splist <- mock_species$Spp
genus <- sapply(splist, FUN=function(x){strsplit(x, " ")[[1]][1]}) %>% unique

#Species with names in sequenced species
mock_species2 <- read_tsv("~/projects/pollen/METADATA/mock_composition2.tsv")
spnames <- read_tsv("~/projects/pollen/METADATA/species_name.tsv")
barcodes <- read_csv("~/projects/pollen/METADATA/barcodes_aerialsamples.tsv")

barcodes$SUMM <- gsub("_summary.csv", "",barcodes$SUMM) %>% gsub("_", "-", .)
barcodes$PROJ2 <- barcodes$PROJ
mocksamplenum <- barcodes #read_tsv("~/projects/pollen/METADATA/easi_samplenum.tsv")
#mocksamplenum$rep <- gsub("[0-9]+", "", mocksamplenum$Mock_commiunity)
#mocksamplenum$mix <- paste("mix_", gsub("[A-Z]", "", mocksamplenum$Mock_commiunity), sep="")
mocksamplenum$Mock_commiunity <-mocksamplenum$PROJ2 # barcodes$PROJ2[match(mocksamplenum$EASI_ID, barcodes$ASSAMPLE)]
mocksamplenum$proj2 <-barcodes$SUMM 

setwd("/home/carmoma//projects/pollen/results/")
# Make sure to perform sed -i "s/#//"  before reading tables, there is a # name in header
dirlist <- c( 
  "/home/carmoma/projects/pollen/results/problem_samples_1/binresults01",
  "/home/carmoma/projects/pollen/results/problem_samples_1/binresults05",
  "/home/carmoma/projects/pollen/results/problem_samples_1/binresults10",
  "/home/carmoma/projects/pollen/results/problem_samples_1/binresults15")

PERC_LIM <- 0.01
MIN_PERC_SP <- 0.01

#REMOVE READS WITH LESS THAN 15% OF THEIR LENGTH COVERED

all_res4 <- data.frame()
for(resultsdir in dirlist){
  cat(resultsdir)
  setwd(resultsdir)
  res <- read_all_bin_csvs(recycle=TRUE)
  res$binned_species2 <- mapply(res$ILLUMINA, res$coverage, FUN=function(sp, pc) ifelse(pc < PERC_LIM, "Unassigned", sp))

  #CALCULATE PERCENTAGE OF SPECIES
  res3 <- getSpeciesSummary(res)
  res4 <- spread2matrixshape(res3)
  mattemp <- buildMatrix(res4, mock_species2)
  res4$condition <- resultsdir 
  write.table(res4, "results_matrix.csv", sep="\t", quote=F, row.names = F)
  all_res4 <- rbind(all_res4, res4)
  
  res2write <- as.data.frame(round(mattemp*100, 3))
  write.table(res2write, "results_matrix_100.csv", sep="\t", row.names = T, quote=F)

  #makePvclust(mat)
  getPCA(mat)
  makeAllHeatmaps(mattemp)
  #getDotPlot(res4, mattemp)
  getBarPlot(res3, mattemp)

}


setwd("/home/carmoma//projects/pollen/results/")
#all_prev <- read.table("230101_allresults3.csv", sep="\t", stringsAsFactors=F, head=T)
#all_res4 <- rbind(all_prev, all_res4)
write.table(all_res4, "230107_allresults_samples.csv", sep="\t", quote=F, row.names = F)

#auxsum <- all_res4 %>% 
#  group_by(condition) %>% 
#  select(correlation, 
#        correlation_present, 
#        mean_diff, mean_diff_positive,
 #       mean_diff_negative, 
#        max_diff_positive, 
#        max_diff_negative, 
#        TP, TN, FP, FN, 
#        specificity, sensitivity, PPV, NPV
#    ) %>% 
#  summarise_all(mean)
#3write.table(auxsum, "230107_allsummary2samples.csv", sep="\t", quote=F, row.names = F)
