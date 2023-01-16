library(tidyverse)
library(pheatmap)
library(pvclust)
#library(dendextend)
library(pvclust)
library(ggtext)
library(khroma)

#THIS IS THE MAIN SCRIPT FOR THE ANALYSIS OF REVMET RESULTS
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
  theme(axis.text.x = element_text(size = 11, 
                                   colour = "black", angle = 0, 
                                   face = "bold"))+
  theme(axis.title.x = element_text(vjust = 1, hjust = 0.5, 
                                    size = 16, colour = "black", 
                                    angle = 0, face = "bold")) +
  theme(axis.title.y= element_text(vjust = 1, hjust = 0.5, 
                                   size = 16, colour = "black", 
                                   angle = 90, face = "bold"))

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
    select(mock_sample, mock_mix, dilution, species, pct_species_trimmed) %>% 
    spread(species, pct_species_trimmed)  %>% 
    data.frame %>% select(-ONT)
  return(res4)
}

matMockData <- function(res4, mat, mock_species2){
  matextra <- matrix(numeric(length(unique(res4$mock_mix))*ncol(mat)), nrow=ncol(mat))
  rownames(matextra) <- colnames(mat)
  for (n in mock_species2$Spp){
    matextra[n, ] <- unlist(mock_species2[mock_species2$Spp == n, names(mock_species2) %in% res4$mock_mix])
  }
  colnames(matextra) <- as.factor(as.character(1:10))
  matextra <- as.data.frame(matextra)
  return(matextra)
}

#Esta tiene la matriz de todas las plantas secuenciadas, pero Agrostis agrupadas en una sola. 
matMockData_1 <- function(res4, mat, mock_species){
  spnames <- c("Agrostis", colnames(mat)[!grepl("Agrostis", colnames(mat))])
  present_species <- mock_species$Spp[mock_species$Spp %in% spnames]
  matextra_1 <- matrix(numeric(length(unique(res4$mock_mix))*length(spnames)), nrow=length(spnames))
  rownames(matextra_1) <- spnames
  for (n in present_species){
    matextra_1[n, ] <- unlist(mock_species[mock_species$Spp == n, names(mock_species) %in% res4$mock_mix])
  }
  colnames(matextra_1) <- as.factor(as.character(1:10))
  matextra_1 <- as.data.frame(matextra_1)
  return(matextra_1)
}

getSampleAnnotation <- function(res4, mat){
  ann <- data.frame(mix = as.factor(gsub("mix_", "", res4$mock_mix)), 
                    dilution = res4$dilution)
  names(ann)[1] <- "mock group"
  rownames(ann) <- rownames(mat)  
  return(ann)
}

getMergedMat2 <- function(mat, matextra, ann){
  mat2 <- t(mat) %>% as.data.frame
  colnames(mat2) <- paste("sample ", colnames(mat2), sep="")
  matextra2 <- matextra
  colnames(matextra2) <- paste("mix ", colnames(matextra), sep = "")

  mergedmat <- data.frame(species = rownames(mat2))
  for(g in colnames(matextra2)){
    mergedmat[, g] <- matextra2[, g]
    samples <- which(as.character(ann$`mock group`) == gsub("mix ", "", g))
    mergedmat <- cbind(mergedmat, mat2[, samples])
  }
  mergedmat2 <- as.matrix(mergedmat[, 2:ncol(mergedmat)])
  return(mergedmat2)
}

getGenusSpecies<- function(mergedmat2){
  genus_species <- sapply(rownames(mergedmat2), FUN=function(x){strsplit(x, " ")[[1]][1]})
  ann_species <- data.frame(present = ifelse(rownames(mergedmat2) %in% splist, "species", 
                                           ifelse(genus_species %in% genus, "genus", "Other (false positive)")
  ))
  rownames(ann_species) <- rownames(mergedmat2)
  return(ann_species)
}

buildMatrix <- function(res4, mock_species2){
  mat <- res4 %>% select_if(is.numeric) %>% as.matrix
  mat[is.na(mat)] <- 0
  rownames(mat) <- res4$mock_sample
  
  colnames(mat) <- gsub("\\.", " ", colnames(mat)) %>% gsub(" $", ".", .)
  ann <- getSampleAnnotation(res4, mat)
  
  ord <- order(rownames(mat))
  ord <- ord[c(4:length(ord), 1:3)]
  mat <- mat[ord, ]
  ann <- ann[ord, ]
  matextra <- matMockData(res4, mat, mock_species2)
  matextra_1 <- matMockData_1(res4, mat, mock_species)
  mergedmat2 <- getMergedMat2(mat, matextra, ann)
  ann_species <- getGenusSpecies(mergedmat2)
  
  ann2 <- data.frame(dilution = rep(c("theoretical %", "dil A", "dil B", "dil C"), 
                                    dim(mergedmat2)[2]/4))
  rownames(ann2) <- colnames(mergedmat2)
  matcut <- mat[ , colnames(mat) %in% mock_species2$Spp]
  ann_species_cut <- ann_species[colnames(matcut), ]
  return(list("mat" = mat, "ann"=ann, 
              "matextra"=matextra, 
              "matextra_1"=matextra_1, 
              "mergedmat2"=mergedmat2, 
              "ann_species"=ann_species, 
              "ann2"=ann2,
              "matcut"=matcut,
              "ann_species_cut"=ann_species_cut)
  )
}

makePvclust <- function(mat){
  pvc <- pvclust(t(mat), method.dist = "canberra")
  pdf("pvclust_canberra.pdf")
  plot(pvc)
  dev.off()
  #tryCatch(dev.off(),cat(" device already closed. "))
}

makeAllHeatmaps<-function(mattemp){
  mat <- mattemp[["mat"]]
  ann <- mattemp[["ann"]]
  matextra <- mattemp[["matextra"]]
  mergedmat2 <- mattemp[["mergedmat2"]]
  ann_species <- mattemp[["ann_species"]]
  ann2 <- mattemp[["ann2"]]
  matcut <- mattemp[["matcut"]]
  ann_species_cut <- mattemp[["ann_species_cut"]]
  
  spsitalic <-sapply(
    colnames(mat),
    function(x) bquote(italic(.(x))))
  
  p <- pheatmap(t(mat), 
                annotation_col = ann, 
                annotation_row = ann_species,  
                filename = "heatmap_cluster.png", 
                labels_row = as.expression(spsitalic),
                angle_col=315,
                width = 12, height = 8)
  #dev.off()
  p <- pheatmap(t(mat), 
                annotation_col = ann, 
                annotation_row = ann_species,
                filename = "heatmap_bygroup.png", 
                labels_row = as.expression(spsitalic),
                cluster_cols = F,
                gaps_col = seq(3, nrow(mat)-1, by=3),
                angle_col=315,
                width = 12, height = 8)
  #dev.off()
  #Only species of interest
  
  p <- pheatmap(t(matcut), 
                annotation_col = ann, 
                #annotation_row = ann_species_cut,
                filename = "heatmap_cluster_selectedspecies.png", 
                labels_row = as.expression(spsitalic),
                angle_col=315,
                width = 12, height = 5)
  #dev.off()
  p <- pheatmap(t(matcut), 
                annotation_col = ann, 
                #annotation_row = ann_species_cut,
                filename = "heatmap_bygroup_selectedspecies.png", 
                labels_row = as.expression(spsitalic),
                cluster_cols = F,
                gaps_col = seq(3, nrow(mat)-1, by=3),
                angle_col=315,
                width = 12, height = 5)
  #dev.off()
  
  p <- pheatmap(mergedmat2, 
                annotation_col = ann2, 
                annotation_row = ann_species,
                filename = "heatmap_with_originalmix.png", 
                labels_row = as.expression(spsitalic),
                cluster_cols = F,
                gaps_col = seq(4, nrow(mergedmat2)-1, by=4),
                angle_col=315,
                width = 14, height = 8)
  #dev.off()
}

getCorrelations<-function(res4, mattemp, mock_species, mock_species2){
  mat <- mattemp[["mat"]]
  matextra3 <- mattemp[["matextra_1"]]
  
  matextra3[grepl("Poa", rownames(matextra3)),] <- 0 #These species are not actually there
  
  #Average agrostis
  Agrostis <- mat[, grepl("Agrostis", colnames(mat))] %>% apply(MAR=1, sum)
  mat <- cbind(mat[, ! grepl("Agrostis", colnames(mat))], Agrostis)
  
  #agrostis_theo <- matextra3["Agrostis canina", ]
  #rownames(agrostis_theo) <- "Agrostis"
  #matextra3 <- rbind(matextra3[! grepl("Agrostis", rownames(matextra3)), ], agrostis_theo)
  
  colnames(matextra3) <- paste("mix_", colnames(matextra3), sep="")
  mat <- mat[,rownames(matextra3)]
  #colnames(mat) == rownames(matextra3)
  present_species <- mock_species$Spp[mock_species$Spp %in% rownames(matextra3)]
  species_absent <- rownames(matextra3)[! rownames(matextra3) %in% present_species]
  
  for(g in colnames(matextra3)){
    ind = which(res4$mock_mix == g)
    for(i in ind){
      res4$correlation[i] <- cor(matextra3[, g], mat[res4$mock_sample[i], ])
      res4$correlation_present[i] <- cor(matextra3[present_species, g], mat[res4$mock_sample[i], present_species])
      res4$total_negative[i] <- sum(mat[res4$mock_sample[i], species_absent])      
      res4$mean_diff[i] <-  mean(abs(matextra3[, g] - mat[res4$mock_sample[i], ]))
      res4$mean_diff_positive[i] <- mean(abs(matextra3[present_species, g] - 
                                                  mat[res4$mock_sample[i], present_species]))
      res4$max_diff_positive[i] <- max(abs(matextra3[present_species, g] - 
                                               mat[res4$mock_sample[i], present_species]))

      res4$mean_diff_negative[i] <- mean(abs(matextra3[species_absent, g] - 
                                          mat[res4$mock_sample[i], species_absent]))
      res4$max_diff_negative[i] <- max(abs(matextra3[species_absent, g] - 
                                            mat[res4$mock_sample[i], species_absent]))

      
      all_plants <- 1:nrow(mat)
      plants_found <- which(mat[res4$mock_sample[i], ] > 0)
      plants_in <- which(matextra3[, g] > 0)
      tp = length(which(plants_found %in% plants_in))
      tn = length(which(! all_plants %in% c(plants_found, plants_in)))
      fp = length(which(! plants_found %in% plants_in))
      fn = length(which( all_plants %in% plants_in & ! all_plants %in% plants_found))
      
      res4$TP[i] <- tp
      res4$TN[i] <- tn
      res4$FP[i] <- fp
      res4$FN[i] <- fn
      res4$specificity[i] <-  100 * tn/(tn + fp)
      res4$sensitivity[i] <- 100 * tp/(tp+fn)
      res4$PPV[i] <-100 * tp /(tp+fp)
      res4$NPV[i] <- 100 * tn/(tn+fn)
    }
  }
  return(res4)
}

getDotPlot<- function(res4, mattemp){
  matextra <- mattemp[["matextra"]]
  matextra[grepl("Poa", rownames(matextra)),] <- 0
  rn <- rownames(matextra)

  colnames(matextra) <- paste("mix_", colnames(matextra), sep="")
  
  names(res4) <- gsub("\\.", " ", names(res4)) %>% gsub(" $", ".", .)
  resvert <- res4 %>% select(mock_sample, mock_mix, dilution, rownames(matextra)) %>% 
    gather(species, proportion_found, rownames(matextra))
  resvert$proportion_expected <- mapply(resvert$species, resvert$mock_mix, 
                                        FUN=function(sp, mix)matextra[sp, mix])
  resvert$proportion_found[is.na(resvert$proportion_found)] <- 0
  resvert$mock_mix <- factor(gsub("_", " ", as.character(resvert$mock_mix)))
  write.table(resvert, "results_long_withexpected.csv", sep="\t", row.names = T, quote=F)
  
  lev <- sort(unique(resvert$mock_sample))
  lev <- lev[c(4:length(lev), 1:3)]
  resvert$species_mod <- ifelse(resvert$proportion_expected > 0, 
                                resvert$species, "Other (false positive)")
  resvert$mock_sample <- factor(resvert$mock_sample, 
                                   levels = lev)
  
  lev <- sort(unique(resvert$mock_mix))
  lev <- lev[c(1, 3:length(lev), 2)]
  resvert$mock_mix <- factor(resvert$mock_mix, lev)
  
  #splabs <- paste("*",unique(resvert$species_mod), "*", sep="")
 # splabs <- ifelse(splabs=="*absent*", "Other (false positive)", splabs)
  agr <- resvert %>% filter(grepl("Agrostis", species)) %>% 
    group_by(mock_sample, mock_mix, dilution) %>% 
    summarise(species = "Agrostis total", species_mod = "Agrostis total",
              proportion_found = sum(proportion_found),
              proportion_expected = unique(proportion_expected))
  agr <- agr[, names(resvert)]
  
  g1 <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                      col=species_mod, fill=species_mod))+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_colour_muted() +
    mytheme
  ggsave("corrplot1_sepNone_colSpecies.png", g1, width = 10, height = 7)
  
  g1b <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                            col=mock_mix, fill=mock_mix, shape=dilution))+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    mytheme
  ggsave("corrplot1_sepNone_colSample.png", g1b, width = 10, height = 7)
  
  g2 <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                            col=species_mod, fill=species_mod))+
    facet_wrap(.~ mock_sample)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_colour_muted() +
    #scale_fill_discrete(
    #  "Species",
    #  breaks = c(0:(length(splabs)-1)),
    #  labels = splabs
    #) +
    mytheme +
    theme(axis.text.x = element_text(size = 9, 
                                     colour = "black", angle = 0, 
                                     face = "bold"))+
    theme(axis.text.y = element_text( size = 9, 
                                     colour = "black", angle = 0, 
                                     face = "bold"))
    #theme(legend.text = element_markdown()) 
  ggsave("corrplot2_sepSample_colSpecies.png", g2, width = 12, height = 7)
  
  g2b <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                             col=mock_mix, fill=mock_mix, shape=dilution))+
    facet_wrap(.~ species_mod)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_fill_tofino(discrete=T)+
    mytheme +
    theme(strip.text.y = element_text(size = 12, 
                                      colour = "black", angle = 0, face = "italic")) 
  ggsave("corrplot2b_sepSpecies_colSample.png", g2b, width = 10, height = 7)

  resvertAgr <- rbind(as.data.frame(resvert),as.data.frame(agr)) %>% filter(species_mod != "Other (false positive)")
  g2b2 <- ggplot(resvertAgr, aes(x=proportion_expected, y=proportion_found, 
                             col=mock_mix, fill=mock_mix, shape=dilution))+
    facet_wrap(.~ species_mod)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_fill_tofino(discrete=T)+
    mytheme +
    theme(strip.text.y = element_text(size = 12, 
                                      colour = "black", angle = 0, face = "italic")) 
  ggsave("corrplot2b_sepSpecies_colSample_AgrostisSum.png", g2b2, width = 10, height = 7)
  
    
  g2c <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                             col=mock_mix, fill=mock_mix, shape=dilution))+
    facet_wrap(.~ species)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    mytheme
  ggsave("corrplot2c_sepSpeciesAll_colSample.png", g2c, width = 16, height = 16)
  
  g2d <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                            col=species_mod, fill=species_mod, grou, shape=dilution))+
    facet_wrap(.~ mock_mix, nrow=2)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_colour_muted() +
    mytheme +
    theme(axis.text.x = element_text(size = 10, 
                                     colour = "black", angle = 0, 
                                     face = "bold"))+
    theme(axis.text.y = element_text( size = 10, 
                                     colour = "black", angle = 0, 
                                     face = "bold"))
  ggsave("corrplot2d_sepMix_colSpecies.png", g2d, width = 12, height = 6)
  
  g2e <- ggplot(resvert, aes(x=proportion_expected, y=proportion_found, 
                             col=species_mod, fill=species_mod, grou))+
    facet_wrap(.~ dilution)+
    geom_point() +
    geom_abline(slope = 1)+
    xlab("Expected proportion") +
    ylab("Obtained proportion") +
    xlim(0, max(resvert$proportion_expected+0.05)) +
    ylim(0, max(resvert$proportion_found+0.05)) +
    scale_colour_muted() +
    mytheme+
    theme(axis.text.x = element_text(size = 11, 
                                     colour = "black", angle = 0, 
                                     face = "bold"))+
    theme(axis.text.y = element_text( size = 11, 
                                      colour = "black", angle = 0, 
                                      face = "bold"))
  ggsave("corrplot2e_sepDil_colSpecies.png", g2e, width = 10, height = 4)
  
  #Barplot combining Agrostis
  resvert2 <- resvert %>% 
    gather(type, proportion, proportion_expected, proportion_found) %>% 
    distinct() %>%
    mutate(mock_sample2 = ifelse(type == "proportion_expected", as.character(mock_mix), as.character(mock_sample)),
           dilution = ifelse(type == "proportion_expected", "Exp", 
                             as.character(dilution)),
           )
  #Yes, the operation is repeated
  rag <- resvert2[grepl("Agrostis", resvert2$species_mod), ]
  rag$species <- "Agrostis"
  rag$species_mod <- "Agrostis"
  
  rag_exp <- rag %>% filter(type=="proportion_expected") %>% distinct()
  rag_obs<- rag %>% filter(type!="proportion_expected") %>% 
    group_by(mock_sample, mock_mix, dilution, species, species_mod, type, mock_sample2) %>% 
    summarise(proportion=sum(proportion))
  rag_obs <- rag_obs[, names(rag)]
  rag <- rbind(as.data.frame(rag_obs), as.data.frame(rag_exp))
  resvert3 <- rbind(resvert2[!grepl("Agrostis", resvert2$species_mod), ], rag)
  #Re-normalize expected pcts
  rev_allexp <- resvert3[resvert3$type == "proportion_expected", ] %>% 
    group_by(mock_mix) %>% 
    mutate(proportion = proportion/(sum(proportion))) %>% 
    as.data.frame()
  if(any(is.na(rev_allexp$proportion))){
    rev_allexp$proportion[is.na(rev_allexp$proportion)] <- 0
  }
  resvert4 <- rbind(resvert3[resvert3$type != "proportion_expected", ] , rev_allexp)
  
  g3 <- ggplot(resvert4, aes(y = proportion, x=dilution, fill=species_mod))+
    facet_wrap(. ~ mock_mix, scales = "free_x")+
    geom_bar(stat = "identity") +
    scale_fill_bright()+
    mytheme
  ggsave("barplot.png", g3, width = 10, height = 7)
  
  resvert_f <- resvert4 %>% filter(species_mod == "Other (false positive)" & dilution != "Exp")
  mins <- resvert_f %>% group_by(species) %>% summarise(maxp = max(proportion)) %>% 
    filter(maxp > 0)
  resvert_f <- resvert_f %>% filter(species %in% mins$species)
  
  g4 <- ggplot(resvert_f, aes(y = proportion, x=dilution, fill=species))+
    facet_wrap(. ~ mock_mix, scales = "free_x")+
    geom_bar(stat = "identity") +
    scale_fill_tofino(discrete=T)+
    mytheme
  ggsave("barplot_falsePositives.png", g4, width = 10, height = 7)
  
  #Finally, dot plots:
  present <- resvert$species[resvert$proportion_expected > 0] %>% unique
  genus_present <-  sapply(present, FUN=function(x)strsplit(x, " ")[[1]][1]) %>% unique
  resvert$`genus` <- sapply(resvert$species, FUN=function(x)strsplit(x, " ")[[1]][1])
  resvert$present <- ifelse(resvert$species %in% present & resvert$genus != "Agrostis",
                            "present", 
                            ifelse(resvert$genus %in% genus_present, "genus present", "absent"))
  resvert$`plant present` <- ifelse(resvert$proportion_expected==0 & resvert$proportion_found > 0,
                                   "false positive", "true positive")
  resvert$orderval <- ifelse(resvert$present=="absent",
                             -resvert$proportion_found, resvert$proportion_found)
  resord <- resvert %>% group_by(species) %>% summarise(m=mean(orderval))
  resord <- resord[order(resord$m), ]
  resvert$species <- factor(resvert$species, levels=resord$species)
  
  g5 <- ggplot(resvert, aes(y = species,x = mock_sample)) +     
    facet_grid(.~ mock_mix, scales="free_x")+
    geom_tile(fill="white") + 
    geom_point(aes(colour = `plant present`, 
                   size =proportion_found))  +  
    scale_size(range = c(0, 8))+             ## to tune the size of circles
    scale_colour_manual(values = c("firebrick3", "dodgerblue2")) +
    mytheme +
    theme(axis.text.y = element_text( face = "italic")) +
    theme(
          strip.background = element_blank())
  ggsave("g5_dotplot1_byGroup.pdf", g5, width = 12, height = 8)
  
  resvert_theo <- resvert %>% group_by(mock_sample, mock_mix, dilution, 
                                       species, species_mod, genus, present, 
                                        `plant present`) %>% 
    summarise(proportion_found = unique(proportion_expected), 
              proportion_expected = unique(proportion_expected),
              orderval = mean(orderval)) %>% 
              mutate(`plant present` = "true positive") %>% 
    mutate(mock_sample = mock_mix, type = "theoretical",
           dilution="Exp.")
  #resvert$`plant present`[resvert$proportion_expected == 0] <- "true positive"
  resvert$type <- "retrieved"
  resvert5 <- rbind(as.data.frame(resvert), as.data.frame(resvert_theo))
  
  lev <- unique(resvert5$mock_sample)
  lev <- lev[c(31:40, 1:30)]
  resvert5$mock_sample <- factor(resvert5$mock_sample, levels = lev)
  resvert5$dilution <- factor(resvert5$dilution, levels = c("Exp.", "A", "B", "C"))
  
  g6 <- ggplot(resvert5, aes(y = species,x = dilution, fill=type)) +     
    facet_grid(.~ mock_mix, scales="free_x")+
    geom_tile() + 
    scale_fill_manual(values = c("white", "grey83")) +
    geom_point(aes(colour = `plant present`, 
                   size =proportion_found))  +  
    scale_size(range = c(0, 8))+             ## to tune the size of circles
    scale_colour_manual(values = c("firebrick3", "dodgerblue2")) +
    mytheme +
    theme(axis.text.y = element_text( face = "italic")) +
    theme(strip.background = element_blank())
  ggsave("g6_dotplot1_byGroupWithTheo.pdf", g6, width = 13, height = 9)
  
  datamod <- resvert[as.character(resvert$species) %in% present,c("mock_sample", "mock_mix", 
                                                                  "dilution", "species", 
                                                                  "proportion_found", "proportion_expected")]
  agr <- datamod[grepl("Agrostis", datamod$species),]
  agr <- agr %>% 
    mutate(species = "Agrostis") %>% 
    group_by(mock_sample, mock_mix, dilution, species) %>% 
    summarise(proportion_found = sum(proportion_found),
           proportion_expected = unique(proportion_expected))
  noagr <-  datamod[!grepl("Agrostis", datamod$species),]
  datamod2 <- rbind(as.data.frame(agr), as.data.frame(noagr))
  datamod2$diff <- datamod2$proportion_expected - datamod2$proportion_found
  mod <- lm(diff ~  dilution + species,  data=datamod2)
  
  mod1 <- lm(proportion_found ~ proportion_expected + dilution * species,  data=datamod2) #R^2 0.71 
  mod2 <- lm(proportion_found ~ proportion_expected + dilution + species,  data=datamod2) #R^2 0.71 
  mod3 <- lm(proportion_found ~ proportion_expected  + species,  data=datamod2)#R^2 0.70 Adjusted
  mod4 <- lm(proportion_found ~ proportion_expected ,  data=datamod2) #R^2 0.6
  mod5 <- lm(proportion_found ~ dilution + proportion_expected ,  data=datamod2) #R^2 0.59 Adjusted
  mod6 <- lm(proportion_found ~ species ,  data=datamod2) #R^2 0.14 adjusted
 
  mod7 <- lm(proportion_found ~  dilution + species,  data=datamod2) #R^2 0.13
  mod8 <- lm(proportion_found ~  dilution,  data=datamod2) #R^2 0

  a0<-anova(mod1, mod2)# no sig -> factor dilucion no importante  
  a1<-anova(mod2, mod3)# no sig -> factor dilucion no importante
  a2<-anova(mod3, mod4) # sig -> factor espcie importante
  a3<-anova(mod3, mod6) #sig -> factor esperado importante
  
  datamod2$predicted <- predict(mod1)
  datamod2$diff <- abs(datamod2$proportion_found - datamod2$predicted)
  
  datamod_summary <- datamod2 %>% group_by(mock_sample, mock_mix, dilution) %>% 
    summarise(mean_difference = mean(diff), 
              correlation = cor(proportion_found, predicted))
  write.table(datamod_summary, file="DifferenceAndCorrelationToPrediction.csv", sep="\t", row.names = F, quote=F)
  models <- list(mod1, mod2, mod3, mod4, mod5, mod6, a1, a2, a3, datamod2, datamod_summary)
  mean(datamod_summary$mean_difference) #0.05057
  mean(datamod_summary$correlation) #0.8494

  save(file="RegressionModels.RData", models)
  #mod1: interaction significant for all
}
#setwd("/home/carmoma/projects/pollen/downloaded_bam/tmp/merge")

################################ READ DATA #############################

mock_species <- read_tsv("~/projects/pollen/METADATA/mock_composition.tsv")
splist <- mock_species$Spp
genus <- sapply(splist, FUN=function(x){strsplit(x, " ")[[1]][1]}) %>% unique

#Species with names in sequenced species
mock_species2 <- read_tsv("~/projects/pollen/METADATA/mock_composition2.tsv")
spnames <- read_tsv("~/projects/pollen/METADATA/species_name.tsv")
barcodes <- read_csv("~/projects/pollen/METADATA/barcodes.tsv")
mocksamplenum <- read_tsv("~/projects/pollen/METADATA/easi_samplenum.tsv")
barcodes$PROJ <- gsub("_summary.csv", "",barcodes$PROJ) %>% gsub("_", "-", .)
barcodes$PROJ2 <- paste(barcodes$PROJ, barcodes$YY, sep="-")
mocksamplenum$rep <- gsub("[0-9]+", "", mocksamplenum$Mock_commiunity)
mocksamplenum$mix <- paste("mix_", gsub("[A-Z]", "", mocksamplenum$Mock_commiunity), sep="")
mocksamplenum$proj2 <- barcodes$PROJ2[match(mocksamplenum$EASI_ID, barcodes$ASSAMPLE)]

setwd("/home/carmoma//projects/pollen/results/")
# Make sure to perform sed -i "s/#//"  before reading tables, there is a # name in header
dirlist <- c( 
  "/home/carmoma/projects/pollen/results/nuclear1/binresults01",
  "/home/carmoma/projects/pollen/results/nuclear1/binresults05",
  "/home/carmoma/projects/pollen/results/nuclear1/binresults10",
  "/home/carmoma/projects/pollen/results/nuclear1/binresults15",
  "/home/carmoma/projects/pollen/results/kraken_dust_mapq6/binresults01",
  "/home/carmoma/projects/pollen/results/kraken_dust_mapq6/binresults05",
  "/home/carmoma/projects/pollen/results/kraken_dust_mapq6/binresults10",
  "/home/carmoma/projects/pollen/results/kraken_dust_mapq6/binresults15",
  "/home/carmoma/projects/pollen/results/kraken_dust_noSamFilter/binresults01",
  "/home/carmoma/projects/pollen/results/kraken_dust_noSamFilter/binresults05",
  "/home/carmoma/projects/pollen/results/kraken_dust_noSamFilter/binresults10",
  "/home/carmoma/projects/pollen/results/kraken_dust_noSamFilter/binresults15",
  "/home/carmoma/projects/pollen/results/dustonly_score15/binresults01",
  "/home/carmoma/projects/pollen/results/dustonly_score15/binresults05",
  "/home/carmoma/projects/pollen/results/dustonly_score15/binresults10",
  "/home/carmoma/projects/pollen/results/dustonly_score15/binresults15",
  "/home/carmoma/projects/pollen/results/dustonly1/binresults01",
  "/home/carmoma/projects/pollen/results/dustonly1/binresults05",
  "/home/carmoma/projects/pollen/results/dustonly1/binresults10",
  "/home/carmoma/projects/pollen/results/dustonly1/binresults15",
  "/home/carmoma/projects/pollen/results/kraken_dust_1/binresults01",
  "/home/carmoma/projects/pollen/results/kraken_dust_1/binresults05",
  "/home/carmoma/projects/pollen/results/kraken_dust_1/binresults10",
  "/home/carmoma/projects/pollen/results/kraken_dust_1/binresults15",
  "/home/carmoma/projects/pollen/results/mock_kraken1/binresults01",
  "/home/carmoma/projects/pollen/results/mock_kraken1/binresults05",
  "/home/carmoma/projects/pollen/results/mock_kraken1/binresults10",
  "/home/carmoma/projects/pollen/results/mock_kraken1/binresults15",
  "/home/carmoma/projects/pollen/results/mock_trimgalore1/binresults01",
  "/home/carmoma/projects/pollen/results/mock_trimgalore1/binresults05",
  "/home/carmoma/projects/pollen/results/mock_trimgalore1/binresults10",
  "/home/carmoma/projects/pollen/results/mock_trimgalore1/binresults15",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_samfilt/diff_thresholds/binresults10",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_samfilt/diff_thresholds/binresults01",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_samfilt/diff_thresholds/binresults05",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_samfilt/diff_thresholds/binresults15",
  "/home/carmoma/projects/pollen/results/mock_nofilt1/all_results/binresults01",
  "/home/carmoma/projects/pollen/results/mock_nofilt1/all_results/binresults05",
  "/home/carmoma/projects/pollen/results/mock_nofilt1/all_results/binresults10",
  "/home/carmoma/projects/pollen/results/mock_nofilt1/all_results/binresults15",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_nosamfilt/diff_thresholds/binresults01",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_nosamfilt/diff_thresholds/binresults05",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_nosamfilt/diff_thresholds/binresults10",
  "/home/carmoma/projects/pollen/results/mock_ontfilt_nosamfilt/diff_thresholdsc"
             )
dirlist <- dirlist[1:4]
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
  res4 <- getCorrelations(res4, mattemp, mock_species, mock_species2)
  res4$condition <- resultsdir 
  write.table(res4, "results_with_correlation.csv", sep="\t", quote=F, row.names = F)
  all_res4 <- rbind(all_res4, res4)
  
  res2write <- as.data.frame(round(mattemp[["mat"]]*100, 3))
  write.table(res2write, "results_matrix_as_percent.csv", sep="\t", row.names = T, quote=F)
  write.table(as.data.frame(mattemp[["mat"]]), "results_matrix.csv", sep="\t", row.names = T, quote=F)

  res2write_trim <- res2write[, names(res2write)%in% spnames$Species[spnames$present_mock != "Other (false positive)"]]
  write.table(res2write_trim, "results_matrix_presentonly.csv", sep="\t", row.names = T, quote=F)

  makeAllHeatmaps(mattemp)
  getDotPlot(res4, mattemp)
  makePvclust(mattemp[["mat"]])

}

setwd("/home/carmoma//projects/pollen/results/")
write.table(all_res4, "230102_allresults_1percent.csv", sep="\t", quote=F, row.names = F)

all_prev <- read.table("230102_allresults_1percent.csv", sep="\t", stringsAsFactors=F, head=T)
#all_res4 <- rbind(all_prev, all_res4)

auxsum <- all_res4 %>% 
  group_by(condition) %>% 
  select(correlation, 
        correlation_present, 
        total_negative,
        mean_diff, mean_diff_positive,
        mean_diff_negative, 
        max_diff_positive, 
        max_diff_negative, 
        TP, TN, FP, FN, 
        specificity, sensitivity, PPV, NPV
    ) %>% 
  summarise_all(mean)
write.table(auxsum, "230115_allsummaryPlastids_1percent.csv", sep="\t", quote=F, row.names = F)
view(auxsum)
