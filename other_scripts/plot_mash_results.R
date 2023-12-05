library(tidyverse)
library(pvclust)
library(pheatmap)
library(dendextend)
library(khroma)

mock_species <- read_tsv("~/projects/pollen/METADATA/mock_composition.tsv") %>%
  select(-c(mix_11, mix_12, mix_13))
mock_species2 <- read_tsv("~/projects/pollen/METADATA/mock_composition2.tsv")
spnames <- read_tsv("~/projects/pollen/METADATA/species_name.tsv") %>% mutate(EASI_ID = as.character(EASI_ID))
barcodes <- read_csv("~/projects/pollen/METADATA/barcodes.tsv")
mocksamplenum <- read_tsv("~/projects/pollen/METADATA/easi_samplenum.tsv")
barcodes$PROJ <- gsub("_summary.csv", "",barcodes$PROJ) %>% gsub("_", "-", .)
barcodes$PROJ2 <- paste(barcodes$PROJ, barcodes$YY, sep="-")
mocksamplenum$rep <- gsub("[0-9]+", "", mocksamplenum$Mock_commiunity)
mocksamplenum$mix <- paste("mix_", gsub("[A-Z]", "", mocksamplenum$Mock_commiunity), sep="")
mocksamplenum$proj2 <- barcodes$PROJ2[match(mocksamplenum$EASI_ID, barcodes$ASSAMPLE)]

revmetmat <- read.table("/home/carmoma/projects/pollen/results/kraken_dust_1/binresults05/results_matrix.csv", 
                        head=T, row.names = 1, sep="\t") %>% as.matrix() %>% t
aggr <- revmetmat[grepl("Agrostis", rownames(revmetmat)), ]
matnoagr <- revmetmat[!rownames(revmetmat) %in% rownames(aggr), ]

aggr <- apply(aggr, MAR=2, sum) %>% t
rownames(aggr) <- "Agrostis"
revmetmat <- rbind(as.data.frame(aggr),as.data.frame(matnoagr)) %>% as.matrix()


pheatmap(revmetmat)
trim_barcode <- function(x){
  y <- strsplit(x, "/")[[1]]
  y <-y[length(y)]
  y <- strsplit(y, "\\.")[[1]][1]
  return(y)
}
trim_ill <- function(x){
  y <- strsplit(x, "/")[[1]]
  y <-y[length(y)]
  y <- strsplit(y, "_")[[1]][1] %>% gsub("AS0", "", .)
  return(y) 
}

prepare_input <- function(df, f1, f2, f1a, f1b, f2a, f2b){
  names(df) <- c("file1","file2", "dist", "pval", "matching_hases")
  df <- df %>% mutate(
    sample1 = sapply(file1, f1),
    sample2 = sapply(file2, f2),
    sample1_name = f1b[match(sample1, f1a)],
    sample2_name = f2b[match(sample2, f2a)]
  )
  }

makePvclust <- function(mat, distmet, clustmet){
  pvc <- pvclust(mat, method.dist = distmet, method.hclust = clustmet)
  pvname <- gsub("clmet", clustmet, "pvclust_clmet_distmet_pct15.pdf")
  pvname <- gsub("distmet", distmet, pvname)
  pdf(pvname)
  plot(pvc)
  dev.off()
  return(pvc)
}

buildMatrix <- function(df){
  mat <- df %>% 
    select(sample1_name, sample2_name, dist) %>% 
    spread(sample1_name, dist) 
  rownames(mat) <- mat$sample2_name
  mat <- as.matrix(mat[,-1])
  mat <- mat[order(rownames(mat)), order(colnames(mat))]
}

setwd("/home/carmoma/projects/pollen/results/mash_1")

ont2illf <- "on2ill.dist"
ont2ontf <- "on2on.dist"
ill2illf <- "ill2ill.dist"

ont2ill <- read.table(ont2illf, head=F, sep="\t", stringsAsFactors = F) %>% 
  prepare_input(trim_barcode, trim_ill, 
                mocksamplenum$proj2, mocksamplenum$Mock_commiunity,
                spnames$EASI_ID, spnames$Species)


ont2ont <- read.table(ont2ontf, head=F, sep="\t", stringsAsFactors = F)%>% 
  prepare_input(trim_barcode, trim_barcode, 
                mocksamplenum$proj2, mocksamplenum$Mock_commiunity,
                mocksamplenum$proj2, mocksamplenum$Mock_commiunity)


ill2ill <- read.table(ill2illf, head=F, sep="\t", stringsAsFactors = F) %>% 
  prepare_input(trim_ill, trim_ill, 
                spnames$EASI_ID, spnames$Species,
                spnames$EASI_ID, spnames$Species)




clustmethods <- c("ward.D", "ward.D2", 
                 "single", "complete", 
                 "average", "mcquitty", 
                 "median", "centroid")

distmethod <- "canberra" #canberra


#plot dist m2i
mo2i <- buildMatrix(ont2ill)
mo2inorm <- (1-mo2i) %>% apply(MAR=2, function(x)x/sum(x))

spsitalic <-sapply(
  rownames(mo2inorm),
  function(x) bquote(italic(.(x))))

p <- pheatmap(mo2inorm, 
              filename = gsub("dist" , distmethod, "heatmap_o2i_dist.png"), 
              labels_row = as.expression(spsitalic),
              angle_col=315,
              width = 12, height = 8)

mo2inormdist <- dist(mo2inorm %>% t, method = distmethod)
p <- pheatmap(mo2inormdist, 
              filename = gsub("dist" , distmethod, "heatmap_o2i_dist_dmap.png"), 
              labels_row = labels(mo2inormdist),
              labels_col = labels(mo2inormdist),
              width =8, height = 8)

dev.off()
pdfname <- gsub("dist" , distmethod, "hclust_o2i_dist.pdf")
pdf(pdfname)
hclust_o2i <- list()
for(distm in clustmethods){
  hclust_o2i[[distm]] <- hclust(mo2inormdist, method=distm)
  plot(hclust_o2i[[distm]], main=distm)
}
dev.off()

#Plot dist m2m
mo2o <- buildMatrix(ont2ont)
md <- as.dist(mo2o)
p <- pheatmap(md , 
              filename = "heatmap_o2o_dist.png",
              labels_row = labels(md),
              labels_col = labels(md),
              width =8, height = 8)


pdfname <- "hclust_o2o_dist.pdf"
pdf(pdfname)
hclust_o2o <- list()
for(distm in clustmethods){
  hclust_o2o[[distm]] <- hclust(md, method=distm)
  plot(hclust_o2o[[distm]], main=distm)
}
dev.off()

#plot dist i2i
mi2i <- buildMatrix(ill2ill)
mdi <- as.dist(mi2i)
hclust_i2i <- list()
pdfname <- "hclust_i2i_dist.pdf"
pdf(pdfname)
for(distm in clustmethods){
  hclust_i2i[[distm]] <- hclust(mdi, method=distm)
  plot(hclust_i2i[[distm]], main=distm)
}
dev.off()

#hclust original mixes
omat <- mock_species %>% filter(Spp != "Poa annua") %>% as.data.frame
rownames(omat) <- omat$Spp
names(omat) <- gsub("_", " ", names(omat))
omat <- omat %>% select(-Spp) %>% as.matrix
omatd <- dist(omat %>% t, method=distmethod)

hclust_theo <- list()
pdfname <- gsub("dist" , dismethod2, "hclust_theoreticalMixes_dist.pdf")
pdf(pdfname)
for(distm in clustmethods){
  hclust_theo[[distm]] <- hclust(omatd, method=distm)
  plot(hclust_theo[[distm]], main=distm)
}
dev.off()

dendtheo <- as.dendrogram(hclust_theo[["ward.D2"]])
plot(dendtheo)

#pdf("hclust_theoretical_canberra_wardd2.pdf", width = 8, height = 5)
png("hclust_theoretical_canberra_wardd2.png", width = 656, height = 410)
dendtheo %>%
  #set("nodes_pch", 19) %>% 
  set("leaves_pch", 19) %>% 
  set("nodes_cex", 2) %>% 
  set("leaves_cex", 2) %>% 
  set("leaves_col", "firebrick4") %>% 
  set("branches_lwd", 3) %>% 
  set("branches_col", "firebrick4") %>% 
  set("labels_cex", 2) %>% 
  set("labels_col", as.numeric(gsub("mix ", "", labels(dendtheo)))) %>% 
  plot(main = "Mock Mixes") 
dev.off()


#tree revmet
dismethod2 <- "canberra"
#BEST METHOD FOR REVMET DATA: CANBERRA DIST + AVERAGE, FOLLOWED BY WARD.D2
revmetdist <- dist(revmetmat %>% t, method = dismethod2)
p <- pheatmap(revmetdist, 
              filename = gsub("dist" , dismethod2, "heatmap_rvmet_dist_dmap_percent5.png"), 
              labels_row = labels(revmetdist),
              labels_col = labels(revmetdist),
              width =8, height = 8)

dev.off()
pdfname <- gsub("dist" , dismethod2, "hclust_revmet_dist_percent5.pdf")
pdf(pdfname)
hclust_revmet <- list()
for(distm in clustmethods){
  hclust_revmet[[distm]] <- hclust(revmetdist, method=distm)
  plot(hclust_revmet[[distm]], main=distm)
}
dev.off()

pvc <- makePvclust(revmetmat, dismethod2, "ward.D2")

pvcdend <- as.dendrogram(pvc)
pdfname <- gsub("dist" , dismethod2, "pvclust_revmet_wardd.2_dist_percent.png")
png(pdfname, width = 1000, height = 410)
#pdf(pdfname, width = 12, height = 5)
pvcdend %>%
  pvclust_show_signif_gradient(pvc) %>%
  set("labels_col", as.numeric(gsub("[ABC]", "", labels(pvcdend), perl=T))) %>%
  pvclust_show_signif(pvc) %>%
  set("leaves_pch", 19) %>% 
  set("nodes_cex", 2) %>% 
  set("leaves_cex", 2) %>% 
  set("leaves_col", "firebrick4") %>% 
  set("branches_lwd", 2) %>% 
  set("branches_col", "firebrick4") %>% 
  set("labels_cex", 1.5) %>%  
  plot(main = "RevMet results")
pvc %>% text
pvc %>% pvrect(alpha=0.8)
dev.off()

# compare tree m2m and m2i
dend1 <- as.dendrogram(hclust_o2o[["ward.D2"]])
dend2 <- as.dendrogram (hclust_o2i[["ward.D2"]])

dend_list <- dendlist("Mash ONT-ONT"=dend1, "Mash ONT-skims"=dend2)

png("MashOntOnt_vs_MashOntIll_canberra_wardd2.png", height = 800, width = 600)
pdf("MashOntOnt_vs_MashOntIll_canberra_wardd2.pdf", height = 8, width = 6)

dend_list %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = TRUE, # Turn-off dashed lines
    common_subtrees_color_lines = TRUE, # Turn-off line colors
    common_subtrees_color_branches = TRUE # Color common branches 
  )
dev.off()
cor.dendlist(dend_list, method = "cophenetic") #0.06

# compare tree m2m and revmet
#dend1 <- as.dendrogram (hclust_revmet[["ward.D2"]]) 
dend1 <- as.dendrogram(pvc)
dend2 <- as.dendrogram (hclust_o2o[["ward.D2"]])

dend_list <- dendlist("RevMet pipeline" = dend1, "Mash ONT-ONT" = dend2)

png("RevMet_vs_MashOntOnt_canberra_wardd2_percent5.png", height = 800, width = 600)
pdf("RevMet_vs_MashOntOnt_canberra_wardd2_percent5.pdf", height = 8, width = 6)
dend_list %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = TRUE, # Turn-off dashed lines
    common_subtrees_color_lines = TRUE, # Turn-off line colors
    common_subtrees_color_branches = TRUE # Color common branches 
  )
dev.off()
cor.dendlist(dend_list, method = "cophenetic")
cor.dendlist(dend_list, method = "baker") # 0.90
cor_cophenetic(dend1, dend2) #0.797, 0.778 with percent 15


# compare tree o2i and revmet
#dend1 <- as.dendrogram (hclust_revmet[["average"]]) 
dend1 <- as.dendrogram(pvc)
dend2 <- as.dendrogram (hclust_o2i[["ward.D2"]])

dend_list <- dendlist("RevMet pipeline"= dend1, "Mash ONT-skims"= dend2)

png("RevMet_vs_MashOntIll_canberra_wardd2.png", height = 800, width = 600)
pdf("RevMet_vs_MashOntIll_canberra_wardd2.pdf", height = 8, width = 6)

dend_list %>%
  untangle(method = "step1side") %>% # Find the best alignment layout
  tanglegram(
    highlight_distinct_edges = TRUE, # Turn-off dashed lines
    common_subtrees_color_lines = TRUE, # Turn-off line colors
    common_subtrees_color_branches = TRUE # Color common branches 
  )
dev.off()
cor.dendlist(dend_list, method = "cophenetic") #0.009

#Compare distance matrices with dot plot - o2o

#Compare distance matrices with dot plot - o2i





