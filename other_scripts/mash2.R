library(tidyverse)
library(pvclust)
library(pheatmap)
library(dendextend)
library(khroma)
library(gridExtra)


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

trim_barcode <- function(x){
  #y <- strsplit(x, "/")[[1]]
  #y <-y[length(y)]
  y <- strsplit(x, "\\.")[[1]][1]
  return(y)
}
trim_ill <- function(x){
  y <- strsplit(x, "/")[[1]]
  y <-y[length(y)]
  y <- strsplit(y, "_")[[1]][1] %>% gsub("AS0", "", .)
  return(y) 
}

prepare_input <- function(df, f1, f2, f1a, f1b, f2a, f2b){
  names(df) <- c("file1","file2", "dist", "matching_hases", "x", "pval")

    df$sample1 = sapply(df$file1, f1)
    df$sample2 = sapply(df$file2, f2)
    df$sample1_name = f1b[match(df$sample1, f1a)]
    df$sample2_name = f2b[match(df$sample2, f2a)]
    return(df)
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
  mat[is.na(mat)] <- 0
  mat <- mat[order(rownames(mat)), order(colnames(mat))]
}

##
setwd("/home/carmoma/projects/pollen/mash2")
files <- list.files(pattern="*.dist")

res <- data.frame()
for(f in files){
  aux <- read.table(f, sep="\t", stringsAsFactors = F)
  aux$V7 <- f
  aux <- aux[, c(7, 5, 1, 2, 3, 4)]
  res <- rbind(res, aux)
}

ALPHA <- 0.000000001
MIN_HASHES <- 3
ont2ill <- res %>% 
  prepare_input(trim_barcode, trim_ill, 
                mocksamplenum$proj2, mocksamplenum$Mock_commiunity,
                spnames$EASI_ID, spnames$Species)

ont2ill$similarity <- ont2ill$dist
ont2ill$dist <- gsub("/1000", "", ont2ill$matching_hases) %>% as.numeric
ont2ill$dist[ont2ill$pval > ALPHA ] <-0
ont2ill$dist[ont2ill$dist < MIN_HASHES ] <-0

ont2ill$species <- ont2ill$sample2_name
ont2ill$dilution <- gsub("[0-9]+", "", ont2ill$sample1_name, perl=T)
ont2ill$mock_mix <- gsub("[ABC]+", "", ont2ill$sample1_name, perl=T)
ont2ill$sample <- paste(ont2ill$mock_mix, ont2ill$dilution, sep="")

data <- ont2ill %>% group_by(mock_mix, dilution) %>% 
  mutate(proportion_found = dist/sum(dist)) %>% 
  select(dilution, mock_mix, sample, species, proportion_found, pval, similarity) %>% 
  as.data.frame()
sp = unique(data$species)
for(v in unique(data$sample)){
  spu <- data$species[data$sample == v]
  spnot <- sp[! sp %in% spu]
  if(length(spnot) > 0){
    dil <- gsub("[0-9]+", "", v)
    mix <- gsub("[ABC]", "", v)
    aux <- data.frame(dilution = dil, mock_mix = mix, 
                      sample = v, 
                    species = spnot, proportion_found = 0,
                    pval = 0, similarity = 0)
    data <- rbind(data, aux)
  }
}

data$type <- "retrieved"
mock_species_data <- data %>% filter(dilution == "A") %>% 
  mutate(proportion_found = 0,
         pval = 1, 
         similarity = 1,
         dilution = "Exp.",
         type = "theoretical",
         proportion_expected = proportion_found,
         `plant present` = "true positive")
mat <- mock_species2 %>% select(-Spp) %>% as.matrix()

mat[grepl("Poa", rownames(mat)), ] <- 0
rownames(mat) <- mock_species2$Spp
colnames(mat) <- gsub("mix_", "", colnames(mat))

rownames(mat) %in% mock_species_data$species
mock_species_data$proportion_found <- mapply(mock_species_data$species, 
                                             mock_species_data$mock_mix, 
                                             FUN=function(sp, mix){
                                               if(sp %in% rownames(mat)){
                                                 mat[sp, mix]
                                               }else{
                                                 return(0)
                                               }
                                               
                                               })

for(i in 1:nrow(data)){
    data[i, "proportion_expected"] <- mock_species_data[mock_species_data$mock_mix == data[i, "mock_mix"] & mock_species_data$species == data[i, "species"], "proportion_found"]
}
data$`plant present` <- ifelse(data$proportion_found > 0 & data$proportion_expected == 0, "false positive", "true positive")
resvert5 <- rbind(mock_species_data, data)

resvert5$dilution <- factor(resvert5$dilution, levels = c("Exp.", "A", "B", "C"))

plants <- unique(data$species)
not_in <- plants[! plants %in% mock_species2$Spp & ! grepl("Poa ", plants)]
resvert5$factor <- ifelse(resvert5$species %in% not_in, -1*resvert5$proportion_found, resvert5$proportion_found)

ordersp <- resvert5 %>% group_by(species) %>% summarise(media = mean(factor))
ordersp <- ordersp$species[order(ordersp$media)]
resvert5$species <- factor(resvert5$species, levels = ordersp)
resvert5$mock_mix <- paste("mix ", resvert5$mock_mix, sep="") %>% 
  factor(levels =paste("mix ", as.character(1:10), sep=""))

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
ggsave("g6_MASH_byGroupWithTheo_pminus9_3hash.pdf", g6, width = 13, height = 9)

## Corrrelation
## Read again with separated Agrostis
revmetmatsep <- read.table("/home/carmoma/projects/pollen/results/kraken_dust_1/binresults05/results_matrix.csv", 
                        head=T, row.names = 1, sep="\t") %>% as.matrix() %>% t
rownames(revmetmatsep) <- gsub("\\.", " ", rownames(revmetmatsep))
rownames(revmetmatsep) <- gsub(" L ", " L.", rownames(revmetmatsep))

resvert6 <- resvert5 %>% filter(dilution != "Exp.")
resvert6$revmet <- mapply(as.character(resvert6$species), as.character(resvert6$sample), FUN=function(sp, sample) revmetmatsep[sp, sample])
plot(resvert6$proportion_expected, resvert6$revmet)
resvert6$`plant present` <- ifelse(resvert6$proportion_expected > 0, "present", "absent")

g7 <- ggplot(resvert6, aes(y = proportion_found,x = proportion_expected, col=`plant present`)) +   
  geom_abline(slope = 1,linetype=4, col="grey", size=1.5)+
  geom_point()  +  
  #scale_size(range = c(0, 8))+             ## to tune the size of circles
  scale_colour_manual(values = c("firebrick3", "dodgerblue2")) +
  xlab("theoretical proportion")+
  ylab("MASH screen retrieved proportion")+
  xlim(0, 0.55)+
  ylim(0,1)+
  mytheme +
  theme(legend.position = "none")

g8 <- ggplot(resvert6, aes(y = revmet,x = proportion_expected, col=`plant present`)) + 
  geom_abline(slope = 1,linetype=4, col="grey", size=1.5)+
  geom_point()  +  
  #scale_size(range = c(0, 8))+             ## to tune the size of circles
  scale_colour_manual(values = c("firebrick3", "dodgerblue2")) +
  xlab("theoretical proportion")+
  ylab("RevMet retrieved proportion")+
  xlim(0, 0.55)+
  ylim(0,1)+
  mytheme +
  theme(legend.position = "none")

g9 <- ggplot(resvert6, aes(x = revmet,y = proportion_found, col=`plant present`)) + 
  geom_abline(slope = 1,linetype=4, col="grey", size=1.5)+
  geom_point()  +  
  #scale_size(range = c(0, 8))+             ## to tune the size of circles
  scale_colour_manual(values = c("firebrick3", "dodgerblue2")) +
  ylab("MASH screen retrieved proportion")+
  xlab("RevMet retrieved proportion")+
  xlim(0, 1)+
  ylim(0,1)+
  mytheme +
 theme(legend.position = "none")

gex <- grid.arrange(g7, g8, g9, ncol=3)
ggsave("Correlations2.pdf", gex, width = 16, height = 7)
cor(resvert6$proportion_found, resvert6$revmet) # 0.5931

#### # heatmap, dend, etc
mo2i <- buildMatrix(ont2ill)
mo2i[mo2i < 5] <- 0
mo2inorm <- mo2i %>% apply(MAR=2, function(x)x/sum(x))
pheatmap(mo2inorm)
hclust(mo2inorm %>% t %>% dist(method="euclidean")) %>% plot


aggr <- mo2inorm[grepl("Agrostis", rownames(mo2inorm)), ]
matnoagr <- mo2inorm[!rownames(mo2inorm) %in% rownames(aggr), ]
aggr <- apply(aggr, MAR=2, sum) %>% t
rownames(aggr) <- "Agrostis"
mo2inorm2 <- rbind(as.data.frame(aggr),as.data.frame(matnoagr)) %>% as.matrix()


mo2inorm2 <- mo2inorm2[rownames(revmetmat), ]
spsitalic <-sapply(
  rownames(mo2inorm),
  function(x) bquote(italic(.(x))))
mo2inorm <- mo2inorm[,c(4:30, 1:3)]

p <- pheatmap(mo2inorm, 
              #annotation_col = ann, 
              #annotation_row = ann_species,
              filename = "heatmap_bygroup_mash2.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(3, ncol(mat)-1, by=3),
              angle_col=315,
              width = 12, height = 8)

pvcmash <- pvclust(mo2inorm, method.dist = "canberra", method.hclust = "ward.D2")
pvc <- pvclust(mat, method.dist = "canberra", method.hclust = "ward.D2")


dend1 <- as.dendrogram(pvc)
dend2 <- as.dendrogram (pvcmash)

dend_list <- dendlist("RevMet pipeline" = dend1, "Mash ONT-skims" = dend2)

png("RevMet_vs_MashOntIll_canberra_wardd2_percent5_mash2.png", height = 800, width = 600)
pdf("RevMet_vs_MashOntIll_canberra_wardd2_percent5_mash2.pdf", height = 8, width = 6)
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
