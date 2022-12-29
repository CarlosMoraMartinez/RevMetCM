library(tidyverse)
library(pheatmap)
library(pvclust)

#setwd("/home/carmoma/projects/pollen/downloaded_bam/tmp/merge")

bad_samples <- c("FAR74611-1-NB02", "FAR74611-1-NB03", "FAR76967-1-NB01")
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

#setwd("~/projects/pollen/results/mock_nofilt1/mock_nofilt_all")
#setwd("/home/carmoma/projects/pollen/results/mock_ontfilt_samfilt/")
#setwd("/home/carmoma/projects/pollen/results/mock_ontfilt_nosamfilt/coverage")
#setwd("/home/carmoma/projects/pollen/results/mock_nofilt1/mock_nofilt_all")
#setwd("~/projects/pollen/results/mock_ontfilt_nosamfilt/diff_thresholds/binresults01")
setwd("/home/carmoma//projects/pollen/results/mock_ontfilt_samfilt/diff_thresholds/binresults10")

PERC_LIM <- 10
MIN_PERC_SP <- 0.01

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
#res <- read.table("ALL_READS_BINNED.csv", head=T, sep="\t", stringsAsFactors=F)
#res2 <- res %>% 
#  group_by(ont_sample, read_id) %>% 
#  top_n(n = 1, wt = pct_covered)

#REMOVE READS WITH LESS THAN 15% OF THEIR LENGTH COVERED

res$binned_species2 <- mapply(res$ILLUMINA, res$coverage, FUN=function(sp, pc) ifelse(pc < PERC_LIM, "Unassigned", sp))
#res$binned_speces <- res$species
#CALCULATE PERCENTAGE OF SPECIES
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

head(res3)
write.table(res3, "percents_10_1_norm.csv", sep="\t", row.names = F, quote = F)
#res3<- read.table("res3_partial_summary_mock.csv", head=T, sep="\t")

#Only species in mock samples
#res3 <- res3 %>% filter(species %in% spnames$Species[spnames$present_mock != "absent"])

res4 <- res3 %>% 
  select(mock_sample, mock_mix, dilution, species, pct_species_trimmed) %>% 
  spread(species, pct_species_trimmed)  %>% 
  data.frame %>% select(-ONT)

mat <- res4 %>% select_if(is.numeric) %>% as.matrix
mat[is.na(mat)] <- 0
rownames(mat) <- res4$mock_sample

colnames(mat) <- gsub("\\.", " ", colnames(mat)) %>% gsub(" $", ".", .)

#apply(mat, MAR=1, sum)
#mat2 <- mat[, apply(mat, MAR=2, FUN=function(x)any(x>0))]
#mat2 <- apply(mat, MAR=1, FUN=function(x){x/sum(x)})
#mat2 <- t(mat2)
#apply(mat2, MAR=1, sum)
#mat <- mat2

ann <- data.frame(mix = as.factor(gsub("mix_", "", res4$mock_mix)), 
                  dilution = res4$dilution)
names(ann)[1] <- "mock group"
rownames(ann) <- rownames(mat)
#annsp <- data.frame(presence = spnames$present_mock[match(colnames(mat), spnames$Species) ])
#names(annsp)[1] <- "species presence"
#rownames(annsp) <- colnames(mat)
#
ord <- order(rownames(mat))
ord <- ord[c(4:length(ord), 1:3)]
mat <- mat[ord, ]
ann <- ann[ord, ]

#pdf("mock_ontfilt1.pdf")
#x <-hclust(t(dist(mat, method = "canberra")))
pvc <- pvclust(t(mat), method.dist = "canberra")
pdf("pvclust_mock.pdf")
plot(pvc)
dev.off()

dev.off()

spsitalic <-sapply(
  colnames(mat),
  function(x) bquote(italic(.(x))))

#Crear matriz adyacente con los datos mock
matextra <- matrix(numeric(length(unique(res4$mock_mix))*ncol(mat)), nrow=ncol(mat))
rownames(matextra) <- colnames(mat)
for (n in mock_species2$Spp){
  matextra[n, ] <- unlist(mock_species2[mock_species2$Spp == n, names(mock_species2) %in% res4$mock_mix])
}
colnames(matextra) <- as.factor(as.character(1:10))
matextra <- as.data.frame(matextra)


p <- pheatmap(t(mat), 
         annotation_col = ann, 
#         annotation_row = annsp,  
         filename = "partial_mock_samples_all_15_1_norm_b.png", 
         labels_row = as.expression(spsitalic),
         angle_col=315,
         width = 12, height = 8)
dev.off()
p <- pheatmap(t(mat), 
              annotation_col = ann, 
#              annotation_row = matextra,
              filename = "partial_mock_samples_all_10_1_norm_b_noclust.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(3, nrow(mat)-1, by=3),
              angle_col=315,
              width = 12, height = 8)

dev.off()
p
#Only species of interest
matcut <- mat[ , colnames(mat) %in% mock_species2$Spp]

p <- pheatmap(t(matcut), 
         annotation_col = ann, 
        # annotation_row = annsp,
         filename = "partial_mock_samples_all_15_1_norm_b_selectedspecies.png", 
         labels_row = as.expression(spsitalic),
         angle_col=315,
         width = 12, height = 8)
dev.off()
p
p <- pheatmap(t(matcut), 
              annotation_col = ann, 
              #annotation_row = annsp,
              filename = "partial_mock_samples_all_15_1_norm_b_selectedspecies_noclust.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(3, nrow(mat)-1, by=3),
              angle_col=315,
              width = 12, height = 8)
dev.off()
p
res2write <- as.data.frame(round(mat*100, 3))
write.table(res2write, "Mock_results_1_cov15_sp1.csv", sep="\t", row.names = T, quote=F)

res2write_trim <- res2write[, names(res2write)%in% spnames$Species[spnames$present_mock != "absent"]]
write.table(res2write_trim, "Mock_results_1_cov15_sp1_trim.csv", sep="\t", row.names = T, quote=F)

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

genus_species <- sapply(rownames(mergedmat2), FUN=function(x){strsplit(x, " ")[[1]][1]})
ann_species <- data.frame(present = ifelse(rownames(mergedmat2) %in% splist, "species", 
                          ifelse(genus_species %in% genus, "genus", "absent")
                          ))


rownames(ann_species) <- rownames(mergedmat)

ann2 <- data.frame(dilution = rep(c("theoretical %", "dil A", "dil B", "dil C"), dim(mergedmat2)[2]/4))
rownames(ann2) <- colnames(mergedmat2)


p <- pheatmap(mergedmat2, 
              annotation_col = ann2, 
              annotation_row = ann_species,
              filename = "heatmap_with_mix.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(4, nrow(mergedmat2)-1, by=4),
              angle_col=315,
              width = 14, height = 8)
dev.off();p
# Correlations 
matextra3 <- matextra
colnames(matextra3) <- paste("mix_", colnames(matextra3), sep="")

for(g in colnames(matextra3)){
  ind = which(res4$mock_mix == g)
  for(i in ind){
  res4$correlation[i] <- cor(matextra3[, g], mat[res4$mock_sample[i], ])
  res4$correlation_present[i] <- cor(matextra3[mock_species2$Spp, g], mat[res4$mock_sample[i], mock_species2$Spp])
  }
}
write.table(res4, "results_with_correlation.csv", sep="\t", quote=F, row.names = F)
res_nofilt <- res4
mean(res4$correlation)

#0.5349616 ont filt + no sam filt 1% covered
#0.7569506  ont filt + no sam filt 5% covered
#0.7550329 ont filt + no sam filt 10% covered
#0.7517308 ont filt + no sam filt 15% covered

# 0.47798 filt ont filt +  sam filt 1% covered
# 0.7576867 filt ont filt +  sam filt 5% covered
# 0.7544688 filt ont filt +  sam filt 10% covered
#0.750 filt ont filt + sam filt 15% covered

#0.751728 no ont filt + no sam filt 15% covered

mean(res4$correlation_present)
# 0.631657  ont filt + no sam filt 1% covered
#0.6527326 ont filt + no sam filt 5% covered
#0.6490048 ont filt + no sam filt 10% covered
#0.645 ont filt + no sam filt 15% covered

# 0.6538847 ont filt +sam filt 1% covered
# 0.652368 ont filt +  sam filt 5% covered
#0.6466251 ont filt +  sam filt 10% covered
# 0.64 ont filt +sam filt 15% covered

# 0.6463454 no ont filt + no sam filt 15% covered

##Dendograms
library(dendextend)

dend <- hclust(dist(as.matrix(t(matextra))))
dend2 <- hclust(dist(as.matrix(t(mat2))))
dl <- dendlist(dend, dend2)


library(pvclust)
pvc <- pvclust(mat, method.dist = "canberra")
