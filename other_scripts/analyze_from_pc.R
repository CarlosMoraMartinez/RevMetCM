library(tidyverse)
library(pheatmap)
library(pvclust)

#setwd("/home/carmoma/projects/pollen/downloaded_bam/tmp/merge")

bad_samples <- c("FAR74611-1-NB02", "FAR74611-1-NB03", "FAR76967-1-NB01")
mock_species <- read_tsv("~/projects/pollen/METADATA/mock_composition.tsv")
spnames <- read_tsv("~/projects/pollen/METADATA/species_name.tsv")

setwd("/home/carmoma/projects/pollen/results/mock2_v2/mock2_v2_all")
PERC_LIM <- 15
MIN_PERC_SP <- 0.01

f <- list.files(pattern="_all.csv")

res <- data.frame()
for(ont in f){
  cat(ont, "\n")
  a <- read.table(ont, sep="\t", head=T, stringsAsFactors = F)
  names(a) <- c("species", "read_id", "pct_covered")
  a <- a %>% group_by(read_id) %>% 
    top_n(n = 1, wt = pct_covered) %>% 
    mutate(ont_sample = gsub("_all.pc", "", ont)) %>% 
    ungroup
  b <- table(a$read_id)
  b <- b[b > 1]
  a <- a[!a$read_id %in% names(b), ]
  res <- rbind(res, a)
}

write.table(res, "ALL_READS_BINNED.csv", sep="\t", row.names = F, quote=F)
#res2 <- res %>% 
#  group_by(ont_sample, read_id) %>% 
#  top_n(n = 1, wt = pct_covered)

#REMOVE READS WITH LESS THAN 15% OF THEIR LENGTH COVERED

res$binned_speces <- mapply(res$species, res$pct_covered, FUN=function(sp, pc) ifelse(pc < PERC_LIM, "Unassigned", sp))
#res$binned_speces <- res$species
#CALCULATE PERCENTAGE OF SPECIES
res3 <- res %>% group_by(ont_sample) %>% 
  mutate(n_ont = n()) %>% 
  filter(binned_speces != "Unassigned") %>% 
  group_by(ont_sample, binned_speces) %>% 
  summarise(n_species = n(),
            n_ont = unique(n_ont)) %>% 
  mutate(species = binned_speces,
         pct_species = n_species / n_ont) %>% 
  mutate(pct_species = pct_species/sum(pct_species),
         pct_species_trimmed = ifelse(pct_species < MIN_PERC_SP, 0, pct_species),
         species_ref = species,
         species_num = as.numeric(gsub("AS0", "", species)),
         species = spnames$Species[match(species_num, spnames$EASI_ID)]
         )

head(res3)
write.table(res3, "res3_partial_summary_mock_15_1_norm.csv", sep="\t", row.names = F, quote = F)
res3<- read.table("res3_partial_summary_mock.csv", head=T, sep="\t")

#Only species in mock samples
#res3 <- res3 %>% filter(species %in% spnames$Species[spnames$present_mock != "absent"])

res4 <- res3 %>% 
  select(ont_sample, species, pct_species_trimmed) %>% 
  spread(species, pct_species_trimmed) %>% data.frame

mat <- as.matrix(res4[, -1])
mat[is.na(mat)] <- 0


rownames(mat) <- res4$ont_sample
colnames(mat) <- gsub("\\.", " ", colnames(mat)) %>% gsub(" $", ".", .)

apply(mat, MAR=1, sum)
#mat2 <- mat[, apply(mat, MAR=2, FUN=function(x)any(x>0))]
mat2 <- apply(mat, MAR=1, FUN=function(x){x/sum(x)})
mat2 <- t(mat2)
apply(mat2, MAR=1, sum)
mat <- mat2

ann <- data.frame(sapply(rownames(mat), FUN=function(x)strsplit(x, "-")[[1]][1]))
names(ann)[1] <- "mock group"
annsp <- data.frame(presence = spnames$present_mock[match(colnames(mat), spnames$Species) ])
names(annsp)[1] <- "species presence"
rownames(annsp) <- colnames(mat)
#
spsitalic <-sapply(
  colnames(mat),
  function(x) bquote(italic(.(x))))

p <- pheatmap(t(mat), 
         annotation_col = ann, 
         annotation_row = annsp,  
         filename = "partial_mock_samples_all_15_1_norm_b.png", 
         labels_row = as.expression(spsitalic),
         angle_col=315,
         width = 12, height = 8)

p <- pheatmap(t(mat), 
              annotation_col = ann, 
              annotation_row = annsp,
              filename = "partial_mock_samples_all_15_1_norm_b_noclust.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(3, nrow(mat)-1, by=3),
              angle_col=315,
              width = 12, height = 8)

#Only species of interest
matcut <- mat[ , colnames(mat) %in% spnames$Species[spnames$present_mock != "absent"]]

pheatmap(t(matcut), 
         annotation_col = ann, 
         annotation_row = annsp,
         filename = "partial_mock_samples_all_15_1_norm_b_selectedspecies.png", 
         labels_row = as.expression(spsitalic),
         angle_col=315,
         width = 12, height = 8)

p <- pheatmap(t(matcut), 
              annotation_col = ann, 
              annotation_row = annsp,
              filename = "partial_mock_samples_all_15_1_norm_b_selectedspecies_noclust.png", 
              labels_row = as.expression(spsitalic),
              cluster_cols = F,
              gaps_col = seq(3, nrow(mat)-1, by=3),
              angle_col=315,
              width = 12, height = 8)

res2write <- as.data.frame(round(mat*100, 3))
write.table(res2write, "Mock_results_1_cov15_sp1.csv", sep="\t", row.names = T, quote=F)

res2write_trim <- res2write[, names(res2write)%in% spnames$Species[spnames$present_mock != "absent"]]
write.table(res2write_trim, "Mock_results_1_cov15_sp1_trim.csv", sep="\t", row.names = T, quote=F)


apply(round(mat*100, 3), MAR=1, sum)

my_gtable = p$gtable

my_gtable$grobs[[3]]$gp=gpar(col="#ffff00", fontsize=20)# assuming that the xlabels are in the third grob
my_gtable$grobs[[4]]$gp=gpar(col="#ffffff", fontsize=20)# assuming that the ylabels are in the fourth grob
my_gtable$grobs[[1]]$gp=gpar(col="#ffffff", lwd=2) # change the color of the dendrogram and set the linewidth to 2
my_gtable$grobs[[5]]$gp=gpar(col="#ffffff", fontsize="20", just="center") 

