
library(tidyverse)
library(RColorBrewer)

setwd("~/projects/pollen/results/mock_ontfilt_samfilt/alignments_filt")

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
  theme(axis.text.x = element_text(vjust = 1, 
                                   hjust = 1, size = 14, 
                                   colour = "black", angle = 0, 
                                   face = "bold"))+
  theme(axis.title.x = element_text(vjust = 1, hjust = 0.5, 
                                    size = 16, colour = "black", 
                                    angle = 0, face = "bold")) +
  theme(axis.title.y= element_text(vjust = 1, hjust = 0.5, 
                                   size = 16, colour = "black", 
                                   angle = 90, face = "bold"))

f <- list.files();f
 a <- read.table("FAR74611-1-NB01.trim_AS0221_221_EASI_48_7157AF_HFYTVDRXY_1_256UDI-idt-UMI.filt.bam.mapq",
                 head=F)
 names(a)[1] <- "mapq"
 b <- read.table("FAR74611-1-NB01.trim_AS0224_224_EASI_48_7160AF_HFWKLDRXY_1_198UDI-idt-UMI.filt.bam.mapq",
                 head=F)
 
 names(b)[1] <- "mapq"
a$alignment <- "Schedonorus pratensis (absent)"
b$alignment <- "Lolium perenne (present)"
df <- rbind(b, a) 

g1 <- ggplot(df, aes(x=mapq, col=alignment, linetype=alignment))+
  geom_density(bw=1, size=1.5) +
  mytheme+ 
  scale_fill_brewer(palette = "Dark2") +
  theme(legend.position = c(0.8, 0.93))
 g1
 
 ggsave("density_FAR74611_NB01_mapq_s.pdf", g1, width=8, height = 5)
 
 
 
 
 