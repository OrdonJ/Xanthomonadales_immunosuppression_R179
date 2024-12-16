# Figure 4

library("tidyverse")
library("rstatix")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Figure 4a

# read data files
peptides <- read.delim("Fig.4a.txt", stringsAsFactors = F)

# add combined column
peptides$peptide_bacteria <- paste0(peptides$Protein.Name,"_", peptides$Bacteria)

# define colors
color <- c("axenic" = "black", "wt" = "#ffc300", "dssAB" = "#386641")

idx <- peptides$Protein.Name %in% c("flg22", "elf18", "AtPep1") &
  TRUE
immune_peptides <- peptides[idx, ]

idx <- peptides$Protein.Name %in% c("dRGF1") &
  TRUE
dRGF1 <- peptides[idx, ]

idx <- peptides$Protein.Name %in% c("PSY1") &
  TRUE
PSY1 <- peptides[idx, ]

# plot
ggplot((data=immune_peptides), aes(x=peptide_bacteria, y=percent_rel_axenic)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("flg22_axenic", "flg22_wt", "flg22_dssAB", 
                            "elf18_axenic", "elf18_wt", "elf18_dssAB",
                            "AtPep1_axenic", "AtPep1_wt", "AtPep1_dssAB")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("Peptide_quantification_1h_XVM2_biggest_precursor_rel.pdf", sep=""), width=5, height=4)

# plot
ggplot((data=dRGF1), aes(x=peptide_bacteria, y=percent_rel_axenic)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("dRGF1_axenic", "dRGF1_wt", "dRGF1_dssAB")) +
  ylim(0, 450)+
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("Peptide_quantification_1h_XVM2_dRGF1_rel.pdf", sep=""), width=2, height=4)

# plot
ggplot((data=PSY1), aes(x=peptide_bacteria, y=percent_rel_axenic)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("PSY1_axenic", "PSY1_wt", "PSY1_dssAB")) +
  ylim(0, 150)+
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("Peptide_quantification_1h_XVM2_PSY1_rel.pdf", sep=""), width=2, height=4)

################################################################################################################################
# Figure 4b

MF79_peptides <- read.delim("Fig.4b.txt", stringsAsFactors = F)

# add combined column
MF79_peptides$peptide_bacteria <- paste0(MF79_peptides$Protein.Name,"_", MF79_peptides$Bacteria)

# define colors
color <- c("axenic" = "black", "MF79" = "#15587A", "MF79 gspD" = "#4F88DB")

######### Plot
ggplot((data=MF79_peptides), aes(x=peptide_bacteria, y=percent_rel_axenic)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.9), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("flg22_axenic", "flg22_MF79", "flg22_MF79 gspD", 
                            "elf18_axenic", "elf18_MF79", "elf18_MF79 gspD",
                            "AtPep1_axenic", "AtPep1_MF79", "AtPep1_MF79 gspD")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("Peptide_quantification_1h_XVM2_biggest_precursor_MF79.pdf", sep=""), width=5, height=4)

