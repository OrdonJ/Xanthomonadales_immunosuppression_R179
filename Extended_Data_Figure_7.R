# Extended Data Figure 7

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)
library(rstatix)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Figure 7a

# read files
RL_supernatant <- read.delim("Ext_Data_Fig_7a.txt", stringsAsFactors = F)

# add combined column
RL_supernatant$strain_treatment <- paste0(RL_supernatant$strain,"_", RL_supernatant$treatment)

# define colors
color <- c("none" = "NA", "wt" = "#ffc300", "wt_HK" ="#FF6200","dssAB" = "#386641")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"


# Kruskal-Wallis test
kw_supernatant <- kruskal.test(strain_treatment ~ root_length, data=RL_supernatant)
print(kw_supernatant)

dn_supernatant <- dunnTest (root_length ~ strain_treatment, data=RL_supernatant, method=p.adjust.method)
dn_supernatant <- dn_supernatant$res
print(dn_supernatant)

l_supernatant <- cldList(comparison=dn_supernatant$Comparison, p.value=dn_supernatant$P.adj, threshold=0.01, remove.zero = FALSE)
print(l_supernatant)

# plot
ggplot((data = RL_supernatant), aes(x=strain_treatment, y=root_length, fill=strain))+
  geom_boxplot(aes(fill = strain, alpha = 0.4), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape=replicate), size= 0.6, width = 0.35, alpha =0.4)+
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  scale_x_discrete(limits=c("none_mock", "wt_mock", "dssAB_mock", "wt_HK_mock", 
                            "none_flg22",  "wt_flg22", "dssAB_flg22", "wt_HK_flg22")) +
  scale_fill_manual(values = color)+
  annotate("text",
           x = l_supernatant$Group,
           label = l_supernatant$Letter,
           y=7.2,
           size= 2.5)+
  labs(title="supernatant") +
  ggsave(paste("Root_length_supernatant_R179wt_dssAB.pdf", sep=""), width=3.3, height=4)

################################################################################################################################
# Extended Data Figure 7b

# read files
flg22_HKR179 <- read.delim("Ext_Data_Fig_7b.txt", stringsAsFactors = F)

# define colors
color <- c("none" = "NA", "wt" = "#ffc300", "HK-wt" ="#FF6200","dssAB" = "#386641")

# plot
ggplot((data=flg22_HKR179), aes(x=Bacteria, y=Total_area.MS1)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("axenic", "wt", "HK-wt")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("flg22_quantification_HK-R179_1h_XVM2_biggest_precursor_rel.pdf", sep=""), width=3, height=4)

################################################################################################################################
# Extended Data Figure 7c

# read files
cleaved_Atpep1 <- read.delim("Ext_Data_Fig_7c.txt", stringsAsFactors = F)

# plot
ggplot((data=cleaved_Atpep1), aes(x=Bacteria, y=Total_area.MS1)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("axenic", "wt", "dssAB")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("cleaved_Atpep1_quantification_XVM2_biggest_precursor_rel.pdf", sep=""), width=3, height=4)


################################################################################################################################
# Extended Data Figure 7d

# read files
XVM2_TSB <- read.delim("Ext_Data_Fig_7d.txt", stringsAsFactors = F)

# add combined column
XVM2_TSB$peptide_bacteria <- paste0(XVM2_TSB$Protein.Name,"_", XVM2_TSB$Bacteria)

# plot
ggplot((data=XVM2_TSB), aes(x=peptide_bacteria, y=percent_rel_axenic)) + 
  geom_jitter(aes(color = Bacteria, alpha = 0.8), size= 2, width = 0.25) + 
  stat_summary(fun=mean, geom="point", color="black", size = 0.3) +
  scale_color_manual(values = color)+
  scale_x_discrete(limits=c("flg22_axenic", "flg22_wt", "flg22_dssAB", 
                            "elf18_axenic", "elf18_wt", "elf18_dssAB",
                            "AtPep1_axenic", "AtPep1_wt", "AtPep1_dssAB",
                            "dRGF1_axenic", "dRGF1_wt", "dRGF1_dssAB",
                            "PSY1_axenic", "PSY1_wt", "PSY1_dssAB")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size = 11, angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black", size = 11)) +
  ggsave(paste("Peptide_quantification_1h_XVM2_TSB_biggest_precursor_rel.pdf", sep=""), width=5, height=4)
