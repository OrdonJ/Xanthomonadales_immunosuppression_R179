# Extended Data Figure 8

library("tidyverse")
library("EnvStats")
library("dunn.test")
library("rcompanion")
library("FSA")
library("rstatix")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Figure 8a

# load data
qPCR_R179_PER5 <- read.delim("Ext_Data_Fig_8a.txt", stringsAsFactors = F)

# add combined column
qPCR_R179_PER5$treatment_timepoint <- paste0(qPCR_R179_PER5$treatment,"_", qPCR_R179_PER5$timepoint)

# define colors
colors <- c(  "mock_1h" = "grey","R179_0.5_1h" = "#ffc300",  "R179_1_1h" = "#ffc300","R179_1.5_1h" = "#ffc300",
              "R179_0.5_HK_1h" = "#FF6200","R179_1_HK_1h" = "#FF6200","R179_1.5_HK_1h" = "#FF6200",
              "mock_5h" = "grey", "R179_0.5_5h" = "#ffc300",  "R179_1_5h" = "#ffc300","R179_1.5_5h" = "#ffc300",
              "R179_0.5_HK_5h" = "#FF6200","R179_1_HK_5h" = "#FF6200","R179_1.5_HK_5h" = "#FF6200")

## statistics

qPCR_R179_PER5 %>%
  group_by(timepoint) %>%
  dunn_test(data =., fold_change ~ treatment_timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter(., group1=="mock_1h") %>% 
  filter()->stats_1h

qPCR_R179_PER5 %>%
  group_by(timepoint) %>%
  dunn_test(data =., fold_change ~ treatment_timepoint) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter(., group1=="mock_5h") %>% 
  filter()->stats_5h

# plot
ggplot(qPCR_R179_PER5, aes(x = treatment_timepoint, y= fold_change)) +
  geom_bar(aes(fill = treatment_timepoint), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("mock_1h","R179_0.5_1h","R179_1_1h","R179_1.5_1h",
                            "R179_0.5_HK_1h","R179_1_HK_1h","R179_1.5_HK_1h", 
                            "mock_5h","R179_0.5_5h","R179_1_5h","R179_1.5_5h",
                            "R179_0.5_HK_5h","R179_1_HK_5h","R179_1.5_HK_5h")) +
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, size=0.8) +
  geom_jitter(aes(color = treatment_timepoint, shape = Replicate), alpha = 0.7, size = 1) +
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  ggsave(paste("qPCR_R179_acute_PER5_dosage.pdf", sep=""), width=3.5, height=4)

################################################################################################################################
# Extended Data Figure 8c

# load data
qPCR_R179_EFR <- read.delim("Ext_Data_Fig_8c.txt", stringsAsFactors = F)

# add combined column
qPCR_R179_EFR$treatment_timepoint <- paste0(qPCR_R179_EFR$treatment,"_", qPCR_R179_EFR$timepoint)

# define colors
color <- c("mock" = "black", "R179_1.5" = "#ffc300", "R179_1.5_HK" ="#FF6200", "flg22" = "#DBC7AF", "elf18" = "#a3887a")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_EFR <- kruskal.test(treatment ~ fold_change, data=qPCR_R179_EFR)
print(kw_EFR)

dn_EFR <- dunnTest (fold_change ~ treatment, data=qPCR_R179_EFR, method=p.adjust.method)
dn_EFR <- dn_EFR$res
print(dn_EFR)

l_EFR<- cldList(comparison=dn_EFR$Comparison, p.value=dn_EFR$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_EFR)

# plot
ggplot(qPCR_R179_EFR, aes(x = treatment, y= fold_change)) +
  geom_bar(aes(fill = treatment), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, size=0.8) +
  geom_jitter(aes(color = treatment, shape = Replicate), size = 1.4, alpha = 1) +
  scale_color_manual(values = color)+
  scale_fill_manual(values = color)+
  scale_x_discrete(limits=c("mock", "flg22", "elf18", "R179_1.5", "R179_1.5_HK")) + 
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  annotate("text",
           x = l_EFR$Group,
           label = l_EFR$Letter,
           y=18,
           size= 3.5)+
  ggsave(paste("qPCR_R179_acute_EFR.pdf", sep=""), width=2.2, height=4.2)

################################################################################################################################
# Extended Data Figure 8d

# load data
qPCR_supernatant_FRK1 <- read.delim("Ext_Data_Fig_8d.txt", stringsAsFactors = F)

# add combined column
qPCR_supernatant_FRK1$treatment_genotype <- paste0(qPCR_supernatant_FRK1$treatment,"_", qPCR_supernatant_FRK1$genotype)

# define colors
colors <- c(  "mock" = "black", "wt" = "#ffc300", "dssAB" = "#386641")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_FRK1_s <- kruskal.test(treatment_genotype ~ fold_change, data=qPCR_supernatant_FRK1)
print(kw_FRK1_s)

dn_FRK1_s <- dunnTest (fold_change ~ treatment_genotype, data=qPCR_supernatant_FRK1, method=p.adjust.method)
dn_FRK1_s <- dn_FRK1_s$res
print(dn_FRK1_s)

l_FRK1_s<- cldList(comparison=dn_FRK1_s$Comparison, p.value=dn_FRK1_s$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_FRK1_s)

# plot supernatant FRK1
ggplot(qPCR_supernatant_FRK1, aes(x = treatment_genotype, y= fold_change)) +
  geom_bar(aes(fill = treatment), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, size=0.8) +
  geom_jitter(aes(color = treatment, shape = replicate), size = 1, alpha = 0.9) +
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  scale_x_discrete(limits=c("mock_Col0", "wt_Col0", "dssAB_Col0")) + 
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  annotate("text",
           x = l_FRK1_s$Group,
           label = l_FRK1_s$Letter,
           y=9.5,
           size= 4.5)+
  ggsave(paste("qPCR_supernatant_FRK1.pdf", sep=""), width=1.8, height=4)

################################################################################################################################
# Extended Data Figure 8e

# load data
qPCR_resuspension_FRK1 <- read.delim("Ext_Data_Fig_8e.txt", stringsAsFactors = F)

# add combined column
qPCR_resuspension_FRK1$treatment_genotype <- paste0(qPCR_resuspension_FRK1$treatment,"_", qPCR_resuspension_FRK1$genotype)

## statistics

kw_FRK1_r <- kruskal.test(treatment_genotype ~ fold_change, data=qPCR_resuspension_FRK1)
print(kw_FRK1_r)

dn_FRK1_r <- dunnTest (fold_change ~ treatment_genotype, data=qPCR_resuspension_FRK1, method=p.adjust.method)
dn_FRK1_r <- dn_FRK1_r$res
print(dn_FRK1_r)
write_tsv(dn_FRK1_r, "Ext_Data_Fig8e_stats.txt")

l_FRK1_r<- cldList(comparison=dn_FRK1_r$Comparison, p.value=dn_FRK1_r$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_FRK1_r)

# plot
ggplot(qPCR_resuspension_FRK1, aes(x = treatment_genotype, y= fold_change)) +
  geom_bar(aes(fill = treatment), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, size=0.8) +
  geom_jitter(aes(color = treatment, shape = replicate), size = 1, alpha = 0.9) +
  scale_color_manual(values = colors)+
  scale_fill_manual(values = colors)+
  scale_x_discrete(limits=c("mock_Col0", "wt_Col0", "dssAB_Col0",
                            "mock_efs", "wt_efs", "dssAB_efs")) + 
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  annotate("text",
           x = l_FRK1_r$Group,
           label = l_FRK1_r$Letter,
           y=8,
           size= 3.5)+
  ggsave(paste("qPCR_resuspension_FRK1.pdf", sep=""), width=2.6, height=4)

