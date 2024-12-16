# Extended Data Figure 1

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Figure 1e

# read data files
RL_L131 <- read.delim("Ext_Data_Fig.1e.txt", stringsAsFactors = F)

# add combined column
RL_L131$SynCom_strain_treatment <- paste0(RL_L131$SynCom,"_", RL_L131$strain,"_", RL_L131$treatment)

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

kw_L131 <- kruskal.test(SynCom_strain_treatment ~ root_length, data=RL_L131)
print(kw_L131)

dn_L131 <- dunnTest (root_length ~ SynCom_strain_treatment, data=RL_L131, method=p.adjust.method)
dn_L131 <- dn_L131$res
print(dn_L131)

l_L131 <- cldList(comparison=dn_L131$Comparison, p.value=dn_L131$P.adj, threshold=0.001, remove.zero = FALSE)
print(l_L131)

#set color code
colors <- c("mock" = "#434E56","flg22" = "#DBC7AF")

# plot 
ggplot((data = RL_L131), aes(x=SynCom_strain_treatment, y=root_length))+
  geom_boxplot(aes(fill= treatment, alpha = 0.4), outlier.shape = NA,width = 0.7) +
  geom_jitter(aes(color= treatment, shape=replicate),alpha = 0.4, size = 0.6, width= 0.35) +
  scale_x_discrete(limits=c("none_none_mock", "none_none_flg22", "NS_none_mock",
                            "NS_none_flg22","NS_L131_mock", "NS_L131_flg22")) +
  theme_classic()+
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(color="Black", size=10, angle = 90))+
  theme(axis.text.y = element_text(color="Black")) +
  annotate("text",
           x = l_L131$Group,
           label = l_L131$Letter,
           y=9.5,
           size= 3.0)+
  ggsave(paste("Root_length_L131_NS.pdf", sep=""), width= 3, height=4.3)

