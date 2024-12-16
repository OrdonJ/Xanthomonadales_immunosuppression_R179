# Figure 2

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)

rm(list = ls())

################################################################################################################################
# Figure 2a

setwd("C:/Users/Jana/Dropbox/Github/R179")

# read data files
flg22_root <- read.delim("Fig.2a.txt", stringsAsFactors = F)

# add combined columns
flg22_root$comb <- paste0(flg22_root$bacteria, "_", flg22_root$treatment)

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_root <- kruskal.test(comb ~ root_length, data=flg22_root)
print(kw_root)

dn_root <- dunnTest (root_length ~ comb, data=flg22_root, method=p.adjust.method)
dn_root <- dn_root$res
print(dn_root)

l_root <- cldList(comparison=dn_root$Comparison, p.value=dn_root$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_root)

# plot 
ggplot((data = flg22_root), aes(x=comb, y=root_length, fill=bacteria))+
  geom_boxplot(aes(fill= bacteria, alpha = 0.4), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape = replicate), size= 0.6, width = 0.35, alpha =0.4) +
  scale_x_discrete(limits=c("axenic_mock", "R179_mock", 
                            "axenic_flg22", "R179_flg22")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=4, angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_fill_manual(values = c("black", "#ffc300"))+
  annotate("text",
           x = l_root$Group,
           label = l_root$Letter,
           y=9.7,
           size= 4)+
  ggsave(paste("flg22_roots_pWERE.pdf", sep=""), width=3, height=4.3)

################################################################################################################################
# Figure 2b

# read data files
Atpep1_root <- read.delim("Fig.2b.txt", stringsAsFactors = F)

# add combined columns
Atpep1_root$comb <- paste0(Atpep1_root$bacteria, "_", Atpep1_root$treatment)

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_root <- kruskal.test(comb ~ root_length, data=Atpep1_root)
print(kw_root)

dn_root <- dunnTest (root_length ~ comb, data=Atpep1_root, method=p.adjust.method)
dn_root <- dn_root$res
print(dn_root)

l_root <- cldList(comparison=dn_root$Comparison, p.value=dn_root$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_root)

# plot
ggplot((data = Atpep1_root), aes(x=comb, y=root_length, fill=bacteria))+
  geom_boxplot(aes(fill= bacteria, alpha = 0.4), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape = replicate), size= 0.6, width = 0.35, alpha =0.4) +
  scale_x_discrete(limits=c("axenic_none", "R179_none", 
                            "axenic_Atpep1", "R179_Atpep1")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=4, angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_fill_manual(values = c("black", "#ffc300"))+
  annotate("text",
           x = l_root$Group,
           label = l_root$Letter,
           y=6,
           size= 4)+
  ggsave(paste("_Atpep1_roots_Col0.pdf", sep=""), width=2.7, height=4.3)

################################################################################################################################
# Figure 2c

# read data files
elf18_shoot <- read.delim("Fig.2c.txt", stringsAsFactors = F)

# add combined columns
elf18_shoot$comb <- paste0(elf18_shoot$bacteria, "_", elf18_shoot$treatment)

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

kw_shoot <- kruskal.test(comb ~ fresh_weight, data=elf18_shoot)
print(kw_shoot)

dn_shoot <- dunnTest (fresh_weight ~ comb, data=elf18_shoot, method=p.adjust.method)
dn_shoot <- dn_shoot$res
print(dn_shoot)

l_shoot <- cldList(comparison=dn_shoot$Comparison, p.value=dn_shoot$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_shoot)

# plot
ggplot((data = elf18_shoot), aes(x=comb, y=fresh_weight, fill=bacteria))+
  geom_boxplot(aes(fill= bacteria, alpha = 0.4), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape = replicate), size= 0.6, width = 0.35, alpha =0.4) +
  scale_x_discrete(limits=c("axenic_mock", "R179_mock", 
                            "axenic_elf18_std", "R179_elf18_std")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=4, angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_fill_manual(values = c("black", "#ffc300"))+
  annotate("text",
           x = l_shoot$Group,
           label = l_shoot$Letter,
           y=0.0053,
           size= 4)+
  ggsave(paste("elf18_shoots_Col0.pdf", sep=""), width=3, height=4.3)

################################################################################################################################
# Figure 2d

# read data files
RL<- read.delim("Fig.2d.txt", stringsAsFactors = F)

# add combined column
RL$SynCom_strain_treatment <- paste0(RL$SynCom,"_", RL$bacteria,"_", RL$treatment)


## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

kw_R179 <- kruskal.test(SynCom_strain_treatment ~ root_length, data=RL)
print(kw_R179)

dn_R179 <- dunnTest (root_length ~ SynCom_strain_treatment, data=RL, method=p.adjust.method)
dn_R179 <- dn_R179$res
print(dn_R179)

l_R179 <- cldList(comparison=dn_R179$Comparison, p.value=dn_R179$P.adj, threshold=0.001, remove.zero = FALSE)
print(l_R179)

# plot 
ggplot((data = RL), aes(x=SynCom_strain_treatment, y=root_length))+
  geom_boxplot(aes(alpha = 0.4), outlier.shape = NA,width = 0.7) +
  geom_jitter(aes(shape=replicate),alpha = 0.4, size = 0.6, width= 0.35) +
  scale_x_discrete(limits=c("none_none_mock", "NS_none_mock","NS_R179_mock",
                            "none_none_flg22", "NS_none_flg22", "NS_R179_flg22")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=12, angle = 90))+
  theme(axis.text.y = element_text(color="Black")) +
  annotate("text",
           x = l_R179$Group,
           label = l_R179$Letter,
           y=9.5,
           size= 3.0)+
  ggsave(paste("_Root_length_R179_NS.pdf", sep=""), width= 3.2, height=4.3)
