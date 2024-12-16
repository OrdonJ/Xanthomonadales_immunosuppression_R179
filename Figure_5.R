# Figure 5

library("tidyverse")
library("dunn.test")
library("FSA")
library("EnvStats")
library("rcompanion")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Figure 5a

# load data
qPCR_R179 <- read.delim("Fig.5a.txt", stringsAsFactors = F)

# add combined column
qPCR_R179$treatment_genotype <- paste0(qPCR_R179$treatment,"_", qPCR_R179$genotype)

colors <- c(  "mock" = "black", "R179(OD1.5)" = "#ffc300", "R179_HK" = "#FF6200")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_FRK1 <- kruskal.test(treatment_genotype ~ fold_change, data=qPCR_R179)
print(kw_FRK1)

dn_FRK1 <- dunnTest (fold_change ~ treatment_genotype, data=qPCR_R179, method=p.adjust.method)
dn_FRK1 <- dn_FRK1$res
print(dn_FRK1)

l_FRK1<- cldList(comparison=dn_FRK1$Comparison, p.value=dn_FRK1$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_FRK1)

# plot
ggplot(qPCR_R179, aes(x = treatment_genotype, y= fold_change)) +
  geom_bar(aes(fill = treatment), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("mock_Col0", "R179(OD1.5)_Col0", "R179_HK_Col0",
                            "mock_fls2", "R179(OD1.5)_fls2", "R179_HK_fls2",
                            "mock_efr", "R179(OD1.5)_efr", "R179_HK_efr",
                            "mock_sobir", "R179(OD1.5)_sobir", "R179_HK_sobir",
                            "mock_exfxs", "R179(OD1.5)_exfxs", "R179_HK_exfxs")) +
  scale_fill_manual(values = colors)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=0.6) +
  geom_jitter(aes(color = treatment, shape = replicate), size = 0.8, alpha = 0.9) +
  scale_color_manual(values = colors)+
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=10, color= "black", angle = 90, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  coord_cartesian(ylim = c(0, 30)) +
  annotate("text",
           x = l_FRK1$Group,
           label = l_FRK1$Letter,
           y=28,
           size= 2.5,
           angle = 45)+
  ggsave(paste("qPCR_R179_acute_FRK1_PRR.pdf", sep=""), width=4, height=4)

################################################################################################################################
# Figure 5b

# load data
root_length <- read.delim("Fig.5b.txt", stringsAsFactors = F)

# add combined column
root_length$comb <- paste0(root_length$genotype,"_", root_length$treatment)

colors <- c(  "mock" = "black", "HK_R179" = "#FF6200")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_root <- kruskal.test(comb ~ rootlength, data=root_length)
print(kw_root)

dn_root <- dunnTest (rootlength ~ comb, data=root_length, method=p.adjust.method)
dn_root <- dn_root$res
print(dn_root)

l_root <- cldList(comparison=dn_root$Comparison, p.value=dn_root$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_root)

# plot root
ggplot((data = root_length), aes(x=comb, y=rootlength))+
  geom_boxplot(aes(fill= treatment, alpha = 0.5), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape=replicate), size= 0.8, width = 0.35, alpha =0.6)+
  scale_x_discrete(limits= c("Col0_mock", "Col0_HK_R179",
                             "efr_mock", "efr_HK_R179",
                             "sobir_mock", "sobir_HK_R179",
                             "sxfxe_mock", "sxfxe_HK_R179")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust= 0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_fill_manual(values = colors)+
  annotate("text",
           x = l_root$Group,
           label = l_root$Letter,
           y=8,
           size= 3.0)+
  ggsave(paste("Root_length_chronic_HK_R179_sub.pdf", sep=""), width=5, height=3.9)

################################################################################################################################
# Figure 5c

# read data files
ethylene_R179 <- read.delim("Fig.5c.txt", stringsAsFactors = F)

# add combined column
ethylene_R179$gXb <- paste0(ethylene_R179$genotype,"_", ethylene_R179$bacteria)

colors <- c(  "mock" = "black", "R179" = "#ffc300", "R179_HK" = "#FF6200")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_ethylene_R179 <- kruskal.test(gXb ~ ethylene.pmol.ml., data=ethylene_R179)
print(kw_ethylene_R179)

dn_ethylene_R179 <- dunnTest (ethylene.pmol.ml. ~ gXb , data=ethylene_R179, method=p.adjust.method)
dn_ethylene_R179 <- dn_ethylene_R179$res
print(dn_ethylene_R179)

l_ethylene_R179 <- cldList(comparison=dn_ethylene_R179$Comparison, p.value=dn_ethylene_R179$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_ethylene_R179)

# plot
ggplot(ethylene_R179, aes(x=gXb, y=ethylene.pmol.ml.)) +
  geom_bar(aes(fill = bacteria), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("Col0_mock", "Col0_R179", "Col0_R179_HK",
                            "fls2_mock", "fls2_R179", "fls2_R179_HK", 
                            "efr_mock", "efr_R179", "efr_R179_HK",
                            "efrxfls2xsobir1_mock", "efrxfls2xsobir1_R179", "efrxfls2xsobir1_R179_HK")) +
  scale_fill_manual(values = colors)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=0.6) +
  geom_jitter(aes(color = bacteria), size = 0.8, alpha = 0.9) +
  scale_color_manual(values = colors)+
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=12, color= "black"),
        axis.text.x = element_text( size=10, color= "black", angle = 90, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  annotate("text",
           x = l_ethylene_R179$Group,
           label = l_ethylene_R179$Letter,
           y=2.8,
           size= 2.5,
           angle = 45)+
  ggsave(paste("R179_ethylene_bar.pdf", sep=""), width=3.3, height=4.5)

################################################################################################################################
# Figure 5d

# read data files
cfu_efs <- read.delim("Fig.5d.txt", stringsAsFactors = F)

# add combined column
cfu_efs$g_b <- paste0(cfu_efs$genotype,"_", cfu_efs$bacteria)

#set color code
colors <- c("R179_wt" = "#ffc300", "R179_dssAB" = "#386641")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_cfu <- kruskal.test(g_b ~ LOG10_cfu_g, data=cfu_efs)
print(kw_cfu)

dn_cfu <- dunnTest (LOG10_cfu_g ~ g_b, data=cfu_efs, method=p.adjust.method)
dn_cfu <- dn_cfu$res
print(dn_cfu)

l_cfu <- cldList(comparison=dn_cfu$Comparison, p.value=dn_cfu$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_cfu)

# plot
ggplot((data = cfu_efs), aes(x=g_b, y=LOG10_cfu_g))+
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = bacteria, shape = replicate, alpha = 0.4), size = 1) +
  scale_x_discrete(limits=c("Col0_R179_wt", "Col0_R179_dssAB", "fls2_R179_wt","fls2_R179_dssAB",
                            "efr_R179_wt", "efr_R179_dssAB", "sobir_R179_wt", "sobir_R179_dssAB",
                            "fxexs_R179_wt", "fxexs_R179_dssAB")) +
  coord_cartesian(ylim = c(3, 8)) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=4))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_color_manual(values = colors)+
  ggsave(paste("R179wt_dssAB_cfu_PRR_mut_root.pdf", sep=""), width=5.6, height=3.0)

################################################################################################################################
# Figure 5e

# read data files
cfu_MAMPs <- read.delim("Fig.5e.txt", stringsAsFactors = F)

# add combined column
cfu_MAMPs$treatment_strain <- paste0(cfu_MAMPs$strain,"_", cfu_MAMPs$treatment)

## statistics 

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_cfu <- kruskal.test(treatment_strain ~ cfu_log10, data=cfu_MAMPs)
print(kw_cfu)

dn_cfu <- dunnTest (cfu_log10 ~ treatment_strain, data=cfu_MAMPs, method=p.adjust.method)
dn_cfu <- dn_cfu$res
print(dn_cfu)

l_cfu <- cldList(comparison=dn_cfu$Comparison, p.value=dn_cfu$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_cfu)

#set color code
colors <- c("R179_wt" = "#ffc300", "R179_dssAB" = "#386641")

# plot 
ggplot((data = cfu_MAMPs), aes(x=treatment_strain, y=cfu_log10))+
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(aes(color = strain, shape = replicate, alpha = 0.4), size = 1) +
  coord_cartesian(ylim = c(3, 8)) +
  scale_x_discrete(limits=c("R179_wt_mock", "R179_dssAB_mock","R179_wt_flg22",
                            "R179_dssAB_flg22", "R179_wt_elf18", "R179_dssAB_elf18")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=4))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_color_manual(values = colors)+
  annotate("text",
           x = l_cfu$Group,
           label = l_cfu$Letter,
           y=8,
           size= 4)+
  ggsave(paste("Colonization_R179_wt_dssAB_flg22_elf18_ppt.pdf", sep=""), width=4.2, height=3.0)
