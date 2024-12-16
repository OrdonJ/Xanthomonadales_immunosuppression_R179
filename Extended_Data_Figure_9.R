# Extended Data Figure 9

library("tidyverse")
library("dunn.test")
library("FSA")
library("rcompanion")
library("EnvStats")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Figure 9a

# read data
MAMPs <- read.delim("Ext_Data_Fig_9a.txt", stringsAsFactors = F)

# add combined column
MAMPs$gXb <- paste0(MAMPs$genotype,"_", MAMPs$bacteria)

#define colors
colors <- c(  "mock" = "black", 
              "R179_1" = "#ffc300", 
              "R179_0.1" = "#ffc300",
              "R179_0.01" = "#ffc300",
              "HK_R179_1" = "#FF6200",
              "HK_R179_0.1" = "#FF6200",
              "HK_R179_0.01" = "#FF6200")

# plot
ggplot(MAMPs, aes(x=gXb, y=ethylene.pmol.ml.)) +
  geom_bar(aes(fill = bacteria), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("Col0_mock", "Col0_R179_0.01", "Col0_R179_0.1", "Col0_R179_1",
                            "Col0_HK_R179_0.01", "Col0_HK_R179_0.1", "Col0_HK_R179_1",
                            "fls2_mock", "fls2_R179_0.01", "fls2_R179_0.1", "fls2_R179_1",
                            "fls2_HK_R179_0.01", "fls2_HK_R179_0.1", "fls2_HK_R179_1", 
                            "efr_mock", "efr_R179_0.01", "efr_R179_0.1", "efr_R179_1", 
                            "efr_HK_R179_0.01", "efr_HK_R179_0.1", "efr_HK_R179_1",
                            "efrxfls2xsobir_mock", "efrxfls2xsobir_R179_0.01", "efrxfls2xsobir_R179_0.1", "efrxfls2xsobir_R179_1",
                            "efrxfls2xsobir_HK_R179_0.01", "efrxfls2xsobir_HK_R179_0.1", "efrxfls2xsobir_HK_R179_1")) +
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
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  ggsave(paste("R179_ethylene_bar_dosage.pdf", sep=""), width=5, height=4.5)

################################################################################################################################
# Extended Data Figure 9c

# read data
MAMPs_shoot <- read.delim("Ext_Data_Fig_9c.txt", stringsAsFactors = F)

# add combined columns
MAMPs_shoot$comb <- paste0(MAMPs_shoot$genotype,"_", MAMPs_shoot$bacteria, "_", MAMPs_shoot$treatment)

#define colors
colors <- c(  "axenic" = "black", 
              "R179" = "#ffc300") 
              
## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_shoot <- kruskal.test(comb ~ fresh_weight, data=MAMPs_shoot)
print(kw_shoot)

dn_shoot <- dunnTest (fresh_weight ~ comb, data=MAMPs_shoot, method=p.adjust.method)
dn_shoot <- dn_shoot$res
print(dn_shoot)

l_shoot <- cldList(comparison=dn_shoot$Comparison, p.value=dn_shoot$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_shoot)

# plot
ggplot((data = MAMPs_shoot), aes(x=comb, y=fresh_weight, fill=bacteria))+
  geom_boxplot(aes(fill= bacteria, alpha = 0.5), outlier.shape = NA,
               position = position_dodge(0.6)) +
  geom_jitter(aes(shape = replicate), size= 0.8, width = 0.35, alpha =0.6) +
  scale_x_discrete(limits=c("Col0_axenic_mock", "Col0_R179_mock", "Col0_axenic_elf18_R179", 
                            "Col0_R179_elf18_R179","efr_axenic_mock", "efr_R179_mock",
                            "efr_axenic_elf18_R179", "efr_R179_elf18_R179")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14, angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_fill_manual(values = c("black", "#ffc300"))+
  annotate("text",
           x = l_shoot$Group,
           label = l_shoot$Letter,
           y=0.008,
           size= 4)+
  ggsave(paste("elf18_shoots_Col0_efr.pdf", sep=""), width=3.5, height=4.3)

################################################################################################################################
# Extended Data Figure 9d

# read data files
MAMPs_Atpep1 <- read.delim("Ext_Data_Fig_9d.txt", stringsAsFactors = F)

# add combined column
MAMPs_Atpep1$gXb <- paste0(MAMPs_Atpep1$genotype,"_", MAMPs_Atpep1$bacteria)

# define colors
colors2 <- c("mock" = "black", "AtPep1" = "#DBC7AF")

# plot
ggplot(MAMPs_Atpep1, aes(x=gXb, y=ethylene.pmol.ml.)) +
  geom_bar(aes(fill = bacteria), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("Col0_mock", "Col0_AtPep1",
                            "fls2_mock", "fls2_AtPep1", 
                            "efr_mock", "efr_AtPep1",
                            "efrxfls2xsobir1_mock", "efrxfls2xsobir1_AtPep1")) +
  scale_fill_manual(values = colors2)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=0.6) +
  geom_jitter(aes(color = bacteria), size = 0.8, alpha = 0.9) +
  scale_color_manual(values = colors2)+
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
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  ggsave(paste("R179_ethylene_bar_Atpep1.pdf", sep=""), width=2.8, height=4.5)

################################################################################################################################
# Extended Data Figure 9e

# read files
eMax <- read.delim("Ext_Data_Fig_9e.txt", stringsAsFactors = F)

# add combined column
eMax$gXb <- paste0(eMax$genotype,"_", eMax$bacteria)

# define colors
colors <- c(  "mock" = "black",
              "AtPep1" = "#DBC7AF",
              "R179" = "#ffc300", 
              "Xap" = "#c70039",
              "Xcr" = "#900C3F")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_eMax <- kruskal.test(gXb ~ ethylene.pmol.ml., data=eMax)
print(kw_eMax)

dn_eMax <- dunnTest (ethylene.pmol.ml. ~ gXb , data=eMax, method=p.adjust.method)
dn_eMax <- dn_eMax$res
print(dn_eMax)

l_eMax <- cldList(comparison=dn_eMax$Comparison, p.value=dn_eMax$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_eMax)

# plot
ggplot(eMax, aes(x=gXb, y=ethylene.pmol.ml.)) +
  geom_bar(aes(fill = bacteria), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_x_discrete(limits=c("efrxfls2_mock", "efrxfls2_AtPep1",
                            "efrxfls2_Xap", "efrxfls2_Xcr", "efrxfls2_R179",
                            "efrxfls2xrlp1_mock", "efrxfls2xrlp1_AtPep1",
                            "efrxfls2xrlp1_Xap", "efrxfls2xrlp1_Xcr",
                            "efrxfls2xrlp1_R179")) +
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
        axis.text.x = element_text( size=12, color= "black", angle = 45, vjust=0.5),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  annotate("text",
           x = l_eMax$Group,
           label = l_eMax$Letter,
           y=9.5,
           size= 4)+
  ggsave(paste("R179_ethylene_bar.pdf", sep=""), width=3.3, height=4.5)
