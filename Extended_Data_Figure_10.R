# Extended Data Figure 10

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)
library(rstatix)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Figure 10a

# read data
cfu_efs <- read.delim("Ext_Data_Fig_10a.txt", stringsAsFactors = F)

# combine columns
cfu_efs$g_b <- paste0(cfu_efs$genotype,"_", cfu_efs$bacteria)

# define colors
colors <- c("R179_wt" = "#ffc300", "R179_dssAB" = "#386641")

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test
kw_cfu_shoot <- kruskal.test(g_b ~ LOG10_cfu_g, data=cfu_efs)
print(kw_cfu_shoot)

dn_cfu_shoot <- dunnTest (LOG10_cfu_g ~ g_b, data=cfu_efs, method=p.adjust.method)
dn_cfu_shoot <- dn_cfu_shoot$res
print(dn_cfu_shoot)

l_cfu_shoot <- cldList(comparison=dn_cfu_shoot$Comparison, p.value=dn_cfu_shoot$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_cfu_shoot)

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
  annotate("text",
           x = l_cfu_shoot$Group,
           label = l_cfu_shoot$Letter,
           y=8,
           size= 4)+
  ggsave(paste("R179wt_dssAB_cfu_PRR_mut_shoot.pdf", sep=""), width=4.7, height=4.3)

################################################################################################################################
# Extended Data Figure 10b

# read data files
cfu_Xanthos <- read.delim("Ext_Data_Fig_10b_.txt", stringsAsFactors = F)

# add combined column
cfu_Xanthos$bacterium_treatment <- paste0(cfu_Xanthos$bacterium,"_", cfu_Xanthos$treatment)

# define colors
colors <- c("mock" = "black", "flg22" = "#DBC7AF")

cfu_Xanthos$treatment <- factor(cfu_Xanthos$treatment,
                                levels = c('mock','flg22'),ordered = TRUE)

# plot 
ggplot((data = cfu_Xanthos), aes(x=bacterium, y=log10, c))+
  geom_boxplot(aes(color = treatment), position = position_dodge(0.6),
               outlier.shape = NA, width = 0.5) +
  geom_jitter(aes(color = treatment, shape=replicate), position = position_dodge(0.6),
              alpha = 0.6, size = 0.8, width= ) +
  coord_cartesian(ylim = c(2, 10)) +
  theme_classic()+
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(color="Black", size=12, angle = 90))+
  theme(axis.text.y = element_text(color="Black")) +
  ggsave(paste("Colonization_Xanthos.pdf", sep=""), width=8.6, height=3.0)

## statistics

cfu_Xanthos %>%
  group_by(bacterium) %>%
  wilcox_test(data =., log10 ~ treatment) %>%
  add_significance("p") %>% 
  filter()->stats

################################################################################################################################
# Extended Data Figure 10c

# read data files
RL_Xanthos <- read.delim("Ext_Data_Fig_10c_root_length.txt", stringsAsFactors = F)
suppressive_score <- read.delim("Ext_Data_Fig_10c_suppress_score.txt", stringsAsFactors = F)

#add combined column
cfu_Xanthos$bacterium_treatment <- paste0(cfu_Xanthos$bacterium,"_", cfu_Xanthos$treatment)
RL_Xanthos$strain_treatment <- paste0(RL_Xanthos$strain,"_", RL_Xanthos$treatment)

#calculate means per bacterium
Colonization_mean <- aggregate(cfu_Xanthos[, 4], list(cfu_Xanthos$bacterium), mean)
Colonization_mean <- rename(Colonization_mean, bacterium = Group.1)
Colonization_mean <- rename(Colonization_mean, log10=x) 

RL_mean <- aggregate(RL_Xanthos[, 4], list(RL_Xanthos$strain_treatment), mean)
RL_mean <- rename(RL_mean, strain_treatment=Group.1)
RL_mean <- rename(RL_mean, root_length=x) 

joined <- left_join(Colonization_mean, suppressive_score)

# define color
color <- c("At-R-SPHERE" = "#BD312B", "At-Root-US" = "#E58E8B", 
           "Lj-SPHERE" = "#FFD800", "Chlamy" = "#008541", "pathogen" = "#6F0B56")

# plot
ggplot((data = joined), aes(x=log10, y=X1.max))+
  geom_point(aes(color = origin, alpha = 0.8))+
  geom_smooth(method='lm', formula= y~x, se=FALSE, col="black") +
  stat_regline_equation(label.y = 0.6, label.x = 7.5) +
  stat_cor(label.y = 0.65, label.x = 7.5) +
  scale_color_manual(values = color)+
  ylim(0.5, 1) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position = "none") +
  ggsave(paste("Correlation_cfu_suppressiveness.pdf", sep=""), width=2.6, height=4.8)
