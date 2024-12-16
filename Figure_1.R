# Figure 1

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Figure 1a

# read data files
RL_Xanthos <- read.delim("Fig.1a.txt", stringsAsFactors = F)

# add combined column
RL_Xanthos$strain_treatment <- paste0(RL_Xanthos$strain,"_", RL_Xanthos$treatment)

# calculate means per strain
RL_mean_strain <- aggregate(RL_Xanthos[, 4], list(RL_Xanthos$strain_treatment), mean)
RL_mean_strain <- rename(RL_mean_strain, strain_treatment=Group.1)
RL_mean_strain <- rename(RL_mean_strain, mean_root_length=x) 

RL_mean_strain <- separate(RL_mean_strain, strain_treatment, into=c("strain", "treatment"))

# make plot
ggplot((data = RL_mean_strain), aes(x=treatment, y=strain, size=mean_root_length))+
  geom_point(aes(alpha = 0.8),shape = 21)+
  scale_size_area(max_size = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.2))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position="left") +
  scale_y_discrete(limits=c("Xhp", "Xcv", "Xcc", "Chlamy123", "Chlamy24", 
                            "Lj197", "Lj193", "Lj173", "Lj165", "Lj156", "Lj143",
                            "Lj125", "Lj77", "Lj73", "Lj68", "Lj60", "Lj21", "Lj18",
                            "S773", "S772", "L148", "L131", "L70", "MF178", "MF92", 
                            "R983", "R916", "R690", "R667", "R604", "R559", "R179", 
                            "R96", "R65", "R71", "none"))+
  scale_x_discrete(limits=c("mock", "flg22", "Atpep1"), position = "top")+
  ggsave(paste("_Root_length_Xantho_ballonplot.pdf", sep=""), width=3.8, height=12)

### statistics

# filter per treatment
idx <- RL_Xanthos$treatment %in% c("mock") &  
  T
  RL_Xanthos_mock <- RL_Xanthos[idx,]
  
idx <- RL_Xanthos$treatment %in% c("flg22") &  
  T
  RL_Xanthos_flg22 <- RL_Xanthos[idx,]

idx <- RL_Xanthos$treatment %in% c("Atpep1") &  
    T
  RL_Xanthos_Atpep1 <- RL_Xanthos[idx,]
  
alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

# mock
kw_mock <- kruskal.test(strain_treatment ~ root_length, data=RL_Xanthos_mock)
print(kw_mock)

dunn.test(x=RL_Xanthos_mock$root_length, g=RL_Xanthos_mock$strain_treatment, method="bh")->full
cbind(as.data.frame(full$comparisons), as.data.frame(full$P.adjusted), as.data.frame(full$P), as.data.frame(full$Z)) -> full
separate(full, `full$comparisons`, c("well1", "well2"), " - ")->full
cbind(filter(full, well1=="none_mock" | well2=="none_mock"))->full
filter(full, `full$P` < 0.001)->full
as.data.frame(c(full$well1, full$well2))->full
colnames(full)<-c("dif.mock")
filter(full, dif.mock!="none_mock")->full
write_tsv(full, "_Xantho_dif_mock_0.001.txt")

# flg22
kw_flg22 <- kruskal.test(strain_treatment ~ root_length, data=RL_Xanthos_flg22)
print(kw_flg22)

dunn.test(x=RL_Xanthos_flg22$root_length, g=RL_Xanthos_flg22$strain_treatment, method="bh")->full
cbind(as.data.frame(full$comparisons), as.data.frame(full$P.adjusted), as.data.frame(full$P), as.data.frame(full$Z)) -> full
separate(full, `full$comparisons`, c("well1", "well2"), " - ")->full
cbind(filter(full, well1=="none_flg22" | well2=="none_flg22"))->full
filter(full, `full$P` < 0.001)->full
as.data.frame(c(full$well1, full$well2))->full
colnames(full)<-c("dif.mock")
filter(full, dif.mock!="none_flg22")->full
write_tsv(full, "_Xantho_dif_flg22_0.001.txt")

# Atpep1
kw_Atpep1 <- kruskal.test(strain_treatment ~ root_length, data=RL_Xanthos_Atpep1)
print(kw_Atpep1)

dunn.test(x=RL_Xanthos_Atpep1$root_length, g=RL_Xanthos_Atpep1$strain_treatment, method="bh")->full
cbind(as.data.frame(full$comparisons), as.data.frame(full$P.adjusted), as.data.frame(full$P), as.data.frame(full$Z)) -> full
separate(full, `full$comparisons`, c("well1", "well2"), " - ")->full
cbind(filter(full, well1=="none_Atpep1" | well2=="none_Atpep1"))->full
filter(full, `full$P` < 0.001)->full
as.data.frame(c(full$well1, full$well2))->full
colnames(full)<-c("dif.mock")
filter(full, dif.mock!="none_Atpep1")->full
write_tsv(full, "_Xantho_dif_Atpep1_0.001.txt")


################################################################################################################################
# Figure 1c

# read data files
cfu_soil <- read.delim("Fig.1c.txt", stringsAsFactors = F)

### statistics

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

kw_cfu_soil <- kruskal.test(strain ~ log10_cfu_g, data=cfu_soil)
print(kw_cfu_soil)

dn_cfu_soil <- dunnTest (log10_cfu_g ~ strain, data=cfu_soil, method=p.adjust.method)
dn_cfu_soil <- dn_cfu_soil$res
print(dn_cfu_soil)

l_cfu_soil <- cldList(comparison=dn_cfu_soil$Comparison, p.value=dn_cfu_soil$P.adj, threshold=0.05, remove.zero = FALSE)
print(l_cfu_soil)

# set color code
colors <- c("R179" = "#AB1A0D", "S772" = "#776403", "S773" = "#DBB706")

# plot 
ggplot((data = cfu_soil), aes(x=strain, y=log10_cfu_g))+
  geom_boxplot(alpha =1, size = 0.4, outlier.shape = NA) +
  geom_jitter(aes(color= strain, alpha = 0.8)) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_color_manual(values = colors)+
  annotate("text",
           x = l_cfu_soil$Group,
           label = l_cfu_soil$Letter,
           y=6,
           size = 4.0)+
  ggsave(paste("__Colonization_S772_S773.pdf", sep=""), width=3, height=4)


################################################################################################################################
# Figure 1d

# read data files
RL_Ath_Lj <- read.delim("Fig.1d.txt", stringsAsFactors = F)

# add combined column
RL_Ath_Lj$strain_treatment <- paste0(RL_Ath_Lj$strain,"_", RL_Ath_Lj$treatment)

# calculate means per strain
RL_mean_strain <- aggregate(RL_Ath_Lj[, 4], list(RL_Ath_Lj$strain_treatment), mean)
RL_mean_strain <- rename(RL_mean_strain, strain_treatment=Group.1)
RL_mean_strain <- rename(RL_mean_strain, mean_root_length=x) 

RL_mean_strain <- separate(RL_mean_strain, strain_treatment, into=c("strain", "treatment"))

# make plot
ggplot((data = RL_mean_strain), aes(x=treatment, y=strain, size=mean_root_length))+
  geom_point(aes(alpha = 0.8), shape = 21)+
  scale_size_area(max_size = 8) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  theme(legend.position="left") +
  #scale_fill_manual(values = color)+
  scale_y_discrete(limits=c("none", "R179", "Lj18", "Lj60", "Chlamy24", "Chlamy123", "Xcc", "Xcv"))+
  scale_x_discrete(limits=c("mock", "flg22"))+
  ggsave(paste("_Root_length_on_Lj_ballonplot.pdf", sep=""), width=3.4, height=3.9)

### statistics

# filter per treatment
idx <- RL_Ath_Lj$treatment %in% c("mock") &  
  T
RL_Ath_Lj_mock <- RL_Ath_Lj[idx,]

idx <- RL_Ath_Lj$treatment %in% c("flg22") &  
  T
RL_Ath_Lj_flg22 <- RL_Ath_Lj[idx,]

alpha <- 0.05 
p.adjust.method <- "bh"

# Kruskal-Wallis test

# mock
kw_mock <- kruskal.test(strain_treatment ~ root_length, data=RL_Ath_Lj_mock)
print(kw_mock)

dunn.test(x=RL_Ath_Lj_mock$root_length, g=RL_Ath_Lj_mock$strain_treatment, method="bh")->full
cbind(as.data.frame(full$comparisons), as.data.frame(full$P.adjusted), as.data.frame(full$P), as.data.frame(full$Z)) -> full
separate(full, `full$comparisons`, c("well1", "well2"), " - ")->full
cbind(filter(full, well1=="none_mock" | well2=="none_mock"))->full
filter(full, `full$P` < 0.05)->full
as.data.frame(c(full$well1, full$well2))->full
colnames(full)<-c("dif.mock")
filter(full, dif.mock!="none_mock")->full
write_tsv(full, "Xantho_Lj_dif_mock_0.05.txt")

# flg22
kw_flg22 <- kruskal.test(strain_treatment ~ root_length, data=RL_Ath_Lj_flg22)
print(kw_flg22)

dunn.test(x=RL_Ath_Lj_flg22$root_length, g=RL_Ath_Lj_flg22$strain_treatment, method="bh")->full
cbind(as.data.frame(full$comparisons), as.data.frame(full$P.adjusted), as.data.frame(full$P), as.data.frame(full$Z)) -> full
separate(full, `full$comparisons`, c("well1", "well2"), " - ")->full
cbind(filter(full, well1=="none_flg22" | well2=="none_flg22"))->full
filter(full, `full$P` < 0.05)->full
as.data.frame(c(full$well1, full$well2))->full
colnames(full)<-c("dif.mock")
filter(full, dif.mock!="none_flg22")->full
write_tsv(full, "Xantho_Lj_dif_flg22_0.05.txt")
