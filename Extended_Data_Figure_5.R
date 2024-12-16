# Extended Data Figure 5

library(tidyverse)
library(dunn.test)
library(FSA)
library(rcompanion)
library(EnvStats)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Extended Data Fig. 5b

# read file
Kin_dss <- read.delim("Ext_Data_Fig.5b.txt", stringsAsFactors = F)

#set color code
colors <- c("R179wt" = "#ffc300", "R179dssA" = "#A7C957", "R179dssB" = "#6A994E", "R179dssAB" = "#386641")

# plot
ggplot(Kin_dss) + 
  geom_ribbon(aes(x = time_min, y = OD_mean,  fill=R179_deriv, ymin = OD_mean - OD_se/sqrt(3), ymax = OD_mean + OD_se/sqrt(3)), alpha=0.2) +
  geom_line(aes(x = time_min, y = OD_mean,  color=R179_deriv)) +
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text(size=12, color= "black"),
        axis.text.x = element_text( size=12, color= "black"),
        axis.title=element_text(size=12), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black"))+
  scale_colour_manual(values = colors) +
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(color="Black", hjust = 1))+
  ylab("OD600")+
  xlab("Time [min]")+
  ggsave(paste("Kin_R179wt_dssAB_XVM2.pdf", sep=""), width=5, height=5)

################################################################################################################################
# Extended Data Fig. 5c

# read files
cfu_dss <- read.delim("Ext_Data_Fig.5c.txt", stringsAsFactors = F)

## statistics

alpha <- 0.05 
p.adjust.method <- "bh"


# Kruskal-Wallis test

kw_cfu_dss <- kruskal.test(R179_deriv ~ cfu_log, data=cfu_dss)
print(kw_cfu_dss)

dn_cfu_dss <- dunnTest (cfu_log ~ R179_deriv, data=cfu_dss, method=p.adjust.method)
dn_cfu_dss <- dn_cfu_dss$res
print(dn_cfu_dss)

l_cfu_dss <- cldList(comparison=dn_cfu_dss$Comparison, p.value=dn_cfu_dss$P.adj, threshold=0.01, remove.zero = FALSE)
print(l_cfu_dss)

#set color code
colors <- c("R179wt" = "#ffc300", "R179dssA" = "#A7C957", "R179dssB" = "#6A994E", "R179dssAB" = "#386641")

# plot 
ggplot((data = cfu_dss), aes(x=R179_deriv, y=cfu_log))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = R179_deriv, shape = replicate, alpha = 0.8)) +
  coord_cartesian(ylim = c(1, 7)) +
  scale_x_discrete(limits=c("R179wt", "R179dssA", "R179dssB", "R179dssAB")) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", size=14))+
  theme(axis.text.y = element_text(color="Black")) +
  scale_color_manual(values = colors)+
  ggsave(paste("Colonization_Col0_R179_wt_dss.pdf", sep=""), width=3.0, height=2.65)
