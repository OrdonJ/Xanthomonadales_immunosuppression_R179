# Extended Data Figure 3

library(tidyverse)

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################

# read files
proteomics <- read.delim("Ext_Data_Fig.3.txt", stringsAsFactors = F)

# set color code
colors <- c("NO" = "black", "YES" = "red")

# subset proteomics
idx <- proteomics$replicate %in% c("I") &
  TRUE
repI <- proteomics[idx, ]
repII <- proteomics[!idx, ]

# plot
ggplot(data=repI, aes(x=log2_R179vsmock, y=log.pvalue., color=significant)) + 
  geom_point(alpha=0.2)+
  scale_color_manual(values = colors)+
  xlim(-7, 7) +
  ylim(0, 5.5) + 
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
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
ggsave(paste("proteomics_R179_mock_Col0_proteomics_repI.pdf", sep=""), width=3.5, height=5)

ggplot(data=repII, aes(x=log2_R179vsmock, y=log.pvalue., color=significant)) + 
  geom_point(alpha=0.2)+
  scale_color_manual(values = colors)+
  xlim(-7, 7) +
  ylim(0, 5.5) + 
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
        #legend.position="none",
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
ggsave(paste("proteomics_R179_mock_Col0_proteomics_repII.pdf", sep=""), width=3.5, height=5)
