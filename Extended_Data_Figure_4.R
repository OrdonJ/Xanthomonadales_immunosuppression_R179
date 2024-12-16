# Extended Data Figure 4

library("tidyverse")
library("rcompanion")
library("ggpubr")
library("ggrepel")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################

# read data
corr <- read.delim("Ext_Data_Fig.4.txt", stringsAsFactors = F)

# plot
ggplot((data = corr), aes(x=OD, y=counts))+
  geom_line(stat="smooth", method='lm', formula= y~x, se=FALSE, alpha=0.4) +
  geom_point(alpha = 0.8)+
  stat_regline_equation(size = 4, label.y = 2.4e+06) +
  stat_cor(size = 4) +
  theme_classic()+
  theme(axis.text.x = element_text(color="Black", angle=45, hjust=0.95))+
  theme(axis.text.y = element_text(color="Black")) +
  ggsave(paste("correlation_OD_cfu.pdf", sep=""), width=3.2, height=3.2)
