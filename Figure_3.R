# Figure 3

library("tidyverse")
library("rstatix")
library("scales")
library("grid")
library("vegan")

setwd("C:/Users/Jana/Dropbox/Github/R179")

rm(list = ls())

################################################################################################################################
# Figure 3b

# load data
qPCR_fp <- read.delim("Fig.3b.txt", stringsAsFactors = F)

# set colors
color <- c("none" ="black","R179_wt" = "#ffc300", "R179dssA" = "#A7C957", 
           "R179dssB" = "#6A994E")

# plot
ggplot(qPCR_fp, aes(x = strain, y= fold_change)) +
  geom_bar(aes(fill = strain), position = "dodge", stat="summary", fun="mean", width= 0.8, alpha=0.5)+
  scale_fill_manual(values = color)+
  geom_errorbar(stat = 'summary', position = 'dodge', width=0.3, alpha=0.9, size=1.0) +
  geom_jitter(aes(color= strain, shape = Replicate), alpha = 0.9, size = 1) +
  scale_color_manual(values = color)+
  theme(panel.background = element_rect(F))+
  theme(axis.line = element_line(T))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line(colour = "black", size = 0.5),
        panel.border = element_rect(colour = 'black', size=1, fill = NA),
        axis.text=element_text( size=14, color= "black"),
        axis.text.x = element_text( size=14, color= "black", angle = 45, hjust=0.95),
        axis.title=element_text(size=14), 
        legend.background=element_blank(),
        legend.key=element_blank(),
        strip.background = element_rect(color="black", fill=NA, size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black")) +
  ggsave(paste("qPCR_fp.pdf", sep=""), width=3, height=4)

## statistics

qPCR_fp %>%
  dunn_test(data =., fold_change ~ strain) %>%
  filter()->stats

################################################################################################################################
# Figure 3c

# load plotting functions
source("C:/Users/Jana/Dropbox/Github/R179/plotting_functions.R")
source("C:/Users/Jana/Dropbox/Github/R179/cpcoa.func.R")

# files
design.file <- paste("Fig.3cd_design.txt", sep="")
otu_table.file <- paste("Fig.3cd_ASV.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove R179 reads
idx <- rownames(otu_table) %in% c("R179_wt", "R179_dssA", "R179_dssB", "R179_dssAB", "Root179")
depleted_strains <- otu_table[idx, ]
otu_table <- otu_table[!idx,  ]

# spike normalized otu tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC23"
pBCC23 <- otu_table[idx, ]
pBCC23 <- unlist(pBCC23)
design <- cbind(design, pBCC23)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, pBCC23, `/`)

design$cond <- paste(design$SynCom, design$R179_deriv)

######## CPCoA Bray-Curtis
colors <- data.frame(group=c(
  "NS wt",
  "NS dssAB",
  "S wt",
  "S dssAB"
),
color=c(
  "#ffc300",
  "#386641",
  "#ffc300",
  "#386641"
))

shapes <- data.frame(group=c(
  "NS wt",
  "NS dssAB",
  "S wt",
  "S dssAB"
),
shape=c(
  15,
  15,
  17,
  17
))

sqrt_transform <- T

# CPcoA for root samples wt_dssAB
idx <- design$compartment %in% c("root") &  
  design$SynCom %in% c("S", "NS") &
  design$R179_deriv %in% c("wt", "dssAB") &
  T

d <- design[idx, ]

idx <- colnames(otu_table_norm) %in% d$SampleID
bray_curtis_d <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale(bray_curtis_d ~ SynCom * R179_deriv + Condition(biological.replicate * technical.replicate), data=d, add=F, sqrt.dist=sqrt_transform)

perm_anova.gen <- anova.cca(capscale.gen)

# variability tables, confidence intervals
var_tbl.gen <- variability_table(capscale.gen)
eig <- capscale.gen$CCA$eig
variance <- capscale.gen$CCA$tot.chi / capscale.gen$tot.chi
variance <- var_tbl.gen["constrained", "proportion"]
p.val <- perm_anova.gen[1, 4]

points_cpcoa <- capscale.gen$CCA$wa[,1:2]
colnames(points_cpcoa) <- c("x", "y")
points_cpcoa <- cbind(points_cpcoa, d[match(rownames(points_cpcoa), d$SampleID), ])

# plot CPCo 1 and 2
colors_cpcoa <- colors[colors$group %in% points_cpcoa$cond, ]
points_cpcoa$cond1 <- factor(points_cpcoa$cond, levels=colors_cpcoa$group)

shapes_cpcoa <- shapes[shapes$group %in% points_cpcoa$cond, ]
points_cpcoa$cond2 <- factor(points_cpcoa$cond, levels=shapes_cpcoa$group)

# plot
ggplot(points_cpcoa, aes(x=x, y=y, color=cond1, shape=cond2)) +
  geom_point(alpha=.7, size=1.0) +
  stat_ellipse(aes(color=cond1), linetype=2, type="t", size=0.7, alpha=.5) +
  scale_colour_manual(values=as.character(colors_cpcoa$color)) +
  scale_shape_manual(values=shapes_cpcoa$shape) +
  theme(axis.text.x=element_text(size=12)) +
  theme(axis.text.y=element_text(size=12)) +
  theme(axis.title=element_text(size=13)) +
  labs(x=paste("CPCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("CPCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) +
  ggtitle("CPCoA of root samples", paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="bottom")+
  ggsave(paste("CPCoA_roots_16S_woR179_abs_wt_dssAB.pdf", sep=""),  height=4, width=3)

################################################################################################################################
# Figure 3d

# select strains of interest 
idx <- rownames(depleted_strains) %in% c("R179_wt", "R179_dssAB") 
otu_table_BCs <- depleted_strains[idx,  ]

# spike normalize
otu_table_BCs <- sweep(otu_table_BCs, 2, pBCC23, `/`)

# dataframe with information plus RA
df_BCs <- melt(as.matrix(otu_table_BCs))
colnames(df_BCs) <- c("strain", "SampleID", "RA")

#match sampleID with additional information of design file
df_BCs$compartment <- design$compartment[match(df_BCs$SampleID, design$SampleID)]
df_BCs$deriv <- design$R179_deriv[match(df_BCs$SampleID, design$SampleID)]
df_BCs$SynCom <- design$SynCom[match(df_BCs$SampleID, design$SampleID)]
df_BCs$bio.rep <- design$biological.replicate[match(df_BCs$SampleID, design$SampleID)]
df_BCs$tec.rep <- design$technical.replicate[match(df_BCs$SampleID, design$SampleID)]

# add combined column
df_BCs$strain_deriv_SynCom <- paste0(df_BCs$strain, "_", df_BCs$deriv,"_", df_BCs$SynCom)

# error correct BCs
idx <- df_BCs$strain %in% c("R179_wt")
df_BC_wt <- df_BCs[idx,  ]
df_BC_wt$RA_norm <- df_BC_wt$RA / 1.7

idx <- df_BCs$strain %in% c("R179_dssAB")
df_BC_dssAB <- df_BCs[idx,  ]
df_BC_dssAB$RA_norm <- df_BC_dssAB$RA / 2.1

df_BCs_cor <- rbind(df_BC_wt, df_BC_dssAB)

df_BCs_cor$deriv<- factor(df_BCs_cor$deriv, levels= c("wt", "dssAB"))

#set color code
colors <- c("R179_wt" = "#ffc300", "R179_dssAB" = "#386641")

#Set main theme
main_theme <- theme(axis.line.x=element_line(color="black"),
                    axis.line.y=element_line(color="black"),
                    axis.ticks=element_line(color="black"),
                    axis.text=element_text(colour="black", size=10),
                    axis.text.x=element_text(colour="black", size=6, hjust = 1,vjust = 1,angle=45),
                    #use none to remove legend, can choose top, right, left
                    legend.position="top",
                    legend.title=element_blank(),
                    legend.background=element_blank(),
                    legend.key=element_blank(),
                    text=element_text(family="sans", size=15))

# plot
ggplot(df_BCs_cor, aes(x=strain_deriv_SynCom, y=RA_norm)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  geom_jitter(aes(color = strain, shape = bio.rep), size=1, alpha=0.9) +
  scale_x_discrete(limits=c("R179_wt_wt_NS", "R179_dssAB_dssAB_NS",
                            "R179_wt_wt_S", "R179_dssAB_dssAB_S")) +
  scale_color_manual(values = colors) +
  theme_bw() +
  # theme_gray() +
  main_theme +
  #ylim (0, 0.6) +
  labs(y="Normalized to spike")+
  ggsave(paste("box_abs_BC_roots_wt_dssAB.pdf", sep=""), width=1.8, height=5)

df_BCs_cor %>%
  group_by(SynCom) %>%
  dunn_test(data =., RA_norm ~ strain_deriv_SynCom) %>%
  filter()->stats
