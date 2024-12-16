# Extended Data Figure 6

library("scales")
library("grid")
library("vegan")
library("tidyverse")
library("rstatix")

rm(list = ls())

# load plotting functions
source("C:/Users/Jana/Dropbox/Github/R179/plotting_functions.R")
source("C:/Users/Jana/Dropbox/Github/R179/cpcoa.func.R")

setwd("C:/Users/Jana/Dropbox/Github/R179")

################################################################################################################################
# Extended Data Figure 6a

# files
design.file <- paste("Ext_Data_Fig_6abde_design.txt", sep="")
otu_table.file <- paste("Ext_Data_Fig6abde_ASV.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# spike normalized tables
design$depth <- colSums(otu_table)
idx <- rownames(otu_table) == "pBCC78"
pBCC78 <- otu_table[idx, ]
pBCC78 <- unlist(pBCC78)
design <- cbind(design, pBCC78)
otu_table <- otu_table[!idx,  ]
otu_table_norm <- sweep(otu_table, 2, pBCC78, `/`)

design$cond <- paste(design$SynCom, design$R179_deriv)

######## CPCoA Bray-Curtis
colors <- data.frame(group=c(
  "NS4 wt",
  "NS4 dssAB",
  "S2 wt",
  "S2 dssAB"
),
color=c(
  "#ffc300",
  "#386641",
  "#ffc300",
  "#386641"
))

shapes <- data.frame(group=c(
  "NS4 wt",
  "NS4 dssAB",
  "S2 wt",
  "S2 dssAB"
),
shape=c(
  15,
  15,
  17,
  17
))

sqrt_transform <- T

idx <- design$compartment %in% c("root") &  
  design$SynCom %in% c("S2", "NS4") &
  design$R179_deriv %in% c("wt", "dssAB") &
  T

d <- design[idx, ]

idx <- colnames(otu_table_norm) %in% d$SampleID
bray_curtis_d <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale(bray_curtis_d ~ SynCom * R179_deriv + Condition(bio_rep * tech_rep), data=d, add=F, sqrt.dist=sqrt_transform)

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
  ggsave(paste("CPCoA_roots.pdf", sep=""),  height=4, width=3)

################################################################################################################################
# Extended Data Figure 6b

# remove all strain except for Root179
idx <- rownames(otu_table_norm) %in% c("Root179") 
otu_table_R179 <- otu_table_norm[idx,  ]

# dataframe with additional information
df_R179 <- melt(as.matrix(otu_table_R179))
colnames(df_R179) <- c("strain", "SampleID", "RA")

df_R179$compartment <- design$compartment[match(df_R179$SampleID, design$SampleID)]
df_R179$deriv <- design$R179_deriv[match(df_R179$SampleID, design$SampleID)]
df_R179$SynCom <- design$SynCom[match(df_R179$SampleID, design$SampleID)]
df_R179$bio.rep <- design$bio_rep[match(df_R179$SampleID, design$SampleID)]
df_R179$tec.rep <- design$tech_rep[match(df_R179$SampleID, design$SampleID)]

# add combined column
df_R179$deriv_SynCom <- paste0(df_R179$deriv,"_", df_R179$SynCom)

# subset samples
idx <- df_R179$compartment %in% c("root") &
  T
df_root_R179 <- df_R179[idx, ]

df_root_R179$deriv<- factor(df_root_R179$deriv, levels= c("wt", "dssAB"))

#set color code
colors <- c("wt" = "#ffc300", "dssAB" = "#386641")

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

# plot root
ggplot(df_root_R179, aes(x=deriv_SynCom, y=RA)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  geom_jitter(aes(color = deriv, shape = bio.rep), size=2, alpha=0.9) +
  scale_x_discrete(limits=c("wt_NS4", "dssAB_NS4",
                            "wt_S2", "dssAB_S2")) +
  scale_color_manual(values = colors) +
  theme_bw() +
  main_theme +
  ggsave(paste("box_R179_roots_wt_dssAB.pdf", sep=""), width=1.8, height=5)

## statistics

df_root_R179%>%
  dunn_test(data =., RA ~ deriv_SynCom) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter()->root_stats

################################################################################################################################
# Extended Data Figure 6d

# CPcoA for XVMs samples wt_dssAB
idx <- design$compartment %in% c("XVM2") &
  T

d <- design[idx, ]

idx <- colnames(otu_table_norm) %in% d$SampleID
bray_curtis_d <- vegdist(t(otu_table_norm[, idx]), method="bray")

capscale.gen <- capscale(bray_curtis_d ~ SynCom * R179_deriv + Condition(bio_rep * tech_rep), data=d, add=F, sqrt.dist=sqrt_transform)

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
  ggtitle("CPCoA in medium", paste(format(100 * variance, digits=3), " % of variance; P=", format(p.val, digits=2), sep="")) +
  main_theme +
  theme(legend.position="bottom")+
  ggsave(paste("CPCoA_roots_only_exp_strains_medium.pdf", sep=""),  height=4, width=3)

################################################################################################################################
# Extended Data Figure 6e

# subset samples
idx <- df_R179$compartment %in% c("XVM2") &
  T
df_XVM2_R179 <- df_R179[idx, ]

df_XVM2_R179$deriv<- factor(df_XVM2_R179$deriv, levels= c("wt", "dssAB"))

#set color code
colors <- c("wt" = "#ffc300", "dssAB" = "#386641")

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
ggplot(df_XVM2_R179, aes(x=deriv_SynCom, y=RA)) +
  geom_boxplot(alpha=0.7, outlier.shape = NA) +
  geom_jitter(aes(color = deriv, shape = bio.rep), size=2, alpha=0.9) +
  scale_x_discrete(limits=c("wt_NS4", "dssAB_NS4",
                            "wt_S2", "dssAB_S2")) +
  scale_color_manual(values = colors) +
  theme_bw() +
  # theme_gray() +
  main_theme +
  #ylim (0, 0.6) +
  labs(y="Normalized to spike")+
  ggsave(paste("box_R179_XVM2_wt_dssAB.pdf", sep=""), width=2, height=5)

# statistics

df_XVM2_R179%>%
  dunn_test(data =., RA ~ deriv_SynCom) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  filter()->XVM2_stats

################################################################################################################################
# Extended Data Figure 6c

rm(list = ls())

# load plotting functions
source("C:/Users/Jana/Dropbox/Github/R179/plotting_functions.R")
source("C:/Users/Jana/Dropbox/Github/R179/cpcoa.func.R")

setwd("C:/Users/Jana/Dropbox/Github/R179")

# files
design.file <- paste("Ext_Data_Fig_6c_design.txt", sep="")
otu_table.file <- paste("Ext_Data_Fig_6c_ASV.txt", sep="")

# load data
design <- read.table(design.file, header=T, sep="\t")
otu_table <- read.table(otu_table.file, sep="\t", header=T, check.names=F)

# re-order data matrices
idx <- design$SampleID %in% colnames(otu_table)
design <- design[idx, ]

idx <- match(design$SampleID, colnames(otu_table))
otu_table <- otu_table[, idx]

# remove R179_wt, R179_dssA, R179_dssB, dssAB --> difference in S or NS community plus R179
idx <- rownames(otu_table) %in% c("R179_wt", "R179_dssA", "R179_dssB", "R179_dssAB")
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

# CPCoA Bray-Curtis
colors <- data.frame(group=c(
  "NS wt",
  "NS dssA",
  "NS dssB",
  "NS dssAB",
  "S wt",
  "S dssA",
  "S dssB",
  "S dssAB"
),
color=c(
  "#ffc300",
  "#A7C957",
  "#6A994E",
  "#386641",
  "#ffc300",
  "#A7C957",
  "#6A994E",
  "#386641"
))

shapes <- data.frame(group=c(
  "NS wt",
  "NS dssA",
  "NS dssB",
  "NS dssAB",
  "S wt",
  "S dssA",
  "S dssB",
  "S dssAB"
),
shape=c(
  15,
  15,
  15,
  15,
  17,
  17,
  17,
  17
))

sqrt_transform <- T

# CPcoA for root samples
idx <- design$compartment %in% c("root") &  
  design$SynCom %in% c("S", "NS") &
  design$R179_deriv %in% c("wt", "dssA", "dssB", "dssAB") &
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
  geom_point(alpha=.8, size=2.5) +
  stat_ellipse(aes(color=cond1), linetype=2, type="t", size=1.0, alpha=.5) +
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
  ggsave(paste("CPCoA_roots_wt_dssA_dssB_dssAB.pdf", sep=""),  height=7, width=5)
