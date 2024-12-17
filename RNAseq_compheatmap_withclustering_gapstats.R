#!/usr/bin/Rscript
#Author: Tak Lee, tlee@mpipz.mpg.de
#script used to find the optimal number of clusters using NbClust. Clustering and heatmap visualization by ComplexHeatmap
library(RColorBrewer)
library(NbClust)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(factoextra)
library(circlize)
library(cluster)

dat <- read.table("/netscratch/dep_psl/grp_rgo/taklee/Jana_RNAseq/defense_rootgrowth/clustering/DESeq2_vst_sva6_normalized_counts.txt",sep='\t',header=T,row.names=1,check.names=F)
degs <- scan("/netscratch/dep_psl/grp_rgo/taklee/Jana_RNAseq/defense_rootgrowth/diffs/DEGs/chronic_DEGs.txt",what=character())
sdat <- dat %>% select(
		       contains("axenic-none"),
		       contains("HKR179wthiOD-none"),
		       contains("R179wt-none"),
		       contains("NScomR179wt-none"),
		       contains("NScomHKR179wtlowOD-none"),
		       contains("ScomR179wt-none"),
		       contains("ScomHKR179wtlowOD-none")
		       )

scaleddat <- t(scale(t(sdat)))
filtdat <- scaleddat[degs,]

#determine the number of k with gap statistics:
#gap_stat <- clusGap(filtdat, FUN = kmeans, nstart = 4,K.max = 40, B = 500)
#k <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])

k <- 14

#plot heatmap
lgdcols <- colorRamp2(c(min(filtdat),-0.7,0,0.7,max(filtdat)),c("blue","#0080ff","grey45", "yellow","red"))
pdf("chronic_gapstat_clustering_normalizedcounts.pdf",height=14,width=10)
ht <- Heatmap(
	filtdat,
	col=lgdcols,
	clustering_distance_rows = "euclidean",
	clustering_method_rows="ward.D",
	row_km = k,
	clustering_distance_columns = "euclidean",
	clustering_method_columns = "ward.D",
	show_row_names = FALSE,
	column_names_gp = gpar(fontsize = 10),
	heatmap_legend_param = list(title="normalized counts")
	)

hm <- draw(ht)
dev.off()
klusters <- row_order(hm)

dir.create("clusters")
for (c in seq(1,k)) {
    genes <- rownames(filtdat[klusters[[c]],])
    fn <- paste("clusters/chronic_cluster",c,"_genes.txt",sep="")
    write.table(genes,file=fn,col.names=F,row.names=F,quote=F)
}



