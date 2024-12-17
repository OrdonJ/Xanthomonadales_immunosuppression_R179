#!/usr/bin/Rscript
#Author: Tak Lee, tlee@mpipz.mpg.de
#Rscript used to process raw counts file with DESeq2, find surrogate variables with sva, plot PCA plots and find differentially expressed genes with DESeq2
library(sva)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(limma)
library(foreach)



cnts <- read.table("counts.txt", sep='\t',header=T,row.names=1,check.names=F)
cnts$Length <- NULL
cnts$"R179dssAB-HKR179wthiOD_rep4" <- NULL
coldata <- read.table("coldata.txt", sep='\t',header=T,row.names=1,check.names=F)
coldata <- coldata[!(row.names(coldata) %in% c("R179dssAB-HKR179wthiOD_rep4")),]
ddsMat <- DESeqDataSetFromMatrix(countData = cnts,
                                 colData = coldata,
                                 design = ~ condition)

dds <- estimateSizeFactors(ddsMat)
dat  <- counts(dds, normalized = TRUE)

###############filtering the lowly expressed genes################

geneidx  <- rowMeans(dat) > 1
dat <- dat[geneidx,]

###############finding the surrogate variables################
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1,colData(dds))
#manual detection of surrogate variable
#svseq <- svaseq(dat, mod, mod0, n.sv = 6)

dmicnt <- rep(0,50)
if(FALSE){
#auto detection of surrogate variables
for (i in 1:2){
	svobj <- sva(dat, mod, mod0)
	dmicnt[svobj$n.sv] <- dmicnt[svobj$n.sv] +1
}
svnum <- which.max(dmicnt)
print(svnum)
}
svnum <- 6 
for (k in 6:svnum){
	print(k)
	svseq <- svaseq(dat, mod, mod0, n.sv = k)
	ddssva <- dds 
	fmla <- "~" 
	for (j in 1:k){
		newcol <- svseq$sv[,j]
		colData(ddssva) <- cbind(colData(ddssva), newcol)
		names(colData(ddssva))[length(names(colData(ddssva)))] <- paste0("SV",j)
		if (j == 1){fmla=paste0(fmla,paste0("SV",j))} else {
			fmla=paste(fmla,paste0("SV",j),sep="+")}
	}
	fmla=paste(fmla,"condition",sep="+")
	print(fmla)
	design(ddssva) <- as.formula(fmla)
	samples <- coldata$condition
	design <- model.matrix(~samples)

	################plotting PCA#################
	#vst normalization of data
	vsd <- vst(dds)[geneidx,]
	samples <- coldata$condition
	design <- model.matrix(~samples)
	rldData <- assay(vsd) %>%
		  removeBatchEffect(covariates = svseq$sv, design = design)
	sampleIdx <- 1:length(samples)
	#two different treatments, group1 and group2

	group1idx <- sampleIdx[colnames(rldData) %in% c("R179wt-HKR179wthiOD_rep1","R179wt-HKR179wthiOD_rep2","R179wt-HKR179wthiOD_rep3","R179wt-HKR179wthiOD_rep4","R179dssAB-HKR179wthiOD_rep1","R179dssAB-HKR179wthiOD_rep2","R179dssAB-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep1","axenic-HKR179wthiOD_rep2","axenic-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep4","axenic-MgSO4_rep1","axenic-MgSO4_rep2","axenic-MgSO4_rep3","axenic-MgSO4_rep4","axenic-R179wthiOD_rep1","axenic-R179wthiOD_rep2","axenic-R179wthiOD_rep3","axenic-R179wthiOD_rep4","axenic-R179dssABhiOD_rep1","axenic-R179dssABhiOD_rep2","axenic-R179dssABhiOD_rep3","axenic-R179dssABhiOD_rep4")]

	group2idx <- sampleIdx[!colnames(rldData) %in% c("R179wt-HKR179wthiOD_rep1","R179wt-HKR179wthiOD_rep2","R179wt-HKR179wthiOD_rep3","R179wt-HKR179wthiOD_rep4","R179dssAB-HKR179wthiOD_rep1","R179dssAB-HKR179wthiOD_rep2","R179dssAB-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep1","axenic-HKR179wthiOD_rep2","axenic-HKR179wthiOD_rep3","axenic-HKR179wthiOD_rep4","axenic-MgSO4_rep1","axenic-MgSO4_rep2","axenic-MgSO4_rep3","axenic-MgSO4_rep4","axenic-R179wthiOD_rep1","axenic-R179wthiOD_rep2","axenic-R179wthiOD_rep3","axenic-R179wthiOD_rep4","axenic-R179dssABhiOD_rep1","axenic-R179dssABhiOD_rep2","axenic-R179dssABhiOD_rep3","axenic-R179dssABhiOD_rep4")]

	group <- c("all","subset1","subset2")
	top <- 500
	allgn <- length(rownames(rldData))
	allidx <- list(sampleIdx,group1idx,group2idx)
	for(i in 1:length(allidx)){
		idx <- allidx[[i]]
		pca <- prcomp(t(rldData[, idx]))
		length(rownames(rldData))
		for (n in c(top,allgn)){
		#for top 500 genes that contribute the most to the variance of PC1 and PC2.
			p1 <- order(abs(pca$rotation[,1]),decreasing=TRUE)[1:n]
			p2 <- order(abs(pca$rotation[,2]),decreasing=TRUE)[1:n]
			g1 <- names(pca$rotation[,1][p1])
			g2 <- names(pca$rotation[,2][p2])
			p1p2 <- union(g1,g2)
			pca2 <- prcomp(t(rldData[p1p2, idx]))
			pc1 <- pca2$x[,1]
			pc2 <- pca2$x[,2]
			percentVar <- pca2$sdev^2/sum(pca2$sdev^2)
			percentVar <- round(100 * percentVar)

			pcaData <- data.frame(PC1 = pc1, PC2 = pc2, Group = colData(dds)[idx, 1], ID = rownames(colData(dds))[idx])
			ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
				geom_point(size = 3) +
				xlab(paste0("PC1: ",percentVar[1],"% variance")) +
				ylab(paste0("PC2: ",percentVar[2],"% variance")) +
				#scale_colour_manual(name = 'SynCom',values = cols[colorIdx]) +
				stat_ellipse(aes(x = PC1, y = PC2, group = Group), type = 't', linetype = 2, level = 0.8) +
				coord_fixed(1) +
				theme_classic() +
				geom_text(aes(label=ID),vjust=2) +
				theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
				      legend.text.align = 0,
				      axis.text = element_text(size = 13),
				      axis.title = element_text(size = 14),
				      legend.text=element_text(size= 13),
				      legend.title = element_text(size = 14))
			if(n == allgn){
				outpdf <- paste(group[i],"_PCA_sva",k,".pdf",sep="")
				ggsave(outpdf, width = 30,height=30)
			} else {
				outpdf <- paste(group[i],"_PCA_top",n,"_sva",k,".pdf",sep="")
				ggsave(outpdf, width = 30,height=30)
			}
			#outjpg <- paste(group[i],"_PCA_top500.jpg",sep="")
			#ggsave(outjpg, width = 30,height=30)
		}

	}
}
if(FALSE){
##################Differential Expression of genes###################
ddssva <- DESeq(ddssva[geneidx, ])
contrasts <- list(c("axenic-none","HKR179wtlowOD-none"),
		  c("axenic-none","R179wt-none"),
		  c("axenic-none","R179dssAB-none"),
		  c("R179wt-none","R179dssAB-none"),
		  c("HKR179wtlowOD-none","R179wt-none"),
		  c("HKR179wtlowOD-none","R179dssAB-none"),
		  c("axenic-none","NScomHKR179wtlowOD-none"),
		  c("axenic-none","ScomHKR179wtlowOD-none"),
		  c("NScomHKR179wtlowOD-none","ScomHKR179wtlowOD-none"),
		  c("NScomHKR179wtlowOD-none","NScomR179wt-none"),
		  c("NScomHKR179wtlowOD-none","NScomR179dssAB-none"),
		  c("NScomR179wt-none","NScomR179dssAB-none"),
		  c("ScomHKR179wtlowOD-none","ScomR179wt-none"),
		  c("ScomHKR179wtlowOD-none","ScomR179dssAB-none"),
		  c("ScomR179wt-none","ScomR179dssAB-none"),
		  c("ScomHKR179wtlowOD-none","R179wt-none"),
		  c("ScomHKR179wtlowOD-none","NScomR179wt-none"),
		  c("ScomHKR179wtlowOD-none","NScomR179dssAB-none"),
		  c("ScomR179wt-none","NScomR179wt-none"),
		  c("ScomR179dssAB-none","NScomR179dssAB-none"),
		  c("axenic-none","HKR179wthiOD-none"),
		  c("HKR179wthiOD-none","R179wt-none"),
		  c("HKR179wthiOD-none","HKR179wtlowOD-none"),
		  c("axenic-MgSO4","axenic-none"),
		  c("axenic-MgSO4","axenic-R179wthiOD"),
		  c("axenic-MgSO4","axenic-R179dssABhiOD"),
		  c("axenic-MgSO4","axenic-HKR179wthiOD"),
		  c("axenic-R179wthiOD","axenic-HKR179wthiOD"),
		  c("axenic-R179wthiOD","axenic-R179dssABhiOD"),
		  c("axenic-HKR179wthiOD","axenic-R179dssABhiOD"),
		  c("axenic-MgSO4","R179wt-HKR179wthiOD"),
		  c("axenic-MgSO4","R179dssAB-HKR179wthiOD"),
		  c("R179wt-HKR179wthiOD","R179dssAB-HKR179wthiOD"),
		  c("axenic-HKR179wthiOD","R179wt-HKR179wthiOD"),
		  c("axenic-HKR179wthiOD","R179dssAB-HKR179wthiOD"),
		  c("R179wt-none","R179wt-HKR179wthiOD"),
		  c("R179dssAB-none","R179wt-HKR179wthiOD")
		  )
for(i in 1:length(contrasts)){
	ctrl <- contrasts[[i]][1]
	trt <- contrasts[[i]][2]
	res <- results(ddssva, contrast=c("condition",trt,ctrl))
	outf <- paste("AthTAIR10_Jana_",ctrl,".VS.",trt,"_diff.txt",sep="")
	write.table(res,file=outf,quote=F, sep="\t")
	#print(outf)
	}
}
