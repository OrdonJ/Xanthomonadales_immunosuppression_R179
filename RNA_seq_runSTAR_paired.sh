#!/usr/bin/bash
#Author Tak Lee, tlee@mpipz.mpg.de
#bash script used to align the paired end RNAseq reads to the reference genome using STAR aligner
N=$(wc -l $1 | cut -d' ' -f1)
FILES=$(cat $1)
INPUTS=(${FILES//' '/ })
#echo $FILES
for ((i=0;i<$N;i++));
do
	ulimit -n 10000
        f1="${INPUTS[i]}_1.fq.gz"
	f2="${INPUTS[i]}_2.fq.gz"
	#echo $f1 $f2
	STAR \
	--outReadsUnmapped Fastx \
	--genomeDir /netscratch/dep_psl/grp_rgo/taklee/genomes/gindex_gtf_149 \
	--readFilesCommand zcat \
	--readFilesIn $f1 $f2 \
	--outSAMtype BAM SortedByCoordinate \
	--runThreadN 24 \
	--outFileNamePrefix bam/${INPUTS[i]}_ 

done

