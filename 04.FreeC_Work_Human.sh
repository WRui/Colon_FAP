#!/bin/bash
inpath=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/02.bam
sample=$1

outpath=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV

mkdir -p $outpath

for window in 10 #1 2 #5 10
do
    mkdir -p $outpath/${window}M
	cat $sample|while read sp
	do
		mkdir -p $outpath/${window}M/$sp
		echo "
[general]

chrLenFile = /Share/BP/wangrui/Database/hg19/hg19.genome
ploidy = 2
chrFiles = /Share/BP/wangrui/Database/hg19/hg19_chrSep
window = ${window}000000
maxThreads = 2
outputDir= $outpath/${window}M/$sp
BedGraphOutput = TRUE
samtools = /Share/BP/yanglu/software/samtools-1.2/samtools
[sample]

mateFile = $inpath/$sp/$sp.Mdu.sort.bam
inputFormat = BAM
mateOrientation = 0

			" > $sp.${window}M_config.tmp.txt

		echo "/Share/BP/wangrui/software/FREEC/freec -conf $sp.${window}M_config.tmp.txt " >$sp.${window}M.tmp.sh
		qsub -cwd -l vf=4g,io=0,p=3 -V $sp.${window}M.tmp.sh
	done
done


