#!/bin/bash
Human_ref=/Share/BP/wangrui/Database/bwa_index/hg19/hg19.fa
bwa=/Share/BP/wangrui/software/bwa-0.7.5a/bwa
indir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/01.clean_data
outdir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/02.bam
sample=$1
H_fai=/Share/BP/wangrui/Database/bwa_index/hg19/hg19.fa.fai
SAMTOOLS=/Share/BP/yanglu/software/samtools-1.2/samtools

mkdir -p $outdir

cat $sample|while read sp
do
	mkdir -p $outdir/$sp
	echo "$bwa mem -M -t 2 $Human_ref  $indir/$sp/$sp.R1.clean.fq.gz $indir/$sp/$sp.R2.clean.fq.gz| $SAMTOOLS view -u -b -S -t $H_fai -| $SAMTOOLS sort -m 200000000 - $outdir/$sp/$sp.sort" >$sp.mapping.tmp.sh
	qsub -cwd -l vf=5G,io=0,p=3 -V $sp.mapping.tmp.sh
	#mkdir -p $outdir/$sp
	#mv $indir/$sp/*bam $outdir/$sp
done
