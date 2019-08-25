#!/bin/bash
sample=$1
#indir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/02.bam/Old
indir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/02.bam
cat $sample|while read sp
do 
	echo "java -jar /Share/BP/wangrui/software/picard-tools-1.130/picard.jar MarkDuplicates I=$indir/$sp/$sp.sort.bam O=$indir/$sp/$sp.Mdu.sort.bam M=$indir/$sp/marked_dup_metrics.txt" >$sp.MarkDuplicate.tmp.sh
	qsub -cwd -l vf=3G,io=0,p=1 -V $sp.MarkDuplicate.tmp.sh
done
