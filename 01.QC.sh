#!/bin/bash
script=/Share/BP/wangrui/bin/QC/QC_plus_rm_CAT_primer.pl
indir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/00.raw_data
outdir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/01.clean_data
sample=$1
mkdir -p $outdir
cat $sample|while read sp
do 
	echo "perl $script -indir $indir -outdir $outdir -sample $sp -end 2 -scRNA 1" >$sp.qc.tmp.sh
	qsub -cwd -l vf=3g,io=0,p=1 -V $sp.qc.tmp.sh
done

