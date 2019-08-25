#!/bin/bash
sample=$1
#indir=/WPSnew/wangrui/Project/Cancer_organoise/DNA/03.CNV
indir=/WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV
outfile=$2
for window in 5 10 #1 2 5 10
do

	cut -f 1,2 $indir/${window}M/V15-LN8/*_sample.cpn > $indir/${window}M/${window}M_Merge_read.txt
	cut -f 1,2 $indir/${window}M/V15-LN8/*_sample.cpn |grep -v "Chr" > $indir/${window}M/${window}M_Merge_GC.txt
 
	echo -e -n "Chr\tWindow" > $indir/${window}M/${window}M_read_header.txt
	cat $sample| while read sp
	do
		echo -e -n "\t$sp" >> $indir/${window}M/${window}M_read_header.txt
		paste $indir/${window}M/${window}M_Merge_read.txt <(cut -f 3 $indir/${window}M/$sp/*_sample.cpn) > $indir/${window}M/${window}M_Merge_read_tmp.txt
		mv $indir/${window}M/${window}M_Merge_read_tmp.txt $indir/${window}M/${window}M_Merge_read.txt

		paste $indir/${window}M/${window}M_Merge_GC.txt <(cut -f 3 $indir/${window}M/$sp/*bam_ratio.txt|grep -v 'Ratio' ) > $indir/${window}M/${window}M_Merge_GC.tmp.txt
		mv $indir/${window}M/${window}M_Merge_GC.tmp.txt $indir/${window}M/${window}M_Merge_GC.txt
	done
	echo -e -n "\n" >> $indir/${window}M/${window}M_read_header.txt
	cat $indir/${window}M/${window}M_read_header.txt $indir/${window}M/${window}M_Merge_read.txt >$indir/${window}M/${window}M_Merge_read_Final.txt
	mv $indir/${window}M/${window}M_Merge_read_Final.txt $indir/${window}M/${outfile}_${window}M_Merge_read.txt

	cat $indir/${window}M/${window}M_read_header.txt $indir/${window}M/${window}M_Merge_GC.txt >$indir/${window}M/${window}M_Merge_GC_Final.txt
	mv $indir/${window}M/${window}M_Merge_GC_Final.txt $indir/${window}M/${outfile}_${window}M_Merge_GC.txt


done


