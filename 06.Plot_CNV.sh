#!/bin/bash
for window in 5M 10M
do

	Rscript /Share/BP/wangrui/bin/CNV/bin/CNV_Detected/04.CNV_Plot_FreeC.R  /WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV/${window}/${window}_Merge_GC.txt AllSample.txt  ${window}_NSCLC_FreeC.pdf /WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV/$window/

	Rscript /Share/BP/wangrui/bin/CNV/bin/CNV_Detected/04.CNV_ReadBase_Plot_WithControl.R /WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV/${window}/${window}_Merge_read.txt AllSample.txt ${window}_NSCLC_Control.pdf /WPSnew/wangrui/Project/Lung_Cancer/NSCLC/Bulk_DNA/03.CNV/$window/
done

