#!/bin/bash



#Thresholds, as defined by results obtained on the training set
RNAa=19.28
RNAb=386
smallrnaa=1572
smallrnac=0.97

#RT-PCR file (the new one is the one by Luigi Falginella in Rauscedo)
RTPCR=Docs_val/Val_RT_Actin.csv
#Elisa file
ELISA=Docs_val/Val_ELISA.csv


#This is the function to perform the validation against ELISA and RT-PCR
#Additional parameters are specified in the R function
#use both our gold standard combination procedures

Rscript functions/12_val_RT_ELISA.r \
-B tables_val/val_RNA_bracken_VD_206.txt \
-P Docs_val/Val_RT_Actin.csv \
-E Docs_val/Val_ELISA.csv \
-R tables_val/val_RNAseq_VirDet.txt \
-S tables_val/smallRNAseq_align.txt \
-C tables_val/smallRNA_contigs_VirDet.txt \
-N CO-U1-1-1,DO-U2-1-1,MB-5-1,MN-1-1,RI-6-1 \
-O tables_val/ \
-g either

Rscript functions/12_val_RT_ELISA.r \
-B tables_val/val_RNA_bracken_VD_206.txt \
-P Docs_val/Val_RT_Actin.csv \
-E Docs_val/Val_ELISA.csv \
-R tables_val/val_RNAseq_VirDet.txt \
-S tables_val/smallRNAseq_align.txt \
-C tables_val/smallRNA_contigs_VirDet.txt \
-N CO-U1-1-1,DO-U2-1-1,MB-5-1,MN-1-1,RI-6-1 \
-O tables_val/ \
-g concordant



#Compute summary statistics for validation
TABDIR=tables_val
for GS in either concordant
do
Rscript functions/13_explore_validation.r \
--remove.failed TRUE \
--virus core \
--brackenfile ${TABDIR}/Val_RNAb_${GS}.txt \
--rnaseqfile ${TABDIR}/Val_RNAa_${GS}.txt \
--smallrnafile ${TABDIR}/Val_smalla_${GS}.txt \
--smallrnacontigfile ${TABDIR}/Val_smallc_${GS}.txt \
-O ${TABDIR}/final_stats_core_${GS}.txt
done



#Some plants were analyzed in both Test and Validation set. This make me think that we shouldn't treat the validation set as completely independent
#Also, we can check technical replicability for all 4 tests.

Rscript functions/14_replicability.r \
-O tables_techval/

Rscript functions/15_val_newvir.r


