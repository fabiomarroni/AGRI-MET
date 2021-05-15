#!/bin/bash
#This script generates Tables S4-S7

#!/bin/bash

#RNA a
Rscript --no-save functions/09_allvirtable.r \
  --infile1 tables/RNAseq_VirDet.txt \
  --infile2 tables_val/val_RNAseq_VirDet.txt \
  --namesfile Docs_val/Sample_names.txt \
  --threshold 19.28 \
  --namecol Samplename \
  --outfile tables_val/RNA_a.xlsx 

#RNA b
Rscript --no-save functions/09_allvirtable.r \
  --infile1 tables/RNA_bracken_VD_206.txt \
  --infile2 tables_val/val_RNA_bracken_VD_206.txt \
  --namesfile Docs_val/Sample_names.txt \
  --threshold 385.81 \
  --namecol Bracken_name \
  --outfile tables_val/RNA_b.xlsx 

#Small RNA a
Rscript --no-save functions/09_allvirtable.r \
  --infile1 tables/smallRNAseq_VirDet.txt \
  --infile2 tables_val/smallRNAseq_align.txt \
  --namesfile Docs_val/Sample_names.txt \
  --threshold 1572.03 \
  --namecol smallalign_name \
  --outfile tables_val/smallRNA_a.xlsx 

#Small RNA c
Rscript --no-save functions/09_allvirtable.r \
  --infile1 tables/smallRNA_contigs_VirDet.txt \
  --infile2 tables_val/smallRNA_contigs_VirDet.txt \
  --namesfile Docs_val/Sample_names.txt \
  --threshold 0.97 \
  --namecol smallcontig_name \
  --outfile tables_val/smallRNA_c.xlsx 

