#!/bin/bash
#This script generates graphically "appealing" versions of Tables 6-9 of the paper
#It laso takes care of converting the internal names in the official names

#Heatmaps

#RNAa
Rscript --no-save functions/AGRI-MET/09_plot_heatmaps.r \
  --infile tables/RNAseq_VirDet.txt \
  --threshold 19.28 \
  --orignames cabf_neg_FPKM,cort_gva_merged_FPKM,dolc_rl31_FPKM,malb_lr1_FPKM,maln_lr3_FPKM,mont_neg_merged_FPKM,ries_ara_merged_FPKM \
  --newnames CF,CO,DO,MB,MN,MO,RI \
  --outfile plots/heat_RNAseq_VirDet.pdf 

#RNAb
Rscript --no-save functions/AGRI-MET/09_plot_heatmaps.r \
  --infile tables/RNA_bracken_VD_206.txt \
  --threshold 385.81 \
  --orignames RPM_cabf_neg,RPM_cort_gva_merged,RPM_dolc_rl31,RPM_malb_lr1,RPM_maln_lr3,RPM_mont_neg_merged,RPM_ries_ara_merged \
  --newnames CF,CO,DO,MB,MN,MO,RI \
  --outfile plots/heat_RNA_bracken_VD_206.pdf 

#small-RNA contigs
Rscript --no-save functions/AGRI-MET/09_plot_heatmaps.r \
  --infile tables/smallRNA_contigs_VirDet.txt \
  --threshold 0.97 \
  --orignames Percent_mont_neg_merged,Percent_maln_lr3,Percent_cort_gva_merged \
  --newnames MO,MN,CO \
  --outfile plots/heat_smallRNA_contigs_VirDet.pdf 

#small-RNA alignments
Rscript --no-save functions/AGRI-MET/09_plot_heatmaps.r \
  --infile tables/smallRNAseq_VirDet.txt \
  --threshold 1572.03 \
  --orignames mont_neg_merged,maln_lr3,cort_gva_merged \
  --newnames MO,MN,CO \
  --outfile plots/heat_smallRNAseq_VirDet.pdf 

#Tables
#RNAa
Rscript --no-save functions/AGRI-MET/09b_draw_tables.r \
  --infile tables/RNAseq_VirDet.txt \
  --threshold 19.28 \
  --orignames cabf_neg_FPKM,cort_gva_merged_FPKM,dolc_rl31_FPKM,malb_lr1_FPKM,maln_lr3_FPKM,mont_neg_merged_FPKM,ries_ara_merged_FPKM \
  --newnames CF,U1,U2,MB,MN,MO,RI \
  --outfile plots/heat_RNAseq_VirDet.html 

#RNAb
Rscript --no-save functions/AGRI-MET/09b_draw_tables.r \
  --infile tables/RNA_bracken_VD_206.txt \
  --threshold 385.81 \
  --orignames RPM_cabf_neg,RPM_cort_gva_merged,RPM_dolc_rl31,RPM_malb_lr1,RPM_maln_lr3,RPM_mont_neg_merged,RPM_ries_ara_merged \
  --newnames CF,U1,U2,MB,MN,MO,RI \
  --outfile plots/heat_RNA_bracken_VD_206.html 

#small-RNA contigs
Rscript --no-save functions/AGRI-MET/09b_draw_tables.r \
  --infile tables/smallRNA_contigs_VirDet.txt \
  --threshold 0.97 \
  --orignames Percent_mont_neg_merged,Percent_maln_lr3,Percent_cort_gva_merged \
  --newnames MO,MN,U1 \
  --outfile plots/heat_smallRNA_contigs_VirDet.html 

#small-RNA alignments
Rscript --no-save functions/AGRI-MET/09b_draw_tables.r \
  --infile tables/smallRNAseq_VirDet.txt \
  --threshold 1572.03 \
  --orignames mont_neg_merged,maln_lr3,cort_gva_merged \
  --newnames MO,MN,U1 \
  --outfile plots/heat_smallRNAseq_VirDet.html 
