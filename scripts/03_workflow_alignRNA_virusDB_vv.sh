#!/bin/bash

#############################################################ààà
#
# This script generates the RNA^A results mentioned in the paper
#
#############################################################ààà


#Align reads to the database used for smallRNA
module load sw/aligners/hisat2/2.0.4
module load sw/bio/samtools/0.1.18

#We assume that you have vitis vinifera reference sequence and the corresponding annotation gff


#Generate index for vitis and virusdetect viruses. 
#I use this one because the refseq database didn't have grapevine virus A, which we need, and I think I trust more the virusdetect database. 
hisat2-build vitis_vinifera_12xCHR.fasta,VirusDetect_v1.7/databases/vrl_plant reference/vv_virusdetect


#Prepare to have a gff of all our sequences. Gff for viruses will be "fake", but it's ok.

samtools faidx software/VirusDetect_v1.7/databases/vrl_plant 

cd reference
mv VirusDetect_v1.7/databases/vrl_plant.fai .

FAIDX=reference/vrl_plant.fai
MYGFF=reference/vrl_plant.gff
Rscript functions/05_faidx_to_gff.r -X ${FAIDX} -T ${MYGFF}

FULLGFF=reference/vv_vrl_plant.gff
cat /vitis/genes/annotation/V2.1/V2.1.gff3 $MYGFF > $FULLGFF
#Version for htseq-count
sed -e 's/ID=/gene_id=/g' $FULLGFF > ${FULLGFF/plant.gff/htseq_plant.gff}


###########################################################
#
# Training set
#
###########################################################

INDEX=reference/vv_virusdetect
READDIR=RNAseq/01_reads

NCORES=4
OUTDIR=02a_alignments_virusdb_vv
mkdir -p $OUTDIR
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do
echo $MY_IN1
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz}
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz/.sam/g')
echo $OUT_NAME
aaa=${OUTDIR}/$OUT_NAME
aln=`echo "module load sw/aligners/hisat2/2.0.4; \
hisat2 $INDEX \
--no-unal -1 $MY_IN1 -2 $MY_IN2 -S $aaa --threads $NCORES" | qsub -N aln_$OUT_NAME -l vmem=10G,walltime=24:00:00,nodes=1:ppn=$NCORES`
bbb=${aaa/.sam/.bam}
#echo "Converting $aaa to $bbb"
cbam=`echo "module load sw/bio/samtools/0.1.18; samtools view -bS $aaa -o $bbb" | \
qsub -N cbam_$OUT_NAME -l vmem=8G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$aln`
#rm -f ${aaa/.bam/.sam}
ccc=${bbb/.bam/_sorted}
#echo "Sorting $bbb to $ccc.bam"
sbam=`echo "module load sw/bio/samtools/0.1.18; samtools sort $bbb $ccc" | \
qsub -N cbam_$OUT_NAME -l vmem=16G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$cbam`
#echo "indexing $ccc.bam"
ibam=`echo "module load sw/bio/samtools/0.1.18; samtools index ${ccc}.bam" | \
qsub -N cbam_$OUT_NAME -l vmem=8G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$sbam`
#samtools idxstats ${ccc}.bam | cut -f1,3 > ${ccc}_stats.out
done 



NCORES=4
for aaa in ${OUTDIR}/*_sorted.bam
do
MYCOUNTS=${aaa/_sorted.bam/.tab}
MYGFF=${aaa/_sorted.bam/.gff3}
CC=$(basename $MYGFF)
echo $MYGFF
echo "module load it/assemblers/stringtie/1.3.4d; stringtie $aaa -v -p $NCORES -e -G $FULLGFF -A $MYCOUNTS -o $MYGFF" | \
qsub -N st_$CC -l vmem=10G,walltime=24:00:00,nodes=1:ppn=$NCORES
done

#Load R 3.5.1 on anaconda, because it's the only one able to "read" gzipped files. I will start usind that 
module purge
export PATH=$PATH:/anaconda-3/bin
source activate r_3.5.1
TAXFILE=VirusDetect_v1.7/databases/vrl_genbank.info.gz
for MYGFF in ${OUTDIR}/*gff3
do
OUTFILE=${MYGFF/.gff3/_TPM_FPKM.txt}
Rscript functions/06_virusdetect_TPM_FPKM.r -G ${MYGFF} -T ${TAXFILE} -O ${OUTFILE}
done

OUTDIR=02a_alignments_virusdb_vv
OUTFILE=tables/RNAseq_VirDet.txt
#Merge results in a single table
Rscript functions/07_merge_virusdetect_TPM_FPKM.r -I ${OUTDIR} -O ${OUTFILE}


#Plot RT-PCR versus RNAseq alignment data
module load it/lang/r/3.6.1
BRACKENFILE=tables/RNA_bracken_VD_206.txt
RTFILE=tables/RT-PCR.txt
ELISAFILE=tables/Elisa_res.txt
RNASEQFILE=tables/RNAseq_VirDet.txt
OUTFILE=tables/RNA_smallRNA_RTPCR_gold_std.txt
GRAPHDIR=plots/
CONSMALLRNA=tables/smallRNA_contigs_VirDet.txt
SMALLRNA=tables/smallRNAseq_VirDet.txt
Rscript functions/08_plot_RTPCR_RNAseq.r -B ${BRACKENFILE} -P ${RTFILE} -E ${ELISAFILE} -R ${RNASEQFILE} -O ${OUTFILE} -S ${SMALLRNA} -C ${CONSMALLRNA} -G ${GRAPHDIR} -I TRUE
Rscript functions/08_plot_RTPCR_RNAseq.r -B ${BRACKENFILE} -P ${RTFILE} -E ${ELISAFILE} -R ${RNASEQFILE} -O ${OUTFILE} -S ${SMALLRNA} -C ${CONSMALLRNA} -G ${GRAPHDIR} -I FALSE






###########################################################
#
# Test (validation) set
#
###########################################################

INDEX=reference/vv_virusdetect
READDIR=RNAseq_val/01_reads

NCORES=4
OUTDIR=RNAseq_val/02a_alignments_virusdb_vv_val
mkdir -p $OUTDIR
for MY_IN1 in ${READDIR}/*R1_001.fastq.gz
do
echo $MY_IN1
MY_IN2=${MY_IN1/R1_001.fastq.gz/R2_001.fastq.gz}
OUT_NAME=$(basename $MY_IN1 | sed 's/_R1_001.fastq.gz/.sam/g')
echo $OUT_NAME
aaa=${OUTDIR}/$OUT_NAME
aln=`echo "module load sw/aligners/hisat2/2.0.4; \
hisat2 $INDEX \
--no-unal -1 $MY_IN1 -2 $MY_IN2 -S $aaa --threads $NCORES" | qsub -N aln_$OUT_NAME -l vmem=10G,walltime=24:00:00,nodes=1:ppn=$NCORES`
bbb=${aaa/.sam/.bam}
#echo "Converting $aaa to $bbb"
cbam=`echo "module load sw/bio/samtools/0.1.18; samtools view -bS $aaa -o $bbb" | \
qsub -N cbam_$OUT_NAME -l vmem=8G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$aln`
#rm -f ${aaa/.bam/.sam}
ccc=${bbb/.bam/_sorted}
#echo "Sorting $bbb to $ccc.bam"
sbam=`echo "module load sw/bio/samtools/0.1.18; samtools sort $bbb $ccc" | \
qsub -N cbam_$OUT_NAME -l vmem=16G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$cbam`
#echo "indexing $ccc.bam"
ibam=`echo "module load sw/bio/samtools/0.1.18; samtools index ${ccc}.bam" | \
qsub -N cbam_$OUT_NAME -l vmem=8G,walltime=12:00:00,nodes=1:ppn=1 -W depend=afterok:$sbam`
#samtools idxstats ${ccc}.bam | cut -f1,3 > ${ccc}_stats.out
done 


FULLGFF=reference/vv_vrl_plant.gff
NCORES=4
for aaa in ${OUTDIR}/*_sorted.bam
do
MYCOUNTS=${aaa/_sorted.bam/.tab}
MYGFF=${aaa/_sorted.bam/.gff3}
CC=$(basename $MYGFF)
echo $MYGFF
echo "module load it/assemblers/stringtie/1.3.4d; stringtie $aaa -v -p $NCORES -e -G $FULLGFF -A $MYCOUNTS -o $MYGFF" | \
qsub -N st_$CC -l vmem=10G,walltime=24:00:00,nodes=1:ppn=$NCORES
done

#Summarize results
#Load R 3.5.1 on anaconda, because it's the only one able to "read" gzipped files. I will start using that 
module purge
export PATH=$PATH:/anaconda-3/bin
source activate r_3.5.1
TAXFILE=software/VirusDetect_v1.7/databases/vrl_genbank.info.gz
for MYGFF in ${OUTDIR}/*gff3
do
OUTFILE=${MYGFF/.gff3/_TPM_FPKM.txt}
Rscript functions/06_virusdetect_TPM_FPKM.r -G ${MYGFF} -T ${TAXFILE} -O ${OUTFILE}
done

OUTFILE=tables_val/val_RNAseq_VirDet.txt
#Merge results in a single table
Rscript functions/07_merge_virusdetect_TPM_FPKM.r -I ${OUTDIR} -O ${OUTFILE}

