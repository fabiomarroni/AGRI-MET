#!/bin/bash

###############################################################
#
#This script generates the results that in the paper are called small RNA^C and small RNA^A
#
###############################################################

###########################################################
#
# Training set
#
###########################################################



#Performed some quality trimming on small RNA since they may still contain adapter sequences
module load sw/aligners/erne/1.4.6

READDIR=smallRNA/00_reads
min_length=9
output_dir=smallRNA/01_trimming
for SUBDIR in 32_1 32_12 32_4
do
read1=${READDIR}/${SUBDIR}/01_trimming/*fastq.gz
prefix=$(echo $(basename $read1) | sed -e 's/.fastq.gz//g')
erne-filter \
--query1 ${read1} \
--min-size ${min_length} \
--threads 4 \
--gzip \
--output-prefix ${output_dir}/${prefix} > ${output_dir}/${prefix}_ernefilt.log
done


###############################################################
#
#Use of virusdetect pipeline, generates the results that in the paper are called small RNA^C
#
###############################################################

#I now try to use a pipeline specifically developed for analysing sRNA to identify viruses.
#As a first step I perform the raw analysis. 
module load lang/perl/5.10.1
module load sw/aligners/bwa/0.7.10
module load sw/bio/samtools/0.1.18
module load sw/bio/blast/2.2.27

#Set this local perl lib in which I have CGI.pm
#Users may need to tweak with their perl installation to have this work.
#As an alternative users can use the online version of virusdetect
#http://virusdetect.feilab.net/cgi-bin/virusdetect/index.cgi
export PERL5LIB=/home/user/.local/lib/perl5/:$PERL5LIB
cd smallRNA/02_VirusDetect
for MY_IN in smallRNA/01_trimming/*fastq
do
perl software/VirusDetect_v1.7/virus_detect.pl $MY_IN
done 

#Compact output tables.
OUT_DIR=smallRNA/03_smallRNA_tables/
mkdir -p $OUT_DIR
for MY_DIR in smallRNA/02_VirusDetect/*
do
echo $MY_DIR
MY_OUT=$(basename $MY_DIR)
MY_OUT=${MY_OUT/result_/}
MY_OUT=${MY_OUT/.trimmed_1.fastq/}
Rscript functions/01_small_RNA.r -I ${MY_DIR} -O ${OUT_DIR}/${MY_OUT}.txt \
	-S ${OUT_DIR}/${MY_OUT}_summary.txt -N reference/Kraken/taxonomy/names.dmp \
	-T reference/Kraken/seqid2taxid.map
done

#Merge small RNA results in a single table
PATTERN=_summary.txt
Rscript functions/01a_merge_virusdetect_prop.r -I ${OUT_DIR} -P ${PATTERN} -O tables/smallRNA_contigs_VirDet.txt


###############################################################
#
#End of use of virusdetect pipeline (small RNA^C)
#
###############################################################



#Align smallRNA against the same database against which we aligned RNA-seq, and use the same counting approaches
#This part produces the results that in the paper are called RNA^A

#We need R 3.5.1 for the last step
module purge
export PATH=$PATH:/anaconda-3/bin
source activate r_3.5.1
module load sw/aligners/hisat2/2.0.4
module load sw/bio/samtools/0.1.18
module load it/assemblers/stringtie/1.3.4d
NCORES=4
OUTDIR=smallRNA/05_alignments_virusdb_vv
TAXFILE=software/VirusDetect_v1.7/databases/vrl_genbank.info.gz
FULLGFF=reference/vv_vrl_plant.gff
INDEX=reference/vv_virusdetect
mkdir -p $OUTDIR
for MY_IN in smallRNA/01_trimming/*fastq
do
echo $MY_IN
OUT_NAME=$(basename $MY_IN | sed 's/_1.fastq/.sam/g')
echo $OUT_NAME
aaa=${OUTDIR}/$OUT_NAME
echo $aaa
hisat2 $INDEX --no-unal -U $MY_IN -S $aaa --threads $NCORES
bbb=${aaa/.sam/.bam}
echo "Converting $aaa to $bbb"
samtools view -bS $aaa -o $bbb
ccc=${bbb/.bam/_sorted}
samtools sort $bbb $ccc
echo "indexing $ccc.bam"
samtools index ${ccc}.bam
MYCOUNTS=${ccc/_sorted/.tab}
MYGFF=${ccc/_sorted/.gff3}
echo $MYGFF
stringtie $ccc.bam -v -p $NCORES -e -G $FULLGFF -A $MYCOUNTS -o $MYGFF
OUTFILE=${MYGFF/.gff3/_TPM_FPKM.txt}
Rscript functions/06_virusdetect_TPM_FPKM.r -G ${MYGFF} -T ${TAXFILE} -O ${OUTFILE}
done 

OUTDIR=smallRNA/2018/05_alignments_virusdb_vv
OUTFILE=tables/smallRNAseq_VirDet.txt
#Merge results in a single table
Rscript functions/07_merge_virusdetect_TPM_FPKM.r -I ${OUTDIR} -O ${OUTFILE}






###########################################################
#
# Test (validation) set
#
###########################################################


#I performed quality trimming using erne.
INDIR=smallRNA_val/00_reads/
min_length=9
output_dir=smallRNA_val/01_trimming
for read1 in smallRNA_val/00_reads/*fastq.gz
do
echo $read1
prefix=$(echo $(basename $read1) | sed -e 's/.fastq.gz//g')
echo "module load sw/aligners/erne/1.4.6; erne-filter \
--query1 ${read1} \
--min-size ${min_length} \
--threads 4 \
--gzip \
--output-prefix ${output_dir}/${prefix} > ${output_dir}/${prefix}_ernefilt.log" | qsub -N erne_$read1 -l vmem=32G,walltime=24:00:00,nodes=1:ppn=16
done



###############################################################
#
#Use of virusdetect pipeline
#
###############################################################

#I now try to use a pipeline specifically developed for analysing sRNA to identify viruses.
#As a first step I perform the raw analysis. I might then decide to repeat the analysis after removing rRNA.
module load lang/perl/5.10.1
module load sw/aligners/bwa/0.7.10
module load sw/bio/samtools/0.1.18
module load sw/bio/blast/2.2.27
#module load it/aligners/blast/2.10.1+


#Set this local perl lib in which I have CGI.pm
export PERL5LIB=/home/marroni/.local/lib/perl5/:$PERL5LIB
cd smallRNA_val/02_VirusDetect
for MY_IN in smallRNA_val/01_trimming/*.fastq
do
ERRLOG=$(basename $MY_IN)
echo $ERRLOG
perl software/VirusDetect/virus_detect.pl $MY_IN 2>${ERRLOG/fastq/log}
done 


#Compact output tables.
OUT_DIR=smallRNA_val/03_smallRNA_tables
mkdir -p $OUT_DIR
for MY_DIR in smallRNA_val/02_VirusDetect/result*
do
echo $MY_DIR
MY_OUT=$(basename $MY_DIR)
MY_OUT=${MY_OUT/_1.fastq_temp/}
Rscript functions/01_small_RNA.r -I ${MY_DIR} -O ${OUT_DIR}/${MY_OUT}.txt \
-S ${OUT_DIR}/${MY_OUT}_summary.txt -N reference/Kraken/taxonomy/names.dmp \
	-T reference/Kraken/seqid2taxid.map
done


#Merge small RNA results in a single table
PATTERN=_summary.txt
Rscript functions/01b_merge_virusdetect_prop_val.r -I ${OUT_DIR} -P ${PATTERN} -O tables_val/smallRNA_contigs_VirDet.txt



###############################################################
#
#End of use of virusdetect pipeline
#
###############################################################



#Align smallRNA against the same database against which we aligned RNA-seq, and use the same counting approaches


#We need R 3.5.1 for the last step
module purge
export PATH=$PATH:/anaconda-3/bin
source activate r_3.5.1
module load sw/aligners/hisat2/2.0.4
module load sw/bio/samtools/0.1.18
module load it/assemblers/stringtie/1.3.4d
NCORES=4
OUTDIR=smallRNA_val/05_alignments_virusdb_vv
TAXFILE=software/VirusDetect_v1.7/databases/vrl_genbank.info.gz
FULLGFF=reference/vv_vrl_plant.gff
INDEX=reference/vv_virusdetect
mkdir -p $OUTDIR
for MY_IN in smallRNA_val/01_trimming/*fastq
do
echo $MY_IN
OUT_NAME=$(basename $MY_IN | sed 's/_1.fastq/.sam/g')
echo $OUT_NAME
aaa=${OUTDIR}/$OUT_NAME
echo $aaa
hisat2 $INDEX --no-unal -U $MY_IN -S $aaa --threads $NCORES
bbb=${aaa/.sam/.bam}
echo "Converting $aaa to $bbb"
samtools view -bS $aaa -o $bbb
ccc=${bbb/.bam/_sorted}
samtools sort $bbb $ccc
echo "indexing $ccc.bam"
samtools index ${ccc}.bam
MYCOUNTS=${ccc/_sorted/.tab}
MYGFF=${ccc/_sorted/.gff3}
echo $MYGFF
stringtie $ccc.bam -v -p $NCORES -e -G $FULLGFF -A $MYCOUNTS -o $MYGFF
OUTFILE=${MYGFF/.gff3/_TPM_FPKM.txt}
Rscript functions/06_virusdetect_TPM_FPKM.r -G ${MYGFF} -T ${TAXFILE} -O ${OUTFILE}
done 



OUTDIR=smallRNA_val/05_alignments_virusdb_vv
OUTFILE=tables_val/smallRNAseq_align.txt
#Merge results in a single table
Rscript functions/07_merge_virusdetect_TPM_FPKM.r -I ${OUTDIR} -O ${OUTFILE}


