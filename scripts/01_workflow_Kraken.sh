#!/bin/bash

#############################################################ààà
#
# This script generates the RNA^B results mentioned in the paper
#
#############################################################ààà


#Build kraken database based on VirusDetect 77.1 database
#Needed only once
#For full instructions, please refer to Kraken2 manual for custom database
#https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#custom-databases

#Module load is installation specific. You should retrieve and install blast on your HPC 
module load aligners/blast/latest #Needed for having dustmasker
KRAKBIN=software/kraken2-2.0.6-beta/bin
KDB=reference/Kraken
#Needed to have the taxonomy files
cp -r taxonomy $KDB

gunzip $KDB/taxonomy/nucl_gb.accession2taxid.gz

#Modifiy vitis database
VIT_TAXID=reference/vitis_vinifera_taxid.fasta
cp vitis_vinifera_12xCHR.fasta $VIT_TAXID

sed -i -e 's/>/>X75968 /g' $VIT_TAXID

#download file  https://github.com/kentnf/VirusDetect/blob/master/databases/vrl_plant to the folder reference
$KRAKBIN/kraken2-build --add-to-library reference/vrl_plant --db $KDB
$KRAKBIN/kraken2-build --add-to-library $VIT_TAXID --db $KDB
$KRAKBIN/kraken2-build --build --db $KDB --threads=8


#Build bracken kmer distribution (you will need to download bracken2)
BRACK_BIN=software/Bracken2
$BRACK_BIN/bracken-build -d ${KDB} -t 8 -k 35 -l 125 -x software/kraken2-2.0.6-beta/bin/
$BRACK_BIN/bracken-build -d ${KDB} -t 8 -k 35 -l 150 -x software/kraken2-2.0.6-beta/bin/

###########################################################
#
# Training set
#
###########################################################



#################################
#Run kraken 2.0.6
#################################
KDB=reference/Kraken
READDIR=RNAseq/01_reads
RESDIR=RNAseq/01a_kraken_VD206
NPROCS=8
mkdir -p $RESDIR
for READ1 in ${READDIR}/*_1.fastq.gz
do
READ2=${READ1/_1.fastq/_2.fastq}
SAMPLE=$(basename $READ1)
SAMPLE=${SAMPLE/_1.fastq.gz/}
#We assume you have a moab workload manager available. If not, just remove the echo and the part of command after the pipe "|"
echo "$KRAKBIN/kraken2 --db ${KDB} \
--threads $NPROCS \
--gzip-compressed --paired \
--output ${RESDIR}/${SAMPLE}.kraken \
--report ${RESDIR}/${SAMPLE}.kraken.report.txt \
${READ1} ${READ2} " | qsub -N krak1_${SAMPLE} -l vmem=210G,walltime=100:00:00,nodes=1:ppn=$NPROCS
done

#Bracken
BRACK_BIN=Bracken2
KDB=reference/Kraken
RESDIR=RNAseq/01a_kraken_VD206
THRESHOLD=10
READ_LEN=125
for CLASSIFICATION_LEVEL in S G
do
for READ1 in ${READDIR}/*_1.fastq.gz
do
read1=$(basename $READ1)
read1=${read1/_1.fastq.gz/}
echo $read1
KRES=$RESDIR
#Run bracken on kraken results
echo $THRESHOLD
echo "${BRACK_BIN}/bracken -d ${KDB} -i ${KRES}/${read1}.kraken.report.txt -o ${KRES}/${read1}_${CLASSIFICATION_LEVEL}.bracken.txt -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}" | qsub -N b_${read1} -l vmem=16G,walltime=24:00:00
done
done

MY_OUT=tables/RNA_bracken_VD_206.txt
Rscript functions/04a_summarize_bracken.r -I ${RESDIR} -O ${MY_OUT}


###########################################################
#
# Test (validation) set
#
###########################################################


#################################
#Run kraken 2.0.6
#################################
KRAKBIN=software/kraken2-2.0.6-beta/bin
KDB=reference/Kraken
READDIR=RNAseq_val/01_reads
RESDIR=RNAseq_val/01a_kraken_VD206
NPROCS=8
mkdir -p $RESDIR
for READ1 in ${READDIR}/*_R1_001.fastq.gz
do
READ2=${READ1/_R1_001.fastq/_R2_001.fastq}
SAMPLE=$(basename $READ1)
SAMPLE=${SAMPLE/_L001_R1_001.fastq.gz/}
echo "$KRAKBIN/kraken2 --db ${KDB} \
--threads $NPROCS \
--gzip-compressed --paired \
--output ${RESDIR}/${SAMPLE}.kraken \
--report ${RESDIR}/${SAMPLE}.kraken.report.txt \
${READ1} ${READ2} " | qsub -N krak1_${SAMPLE} -l vmem=210G,walltime=100:00:00,nodes=1:ppn=$NPROCS
done


#Bracken
BRACK_BIN=software/Bracken2
KDB=reference/Kraken
RESDIR=RNAseq_val/01a_kraken_VD206
THRESHOLD=10
#In the validation set we have reads of 150 and not 125, so I use the appropriate kmer distribution in bracken
READ_LEN=150
for CLASSIFICATION_LEVEL in S G
do
for READ1 in ${RESDIR}/*.kraken.report.txt
do
read1=$(basename $READ1)
read1=${read1/.kraken.report.txt/}
echo $read1
KRES=$RESDIR
#Run bracken on kraken results
echo $THRESHOLD
echo "${BRACK_BIN}/bracken -d ${KDB} -i ${KRES}/${read1}.kraken.report.txt -o ${KRES}/${read1}_${CLASSIFICATION_LEVEL}.bracken.txt -r ${READ_LEN} -l ${CLASSIFICATION_LEVEL} -t ${THRESHOLD}" | qsub -N b_${read1} -l vmem=16G,walltime=24:00:00
done
done

MY_OUT=tables_val/val_RNA_bracken_VD_206.txt
Rscript functions/04a_summarize_bracken.r -I ${RESDIR} -O ${MY_OUT}





