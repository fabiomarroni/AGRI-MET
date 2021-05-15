#!/bin/bash



###########################################################
#
# Training set
#
###########################################################


#COUNT RNAseq pairs of reads
unset MYSTATS
declare -a MYSTATS
RNADIR=RNAseq/01_reads
cd $RNADIR
for aaa in *_1.fastq.gz
do
MYLINES=$(zcat $aaa | wc -l)
MYREADS=$(echo $MYLINES/2 | bc)
MYSTATS+=($aaa $MYREADS)
done 
cd ../..
echo ${MYSTATS[@]} > tables/RNAseq_reads.txt

#COUNT smallRNA reads
unset MYSTATS
declare -a MYSTATS
RNADIR=smallRNA/01_trimming
cd $RNADIR
for aaa in *_1.fastq
do
MYLINES=$(cat $aaa | wc -l)
MYREADS=$(echo $MYLINES/4 | bc)
MYSTATS+=($aaa $MYREADS)
done 
cd ../..
echo ${MYSTATS[@]} > tables/smallRNA_reads.txt

#Count TPM and FPKM on viruses for total RNAa (used for Table 1)
ALDIR=RNAseq/02a_alignments_virusdb_vv
unset MYVIRT
unset MYVIRF
declare -a MYVIRT
declare -a MYVIRF
cd $ALDIR
for aaa in *_TPM_FPKM.txt
do
FILENAME=$(basename $aaa)
VIRALT=$(awk -F '\t' '{sum += $9} END {print sum}' $aaa)
MYVIRT+=($FILENAME,$VIRALT)
VIRALF=$(awk -F '\t' '{sum += $8} END {print sum}' $aaa)
MYVIRF+=($FILENAME,$VIRALF)
done
cd ../..
echo ${MYVIRT[@]} > tables/totalRNA_viral_TPM.txt
echo ${MYVIRF[@]} > tables/totalRNA_viral_FPKM.txt

#Count TPM and FPKM on viruses for small RNA (used for Table 1)
ALDIR=smallRNA/05_alignments_virusdb_vv
unset MYVIRT
unset MYVIRF
declare -a MYVIRT
declare -a MYVIRF
cd $ALDIR
for aaa in *_TPM_FPKM.txt
do
FILENAME=$(basename $aaa)
VIRALT=$(awk -F '\t' '{sum += $9} END {print sum}' $aaa)
MYVIRT+=($FILENAME,$VIRALT)
VIRALF=$(awk -F '\t' '{sum += $8} END {print sum}' $aaa)
MYVIRF+=($FILENAME,$VIRALF)
done
cd ../..
echo ${MYVIRT[@]} > tables/smallRNA_viral_TPM.txt
echo ${MYVIRF[@]} > tables/smallRNA_viral_FPKM.txt


###########################################################
#
# Test (validation) set
#
###########################################################


#COUNT RNAseq reads
unset MYSTATS
declare -a MYSTATS
RNADIR=RNAseq_val/01_reads
cd $RNADIR
for aaa in *R1_001.fastq.gz
do
MYLINES=$(zcat $aaa | wc -l)
MYREADS=$(echo $MYLINES/2 | bc)
MYSTATS+=($aaa,$MYREADS)
done 
cd ../..
printf '%s\n' "${MYSTATS[@]}" > tables_val/valid_set_RNAseq_reads.txt

#COUNT smallRNA reads
unset MYSTATS
declare -a MYSTATS
RNADIR=smallRNA_val/01_trimming
cd $RNADIR
for aaa in *_1.fastq
do
MYLINES=$(cat $aaa | wc -l)
MYREADS=$(echo $MYLINES/4 | bc)
MYSTATS+=($aaa $MYREADS)
done 
cd ../..
echo ${MYSTATS[@]} > tables_val/valid_set_smallRNA_reads.txt


#Count TPM and FPKM on viruses for total RNA
ALDIR=RNAseq_val/02a_alignments_virusdb_vv_val
unset MYVIRT
unset MYVIRF
declare -a MYVIRT
declare -a MYVIRF
cd $ALDIR
for aaa in *_TPM_FPKM.txt
do
FILENAME=$(basename $aaa)
VIRALT=$(awk -F '\t' '{sum += $9} END {print sum}' $aaa)
MYVIRT+=($FILENAME,$VIRALT)
VIRALF=$(awk -F '\t' '{sum += $8} END {print sum}' $aaa)
MYVIRF+=($FILENAME,$VIRALF)
done
cd ../..
echo ${MYVIRT[@]} > tables_val/totalRNA_viral_TPM.txt
echo ${MYVIRF[@]} > tables_val/totalRNA_viral_FPKM.txt

#Count TPM on viruses for small RNA
ALDIR=smallRNA_val/05_alignments_virusdb_vv
unset MYVIRT
unset MYVIRF
declare -a MYVIRT
declare -a MYVIRF
cd $ALDIR
for aaa in *_TPM_FPKM.txt
do
FILENAME=$(basename $aaa)
VIRALT=$(awk -F '\t' '{sum += $9} END {print sum}' $aaa)
MYVIRT+=($FILENAME,$VIRALT)
VIRALF=$(awk -F '\t' '{sum += $8} END {print sum}' $aaa)
MYVIRF+=($FILENAME,$VIRALF)
done
cd ../..
echo ${MYVIRT[@]} > tables_val/smallRNA_viral_TPM.txt
echo ${MYVIRF[@]} > smallRNA_viral_FPKM.txt

