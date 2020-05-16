#!/bin/bash
#COUNT RNAseq reads


#Folders 01_reads and/or 01_trimming are the folders in which raw and/or trimmed reads are.
#Reads (trimmed) can be downloaded from SRA, accession PRJNA628964


unset MYSTATS
declare -a MYSTATS
RNADIR=../RNAseq/01_reads
mkdir -p $RNADIR/tables
cd $RNADIR
for aaa in *_1.fastq.gz
do
MYLINES=$(zcat $aaa | wc -l)
MYREADS=$(echo $MYLINES/2 | bc)
MYSTATS+=($aaa $MYREADS)
done 
echo ${MYSTATS[@]} > tables/RNAseq_reads.txt

#COUNT smallRNA reads
unset MYSTATS
declare -a MYSTATS
RNADIR=../smallRNA/2018/01_trimming
mkdir -p $RNADIR/tables
cd $RNADIR
for aaa in *_1.fastq
do
MYLINES=$(cat $aaa | wc -l)
MYREADS=$(echo $MYLINES/4 | bc)
MYSTATS+=($aaa $MYREADS)
done 
echo ${MYSTATS[@]} > tables/smallRNA_reads.txt

