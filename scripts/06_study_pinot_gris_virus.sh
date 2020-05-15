#!/bin/bash

#This scripts perform all the analysis on the identification of GPGV strains.


GPGVDIR=GPGV
mkdir -p $GPGVDIR/reference

#I (manually) entered the GPGV accessions studied by Giulia in her 2019 paper in the perl script and then I retrieve the fasta
#with eutils. 
module load lang/perl/5.10.1
perl scripts/AGRI-MET/use_eutils_tarquini_onlyfasta.pl > $GPGVDIR/reference/ref_tarquini.fasta

#I also downloaded the reference assembly (found here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000894735.2/), and will paste to Giulia's list.
cat $GPGVDIR/reference/ref_tarquini.fasta $GPGVDIR/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.fna > $GPGVDIR/reference/newref_longnames_tarquini.fasta
cut -d" " -f1 $GPGVDIR/reference/newref_longnames_tarquini.fasta > $GPGVDIR/reference/newref_tarquini.fasta

#Explore the sequences
module load sw/aligners/clustalw2/2.0.10
#Since clustalw2 is still interactive, I wrote the parameters in a file.
clustalw2 < scripts/AGRI-MET/clustalw.par

#My idea is then to algn against the reference genome, to build a new reference strain for each grape based on SNPs 
#and then do a phylogenetic analysis of the strains to see if we can assign our viruses to one strain

#GATK recommends using STAR, so we do

#NOTE! Some minor parts of the script for running GATK may be missing due to a system crash which prevented saving the script
#Should be minor things, related to creating a single fasta reference out of the two (beacuse STAR can use two separate fasta, but the 
#stupid picard and gatk cannot!)
#Maybe also some indexing of the reference and creating the dictionary file are missing.

# Build a gtf including both vitis and GPGV genes
cat vitis_vinifera/genes/assembly_build_12xCHR/V2.1/V2.1.gtf $GPGVDIR/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.gtf > $GPGVDIR/reference/vv_gpgv.gtf

#Generate STAR genome indices
module load it/aligners/star/2.7.2b

STAR --runThreadN 8 \
	--runMode genomeGenerate \
	--genomeDir $GPGVDIR/reference \
	--genomeFastaFiles vitis_vinifera_12xCHR.fasta $GPGVDIR/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.fna \
	--sjdbGTFfile $GPGVDIR/reference/vv_gpgv.gtf

#Align using STAR	
READDIR=RNAseq/01_reads
STARDIR=$GPGVDIR/star/
mkdir -p $STARDIR

for MY_IN1 in ${READDIR}/*1.fastq.gz 
do
echo $MY_IN1
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz}
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g')
echo $OUT_NAME
STAR --runThreadN 8 \
	--genomeDir $GPGVDIR/reference \
	--readFilesIn $MY_IN1 $MY_IN2 \
	--readFilesCommand gunzip -c \
	--quantMode GeneCounts \
	--twopassMode Basic \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix $STARDIR/$OUT_NAME
done


#Start using GATK
#mark duplicates and split cigars
for MY_IN1 in ${READDIR}/*1.fastq.gz
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz} 
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}Aligned.sortedByCoord.out.bam
echo $INBAM 
MDOUT=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.bam} 
echo $MDOUT; MDMET=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.txt} 
java -jar /iga/scripts/packages/picard-tools-2.17.11/bin/picard.jar MarkDuplicates \
      I=$INBAM \
      O=$MDOUT \
      M=$MDMET
 SCOUT=${MDOUT/_01_md.bam/_02_sc.bam}
 gatk SplitNCigarReads       -R $GPGVDIR/reference/vv_gpgv.fasta       -I $MDOUT       -O $SCOUT; VCFOUT=${SCOUT/_02_sc.bam/.vcf.gz}; gatk --java-options "-Xmx4g" HaplotypeCaller     -R $GPGVDIR/reference/vv_gpgv.fasta    -I $SCOUT    -O $VCFOUT    -ERC GVCF
done

#add rg 
#Assume picard.jar is in you path
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz}
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}Aligned.sortedByCoord.out.bam
echo $INBAM; MDOUT=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.bam}
echo $MDOUT
MDMET=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.txt} 
SCOUT=${MDOUT/_01_md.bam/_02_sc.bam}
RGOUT=${MDOUT/_01_md.bam/_03_rg.bam}
java -jar picard.jar AddOrReplaceReadGroups \
      I=$SCOUT \
      O=$RGOUT \
	  RGID=1 \
      RGLB=F1 \
      RGPL=illumina \
      RGPU=ilmn \
      RGSM=$OUT_NAME 
done

#We need samtools to index bam
module load sw/bio/samtools/0.1.18

#Run haplotype caller also on vitis data (for Gabbo). 
####################################################
#GABBO ONLY!!!!!
##################################################
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz} 
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}_03_rg.bam
GVCFOUT=${INBAM/_03_rg.bam/_vv.gvcf.gz}
VCFOUT=${INBAM/_03_rg.bam/_vv.vcf.gz}
echo $INBAM
gatk=`echo "module load it/variants/gatk/4.1.4.0; gatk --java-options "-Xmx4g" HaplotypeCaller -R $GPGVDIR/reference/vv_gpgv.fasta -I $INBAM \
-O $GVCFOUT \
-ERC GVCF" | qsub -N g1 -l vmem=12G,walltime=8:00:00,nodes=1:ppn=1`
echo "module load it/variants/gatk/4.1.4.0; gatk --java-options "-Xmx4g" GenotypeGVCFs \
-R $GPGVDIR/reference/vv_gpgv.fasta \
-V $GVCFOUT \
-O $VCFOUT" | qsub -N g2 -l vmem=12G,walltime=8:00:00,nodes=1:ppn=1 -W depend=afterok:$gatk
done

###############################################################
##END OF GABBO ONLY!!!!!
###############################################################


#Run haplotype caller only on the viral sequences (because running also on vitis chromosomes is useless and takes forever!)
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz} 
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}Aligned.sortedByCoord.out.bam
MDOUT=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.bam}
MDMET=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.txt}
SCOUT=${MDOUT/_01_md.bam/_03_rg.bam} 
SMOUT=${MDOUT/_01_md.bam/_04_sm.bam} 
#VCFOUT=${SCOUT/_03_rg.bam/.vcf.gz}
VCFOUT=${SMOUT/_04_sm.bam/_04_sm.vcf.gz}
echo $SCOUT
samtools index $SCOUT
samtools view -b $SCOUT NC_015782.2 > $SMOUT
samtools index $SMOUT
gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R $GPGVDIR/reference/vv_gpgv.fasta \
	-I $SMOUT \
	-O $VCFOUT \
	-ERC GVCF
done


#Consolidate genotyping (we do that sample by sample, because we do not really care about the database)
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz} 
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}Aligned.sortedByCoord.out.bam
MDOUT=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.bam}
SMOUT=${MDOUT/_01_md.bam/_04_sm.bam} 
VCFIN=${SMOUT/_04_sm.bam/_04_sm.vcf.gz}
VCFOUT=${SMOUT/_04_sm.bam/_05_final.vcf.gz}
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R $GPGVDIR/reference/vv_gpgv.fasta \
   -V $VCFIN \
   -O $VCFOUT
done

# I try to use the small reference, hoping it will only return the small viral genome and not all the useless vitis
module load it/tools/picard-tools/2.17.11
samtools faidx GPGV/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.fna
java -jar picard.jar CreateSequenceDictionary \
	R=GPGV/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.fna \
	O=GPGV/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.dict
#Build alternative reference
for MY_IN1 in ${READDIR}/*1.fastq.gz 
do 
MY_IN2=${MY_IN1/1.fastq.gz/2.fastq.gz} 
OUT_NAME=$(basename $MY_IN1 | sed 's/_1.fastq.gz//g') 
INBAM=$STARDIR/${OUT_NAME}Aligned.sortedByCoord.out.bam
MDOUT=${INBAM/Aligned.sortedByCoord.out.bam/_01_md.bam}
SMOUT=${MDOUT/_01_md.bam/_04_sm.bam} 
VCFIN=${SMOUT/_04_sm.bam/_04_sm.vcf.gz}
VCFOUT=${SMOUT/_04_sm.bam/_05_final.vcf.gz}
gatk FastaAlternateReferenceMaker \
   -R $GPGVDIR/reference/ncbi-genomes-2020-01-19/GCF_000894735.2_ViralProj70003_genomic.fna \
   -O $GPGVDIR/star/GPGV_${OUT_NAME}.fasta \
   -L NC_015782.2 \
   -V $VCFOUT
done


for aaa in $GPGVDIR/star/GPGV_*fasta
do
filename=$(basename $aaa)
myname=${filename/GPGV_/}
myname=${myname/.fasta/}
echo $myname
sed -i "s/>.*/>${myname}/g" $aaa
done

#Merge RNAseq derived GPGV sequences in one file
cat $GPGVDIR/star/GPGV_*fasta $GPGVDIR/reference/newref_tarquini.fasta > $GPGVDIR/star/samples_and_ref_GPGV.fasta

#Merge RNAseq derived GPGV sequences (only for samples that are positive to GPGV) in one file together with Giulia's reference and samples
cat $GPGVDIR/star/GPGV_cort_gva_merged.fasta \
$GPGVDIR/star/GPGV_malb_lr1.fasta \
$GPGVDIR/star/GPGV_mont_neg_merged.fasta \
$GPGVDIR/star/GPGV_ries_ara_merged.fasta \
$GPGVDIR/reference/newref_tarquini.fasta > $GPGVDIR/star/pos_samples_and_ref_GPGV.fasta


#Explore the sequences
module load sw/aligners/clustalw2/2.0.10
#Use clustalw to align sequences
clustalw2 -INFILE=$GPGVDIR/star/pos_samples_and_ref_GPGV.fasta -ALIGN 

#Since some sequences are shorter, we delete the alignment parts at the beginning and at the end when there are a lot of gaps
#Otheriwse the gaps would be interpreted as phylogenetic divergence instead of difference in length.



Rscript functions/10_remove_gaps_from_alignments.r \
-I $GPGVDIR/star/pos_samples_and_ref_GPGV.aln \
-O $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.fasta
#Perform new alignment and build new phylogenetic tree without the errors due to different amplicon length
clustalw2 -INFILE=$GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.fasta -ALIGN 
clustalw2 -INFILE=$GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.fasta -TREE
clustalw2 -INFILE=$GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.fasta -BOOTSTRAP






#The bootstrap file is $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.phb
#I will substitute the sequence names to have a pretty output file.
#For Giulia's sequences I will use the names she used and not the accessions, and for our sequence I will use the 2 letter code.
#It's horrible to look, I know!!!
sed -i 's/MH087439.1/fvg-Is1/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087440.1/fvg-Is6/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087441.1/fvg-Is7/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087442.1/fvg-Is8/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087443.1/fvg-Is12/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087444.1/fvg-Is13/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087445.1/fvg-Is14/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087446.1/fvg-Is15/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/MH087447.1/fvg-Is17/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KR528581.1/Tannat/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KM491305.1/Mer/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KU194413.1/BC-1/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KX522755/Riesling/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/FR877530.1/IT/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KU312039.1/FEM01/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KF134124.1/SK01/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KF134125.1/SK13/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KF134123.1/SK30/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KF686810.1/SK30-1/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/KT894101.1/USA/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/malb_lr1/MB/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/ries_ara_merged/RI/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/mont_neg_merged/MO/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*
sed -i 's/cort_gva_merged/U1/g' $GPGVDIR/star/pos_samples_and_ref_nogap_GPGV.*

	



