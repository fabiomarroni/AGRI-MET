# Run with --help flag for help.
# Modified 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-B", "--brackenfile"), type="character", default="tables_val/val_RNA_bracken_VD_206.txt",
              help="Bracken summary file [default= %default]", metavar="character"),
  make_option(c("-P", "--pcrfile"), type="character", default="Docs_val/Val_RT_Actin.csv",
              help="RT-PCR file [default= %default]", metavar="character"),
  make_option(c("-E", "--elisafile"), type="character", default="Docs_val/Val_ELISA.csv",
              help="ELISA file [default= %default]", metavar="character"),
  make_option(c("-R", "--rnaseqfile"), type="character", default="tables_val/val_RNAseq_VirDet.txt",
              help="Fasta index file [default= %default]", metavar="character"),
  make_option(c("-S", "--smallrnafile"), type="character", default="tables_val/smallRNAseq_align.txt",
              help="Small RNA file [default= %default]", metavar="character"),
  make_option(c("-C", "--smallrnacontigfile"), type="character", default="tables_val/smallRNA_contigs_VirDet.txt",
              help="Small RNA file [default= %default]", metavar="character"),
  make_option(c("-g", "--gold_standard"), type="character", default="concordant",
              help="Type of gold standard. Possible choices are 'either' and 'concordant' [default= %default]", metavar="character"),
  make_option(c("-N", "--newgeno"), type="character", default="CO-U1-1-1,DO-U2-1-1,MB-5-1,MN-1-1,RI-6-1",
              help="Comma delimited list of short names of the new genotypes, which can be used in the test set [default= %default]", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="tables_val/",
              help="Output table [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$brackenfile)) {
  stop("WARNING: No brackenfile specified with '-B' flag.")
} else {  cat ("Bracken summary file is ", opt$brackenfile, "\n")
  brackenfile <- opt$brackenfile  
  }

if (is.null(opt$pcrfile)) {
  stop("WARNING: No pcrfile specified with '-P' flag.")
} else {  cat ("PCR file is ", opt$pcrfile, "\n")
  pcrfile <- opt$pcrfile  
  }

  if (is.null(opt$elisafile)) {
  stop("WARNING: No elisafile specified with '-E' flag.")
} else {  cat ("ELISA file is ", opt$elisafile, "\n")
  elisafile <- opt$elisafile  
  }

if (is.null(opt$rnaseqfile)) {
  stop("WARNING: No rnaseqfile specified with '-R' flag.")
} else {  cat ("RNA seq file is ", opt$rnaseqfile, "\n")
  rnaseqfile <- opt$rnaseqfile  
  }

if (is.null(opt$smallrnafile)) {
  stop("WARNING: No smallrnafile specified with '-S' flag.")
} else {  cat ("smallrnafile file is ", opt$smallrnafile, "\n")
  smallrnafile <- opt$smallrnafile  
  }

if (is.null(opt$smallrnacontigfile)) {
  stop("WARNING: No smallrnacontigfile specified with '-C' flag.")
} else {  cat ("smallrnacontigfile file is ", opt$smallrnacontigfile, "\n")
  smallrnacontigfile <- opt$smallrnacontigfile  
  }

  if (is.null(opt$newgeno)) {
  stop("WARNING: No newgeno specified with '-N' flag.")
} else {  cat ("newgeno is ", opt$newgeno, "\n")
  newgeno <- opt$newgeno  
  }

  if (is.null(opt$gold_standard)) {
  stop("WARNING: No gold_standard specified with '-g' flag.")
} else {  cat ("gold_standard is ", opt$gold_standard, "\n")
  gold_standard <- opt$gold_standard  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outdir is ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }


  
 
RTPCR_RNAseq<-function(pcrfile,elisafile,rnaseqfile,smallrnafile,smallrnacontigfile,brackenfile,outdir,newgeno,
						gold_standard,RTthr=0.001,rnathr=19.28,brackthr=386,contigthr=0.97,salignthr=1572)
{
library("data.table")
library(scales)
library(RColorBrewer)
#Read the genotype names that can be used as a real test set.
valgeno<-unlist(strsplit(newgeno,","))
#We assume that the files have the same number of columns and with the same name
#Read results of RTPCR 
PCRdat<-fread(pcrfile,data.table=F,fill=T)
names(PCRdat)<-gsub("\\\xa0","",names(PCRdat))  #Solve a stupid formatting problem
PCRdat$ID<-PCRdat$var<-NULL
rownames(PCRdat)<-PCRdat[,1]
PCRdat[,1]<-NULL
PCRdat<-t(PCRdat)
PCRdat<-PCRdat[sort(row.names(PCRdat)),]
PCRdat<-round(PCRdat,3)
#Only select samples independent from the first round
PCRdat<-PCRdat[,valgeno]

#Read ELISA results
ELISAdat<-fread(elisafile,data.table=F)
names(ELISAdat)<-gsub("\\\xa0","",names(ELISAdat))  #Solve a stupid formatting problem
ELISAdat$ID<-ELISAdat$var<-NULL
rownames(ELISAdat)<-ELISAdat[,1]
ELISAdat[,1]<-NULL
ELISAdat<-t(ELISAdat)
ELISAdat<-ELISAdat[sort(row.names(ELISAdat)),]
#NaN are not positives
ELISAdat[ELISAdat=="NaN"]<-"-"
#Remove notes from ELISA names
colnames(ELISAdat)<-gsub("-\\(CO-sulla-falcon)","",colnames(ELISAdat))
colnames(ELISAdat)<-gsub("-\\(D01-sulla-falcon)","",colnames(ELISAdat))
#Only select samples independent from the first round
ELISAdat<-ELISAdat[,valgeno]
#Build Gold standard (GS) table: as we did in first round, a sample is positive if at least one of ELISA and RT is positive
allVir<-unique(c(row.names(PCRdat),row.names(ELISAdat)))
#Loop over viruses and samples and gather info on status
if(gold_standard=="either")
{
	GS<-matrix(0,nrow=length(allVir),ncol=ncol(ELISAdat),dimnames=list(allVir,colnames(ELISAdat)))
	for(vir in row.names(GS))
	{
		for(samp in colnames(GS))
		{
			#If RT is positive we write down and go to next
			if(PCRdat[vir,samp]>=RTthr) 
			{
				GS[vir,samp]<-1
				next
			}
			#Check if the virus has been sampled in ELISA
			if(!vir%in%row.names(ELISAdat)) next
			if(ELISAdat[vir,samp]=="+") GS[vir,samp]<-1
		}
	}
}


if(gold_standard=="concordant")
{
	GS<-matrix(NA,nrow=length(allVir),ncol=ncol(ELISAdat),dimnames=list(allVir,colnames(ELISAdat)))
	for(vir in row.names(GS))
	{
		if(!vir%in%row.names(ELISAdat)) next
		for(samp in colnames(GS))
		{
			#If one of the two is missing, skip
			if(is.na(PCRdat[vir,samp])|is.na(ELISAdat[vir,samp])) next
			#If RT is positive we write down and go to next
			if(PCRdat[vir,samp]>=RTthr & ELISAdat[vir,samp]=="+") 
			{
				GS[vir,samp]<-1
				next
			} 
			if(PCRdat[vir,samp]<RTthr & ELISAdat[vir,samp]=="-") 
			{
				GS[vir,samp]<-0
				next
			} 
		}
	}
}
write.table(GS,paste0(outdir,gold_standard,".txt"),sep="\t",quote=F,col.names=NA)
#Read results of RNA alignment 

#Change names of Viruses: only viruses with short name corresponding to those present in RT or ELISA files are included in analysis.
RNAdat<-fread(rnaseqfile,data.table=F)
tochange<-!names(RNAdat)%in%c("Name","Total")
names(RNAdat)[tochange]<-unlist(lapply(strsplit(names(RNAdat[tochange]),"_"),"[",3))
names(RNAdat)[tochange]<-unlist(lapply(lapply(strsplit(names(RNAdat[tochange]),"-"),"[",-1),paste,collapse="-"))
#Get short names
RNAdat$short<-NA
RNAdat$short[RNAdat$Name=="Arabis mosaic virus"]<-"ArMV"
RNAdat$short[RNAdat$Name=="Grapevine fleck virus"]<-"GFKV"
RNAdat$short[RNAdat$Name=="Grapevine virus A"]<-"GVA"
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
RNAdat$short[RNAdat$Name=="Grapevine fanleaf virus"]<-"GFLV"
RNAdat$short[RNAdat$Name=="Grapevine yellow speckle viroid 1"]<-"GySVd-1"
RNAdat$short[RNAdat$Name=="Hop stunt viroid"]<-"HSVd"
RNAdat$short[RNAdat$Name=="Grapevine Pinot gris virus"]<-"GPGV"
RNAdat$short[RNAdat$Name=="Grapevine rupestris stem pitting-associated virus"]<-"GRSPaV"
RNAdat$Name<-NULL
#Only check viruses tested with PCR
RNAdat<-RNAdat[!is.na(RNAdat$short),]
row.names(RNAdat)<-RNAdat$short
RNAdat$short<-RNAdat$Total<-NULL
RNAdat<-RNAdat[sort(row.names(RNAdat)),]
#Only select samples independent from the first round
RNAdat<-RNAdat[,valgeno]
RNAval<-RNAdat
#This is horrible, I know.
#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(RNAval))
{
	for(bbb in colnames(RNAval))
	{
	if(is.na(GS[aaa,bbb])) 
	{
		RNAval[aaa,bbb]<-NA
		next
	}
	if(GS[aaa,bbb]>0 & RNAdat[aaa,bbb]<rnathr) RNAval[aaa,bbb]<-"FN"
	if(GS[aaa,bbb]>0 & RNAdat[aaa,bbb]>=rnathr) RNAval[aaa,bbb]<-"TP"
	if(GS[aaa,bbb]<1 & RNAdat[aaa,bbb]<rnathr) RNAval[aaa,bbb]<-"TN"
	if(GS[aaa,bbb]<1 & RNAdat[aaa,bbb]>=rnathr) RNAval[aaa,bbb]<-"FP"
	}
}
write.table(RNAdat,paste0(outdir,"Tab_RNAa_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)
write.table(RNAval,paste0(outdir,"Val_RNAa_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)


#Read results of RNA assignment using bracken
Bdat<-fread(brackenfile,data.table=F)
tochange<-!names(Bdat)%in%c("name","taxonomy_id","Total")
names(Bdat)[tochange]<-unlist(lapply(strsplit(names(Bdat[tochange]),"_"),"[",4))
names(Bdat)[tochange]<-unlist(lapply(lapply(strsplit(names(Bdat[tochange]),"-"),"[",-1),paste,collapse="-"))
setnames(Bdat,"name","Name")
Bdat$short<-NA
Bdat$short[Bdat$Name=="Arabis mosaic virus"]<-"ArMV"
Bdat$short[Bdat$Name=="Grapevine fleck virus"]<-"GFKV"
Bdat$short[Bdat$Name=="Grapevine virus A"]<-"GVA"
Bdat$short[Bdat$Name=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
Bdat$short[Bdat$Name=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
Bdat$short[Bdat$Name=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
Bdat$short[Bdat$Name=="Grapevine fanleaf virus"]<-"GFLV"
Bdat$short[Bdat$Name=="Grapevine yellow speckle viroid 1"]<-"GySVd-1"
Bdat$short[Bdat$Name=="Hop stunt viroid"]<-"HSVd"
Bdat$short[Bdat$Name=="Grapevine Pinot gris virus"]<-"GPGV"
Bdat$short[Bdat$Name=="Grapevine rupestris stem pitting-associated virus"]<-"GRSPaV"
Bdat$Name<-NULL

#Only check viruses tested with PCR
Bdat<-Bdat[!is.na(Bdat$short),]
row.names(Bdat)<-Bdat$short
Bdat$short<-Bdat$Total<-Bdat$taxonomy_id<-NULL
Bdat<-Bdat[sort(row.names(Bdat)),]
#Only select samples independent from the first round
Bdat<-Bdat[,valgeno]
Bval<-Bdat



#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(Bval))
{
	for(bbb in colnames(Bval))
	{
	if(is.na(GS[aaa,bbb])) 
	{
		Bval[aaa,bbb]<-NA
		next
	}
	if(GS[aaa,bbb]>0 & Bdat[aaa,bbb]<brackthr) Bval[aaa,bbb]<-"FN"
	if(GS[aaa,bbb]>0 & Bdat[aaa,bbb]>=brackthr) Bval[aaa,bbb]<-"TP"
	if(GS[aaa,bbb]<1 & Bdat[aaa,bbb]<brackthr) Bval[aaa,bbb]<-"TN"
	if(GS[aaa,bbb]<1 & Bdat[aaa,bbb]>=brackthr) Bval[aaa,bbb]<-"FP"
	}
}

write.table(Bdat,paste0(outdir,"Tab_RNAb_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)
write.table(Bval,paste0(outdir,"Val_RNAb_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)


#Read VD output
smallcontig<-fread(smallrnacontigfile,data.table=F)
#Change sample names to the paper notation
names(smallcontig)<-gsub("55311ID1718COU121S201","CO-U1-2-1",names(smallcontig))
names(smallcontig)<-gsub("55312ID1718DOU211S211","DO-U2-1-1",names(smallcontig))
names(smallcontig)<-gsub("55313ID1718MB21S221","MB-2-1",names(smallcontig))
names(smallcontig)<-gsub("55315ID1718MN11S231","MN-1-1",names(smallcontig))
names(smallcontig)<-gsub("55316ID1718MN21S241","MN-2-1",names(smallcontig))
names(smallcontig)<-gsub("55317ID1718RI11S251","RI-1-1",names(smallcontig))
names(smallcontig)<-gsub("55318ID1718RI61S261","RI-6-1",names(smallcontig))
smallcontig$short<-NA
smallcontig$short[smallcontig$Species=="Arabis mosaic virus"]<-"ArMV"
smallcontig$short[smallcontig$Species=="Grapevine fleck virus"]<-"GFKV"
smallcontig$short[smallcontig$Species=="Grapevine virus A"]<-"GVA"
smallcontig$short[smallcontig$Species=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
smallcontig$short[smallcontig$Species=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
smallcontig$short[smallcontig$Species=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
smallcontig$short[smallcontig$Species=="Grapevine fanleaf virus"]<-"GFLV"
smallcontig$short[smallcontig$Species=="Grapevine yellow speckle viroid 1"]<-"GySVd-1"
smallcontig$short[smallcontig$Species=="Hop stunt viroid"]<-"HSVd"
smallcontig$short[smallcontig$Species=="Grapevine Pinot gris virus"]<-"GPGV"
smallcontig$short[smallcontig$Species=="Grapevine rupestris stem pitting-associated virus"]<-"GRSPaV"
smallcontig$Species<-NULL
smallcontig<-smallcontig[!is.na(smallcontig$short),]
row.names(smallcontig)<-smallcontig$short
smallcontig<-smallcontig[,grep("Percent",names(smallcontig))]
names(smallcontig)<-gsub("Percent_","",names(smallcontig))
#Only select samples independent from the first round
smallcontig<-smallcontig[,names(smallcontig)%in%valgeno]
SCval<-smallcontig
#Record TRUE and false positives in small RNA virus detect
for(aaa in row.names(smallcontig))
{
	for(bbb in colnames(smallcontig))
	{
		if(is.na(GS[aaa,bbb])) 
		{
			SCval[aaa,bbb]<-NA
			next
		}
	if(GS[aaa,bbb]>0 & smallcontig[aaa,bbb]<contigthr) SCval[aaa,bbb]<-"FN"
	if(GS[aaa,bbb]>0 & smallcontig[aaa,bbb]>=contigthr) SCval[aaa,bbb]<-"TP"
	if(GS[aaa,bbb]<1 & smallcontig[aaa,bbb]<contigthr) SCval[aaa,bbb]<-"TN"
	if(GS[aaa,bbb]<1 & smallcontig[aaa,bbb]>=contigthr) SCval[aaa,bbb]<-"FP"
	}
}

write.table(smallcontig,paste0(outdir,"Tab_smallc_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)
write.table(SCval,paste0(outdir,"Val_smallc_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)

#Read smallalignment output
smallal<-fread(smallrnafile,data.table=F)
smallal$Total<-NULL
#Use standard names for viruses and only retain the validated ones
smallal$short<-NA
smallal$short[smallal$Name=="Arabis mosaic virus"]<-"ArMV"
smallal$short[smallal$Name=="Grapevine fleck virus"]<-"GFKV"
smallal$short[smallal$Name=="Grapevine virus A"]<-"GVA"
smallal$short[smallal$Name=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
smallal$short[smallal$Name=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
smallal$short[smallal$Name=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
smallal$short[smallal$Name=="Grapevine fanleaf virus"]<-"GFLV"
smallal$short[smallal$Name=="Grapevine yellow speckle viroid 1"]<-"GySVd-1"
smallal$short[smallal$Name=="Hop stunt viroid"]<-"HSVd"
smallal$short[smallal$Name=="Grapevine Pinot gris virus"]<-"GPGV"
smallal$short[smallal$Name=="Grapevine rupestris stem pitting-associated virus"]<-"GRSPaV"
smallal$Name<-NULL
smallal<-smallal[!is.na(smallal$short),]
row.names(smallal)<-smallal$short
smallal$short<-NULL
#Change sample names
names(smallal)<-unlist(lapply(strsplit(names(smallal),"_"),"[",3))
#Only select samples independent from the first round
smallal<-smallal[,names(smallal)%in%valgeno]
ALval<-smallal

#Record TRUE and false positives in small RNA alignment
for(aaa in row.names(smallal))
{
	for(bbb in colnames(smallal))
	{
		if(is.na(GS[aaa,bbb])) 
		{
			ALval[aaa,bbb]<-NA
			next
		}
	if(GS[aaa,bbb]>0 & smallal[aaa,bbb]<salignthr) ALval[aaa,bbb]<-"FN"
	if(GS[aaa,bbb]>0 & smallal[aaa,bbb]>=salignthr) ALval[aaa,bbb]<-"TP"
	if(GS[aaa,bbb]<1 & smallal[aaa,bbb]<salignthr) ALval[aaa,bbb]<-"TN"
	if(GS[aaa,bbb]<1 & smallal[aaa,bbb]>=salignthr) ALval[aaa,bbb]<-"FP"
	}
}
write.table(smallal,paste0(outdir,"Tab_smalla_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)
write.table(ALval,paste0(outdir,"Val_smalla_",gold_standard,".txt"),quote=F,sep="\t",col.names=NA)

}
RTPCR_RNAseq(brackenfile=brackenfile,pcrfile=pcrfile,elisafile=elisafile,rnaseqfile=rnaseqfile,newgeno=newgeno,
smallrnafile=smallrnafile,smallrnacontigfile=smallrnacontigfile,outdir=outdir,gold_standard=gold_standard)


#https://blog.revolutionanalytics.com/2016/08/roc-curves-in-two-lines-of-code.html
