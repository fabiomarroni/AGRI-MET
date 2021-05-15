# Run with --help flag for help.
# Modified 26/01/2021 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--alignrnatest"), type="character", default="tables_val/val_RNAseq_VirDet.txt",
              help="RNA alignment file of test set [default= %default]", metavar="character"),
  make_option(c("-a", "--alignrnatrain"), type="character", default="tables/RNAseq_VirDet.txt",
              help="RNA alignment results of training set [default= %default]", metavar="character"),
  make_option(c("-B", "--brackentest"), type="character", default="tables_val/val_RNA_bracken_VD_206.txt",
              help="Bracken summary file", metavar="character"),
  make_option(c("-b", "--brackentrain"), type="character", default="tables/RNA_bracken_VD_206.txt",
              help="Bracken summary file", metavar="character"),
  make_option(c("-C", "--smallrnacontigtest"), type="character", default="tables_val/smallRNA_contigs_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-c", "--smallrnacontigtrain"), type="character", default="tables/smallRNA_contigs_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-S", "--smallrnatest"), type="character", default="tables_val/smallRNAseq_align.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-s", "--smallrnatrain"), type="character", default="tables/smallRNAseq_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-P", "--pcrtest"), type="character", default="Docs_val/Val_RT_Actin.csv",
              help="RT-PCR file", metavar="character"),
  make_option(c("-p", "--pcrtrain"), type="character", default="tables/RT-PCR.txt",
              help="RT-PCR file", metavar="character"),
  make_option(c("-E", "--elisatest"), type="character", default="Docs_val/Val_ELISA.csv",
              help="ELISA file", metavar="character"),
  make_option(c("-e", "--elisatrain"), type="character", default="tables/Elisa_res.txt",
              help="ELISA file", metavar="character"),
  make_option(c("-N", "--namesfile"), type="character", default="Docs_val/Sample_names.txt",
              help="Output table", metavar="character"),
  make_option(c("-V", "--virus"), type="character", default="core",
              help="Virus set to analyze. Possible values: 'core','novel','all'. See code for details [default= %default]", metavar="character"),
  make_option(c("-R", "--remove_changed"), type="logical", default=TRUE,
              help="Remove form analysis comparisons for which both ELISA and RT-PCR revealed changes between first and second analysis", metavar="character"),
  make_option(c("-O", "--outdir"), type="character", default="tables_techval/",
              help="Fasta index file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$alignrnatest)) {
  stop("WARNING: No alignrnatest specified with '-A' flag.")
} else {  cat ("RNA seq align file is ", opt$alignrnatest, "\n")
  alignrnatest <- opt$alignrnatest  
  }

if (is.null(opt$alignrnatrain)) {
  stop("WARNING: No training set align RNA fle specified with '-a' flag.")
} else {  cat ("RNA seq training set file is ", opt$alignrnatrain, "\n")
  alignrnatrain <- opt$alignrnatrain  
  }

if (is.null(opt$brackentest)) {
  stop("WARNING: No brackentest specified with '-B' flag.")
} else {  cat ("Bracken summary file is ", opt$brackentest, "\n")
  brackentest <- opt$brackentest  
  }

if (is.null(opt$brackentrain)) {
  stop("WARNING: No brackentrain specified with '-b' flag.")
} else {  cat ("Bracken training set summary file is ", opt$brackentrain, "\n")
  brackentrain <- opt$brackentrain  
  }

if (is.null(opt$smallrnacontigtest)) {
  stop("WARNING: No smallrnacontigtest specified with '-C' flag.")
} else {  cat ("smallrnacontigtest file is ", opt$smallrnacontigtest, "\n")
  smallrnacontigtest <- opt$smallrnacontigtest  
  }

if (is.null(opt$smallrnacontigtrain)) {
  stop("WARNING: No smallrnacontigtrain specified with '-c' flag.")
} else {  cat ("smallrnacontigtrain file is ", opt$smallrnacontigtrain, "\n")
  smallrnacontigtrain <- opt$smallrnacontigtrain  
  }

if (is.null(opt$smallrnatest)) {
  stop("WARNING: No smallrnatest specified with '-S' flag.")
} else {  cat ("smallrnatest file is ", opt$smallrnatest, "\n")
  smallrnatest <- opt$smallrnatest  
  }

if (is.null(opt$smallrnatrain)) {
  stop("WARNING: No smallrnatrain specified with '-s' flag.")
} else {  cat ("smallrnatrain file is ", opt$smallrnatrain, "\n")
  smallrnatrain <- opt$smallrnatrain  
  }

if (is.null(opt$pcrtest)) {
  stop("WARNING: No pcrtest specified with '-P' flag.")
} else {  cat ("PCR test file is ", opt$pcrtest, "\n")
  pcrtest <- opt$pcrtest  
  }

if (is.null(opt$pcrtrain)) {
  stop("WARNING: No pcrtrain specified with '-p' flag.")
} else {  cat ("PCR train file is ", opt$pcrtrain, "\n")
  pcrtrain <- opt$pcrtrain  
  }

  if (is.null(opt$elisatest)) {
  stop("WARNING: No elisatest specified with '-E' flag.")
} else {  cat ("ELISA file is ", opt$elisatest, "\n")
  elisatest <- opt$elisatest  
  }

  if (is.null(opt$elisatrain)) {
  stop("WARNING: No elisatrain specified with '-e' flag.")
} else {  cat ("ELISA train set file is ", opt$elisatrain, "\n")
  elisatrain <- opt$elisatrain  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No outdir specified with '-O' flag.")
} else {  cat ("outdir is ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

  if (is.null(opt$namesfile)) {
  stop("WARNING: No namesfile specified with '-N' flag.")
} else {  cat ("namesfile is ", opt$namesfile, "\n")
  namesfile <- opt$namesfile  
  }

  if (is.null(opt$virus)) {
  stop("WARNING: No virus specified with '-V' flag.")
} else {  cat ("virus is ", opt$virus, "\n")
  virus <- opt$virus  
  }

  if (is.null(opt$remove_changed)) {
  stop("WARNING: No remove_changed specified with '-R' flag.")
} else {  cat ("remove_changed is ", opt$remove_changed, "\n")
  remove_changed <- opt$remove_changed  
  }


#Done: compute repeatability for RNAa, RT-PCR and ELISA
#Todo: 1) compute repeatability for RNAb, smallRNAa, and smallRNAc
#      2) see if we want to perform the comparison only on the 7 viruses (I am now using all the viruses that tested positive in at least one sample,
#		  but this may overestimate TNR) 
compare_duplicate<-function(alignrnatest,alignrnatrain,brackentest,brackentrain,smallrnacontigtest,smallrnacontigtrain,smallrnatest,smallrnatrain,
						pcrtest,pcrtrain,elisatest,elisatrain,namesfile,
						outdir,virus,
						remove_changed,RTthr=0.001,rnathr=19.28,brackthr=386,contigthr=0.97,salignthr=1572,col.list=c("black","grey"))
{
library("data.table")
library(scales)
library(RColorBrewer)
#I know I am going to regret, but I hardcode here the long names of viruses for which RT-PCR or ELISA are available for both rounds.
#Only these are used to check technical replicability if onlycore=TRUE 
longvir<-c("Arabis mosaic virus","Grapevine fanleaf virus","Grapevine fleck virus","Grapevine leafroll-associated virus 1",
			"Grapevine leafroll-associated virus 2","Grapevine leafroll-associated virus 3",
			"Grapevine Pinot gris virus","Grapevine rupestris stem pitting-associated virus","Grapevine virus A")
#Another stupid hardcoding. We may not need this one!
shortvir<-c("ArMV","GFLV","GFKV","GLRaV-1","GLRaV-2","GLRaV-3","GPGV","GRSPaV","GVA")
#Identify samples that were analysed both in training and test set 
mynames<-fread(namesfile,data.table=F)
dup<-mynames$Code[duplicated(mynames$Code)]
reps<-mynames[mynames$Code%in%dup,]
########################################
#
#RNA alignment (RNA a)
#
########################################

#Read RNA results for training set and test set
artest<-fread(alignrnatest,data.table=F)
artrain<-fread(alignrnatrain,data.table=F)
artest$Total<-artrain$Total<-NULL
names(artest)<-gsub("_FPKM","",names(artest))
names(artrain)<-gsub("_FPKM","",names(artrain))
#Merge train and test set and only keep virus names and replicated samples
dupdat<-merge(artest,artrain,by="Name",all=T)
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$Samplename])]
if(virus=="core") dupdat<-dupdat[dupdat$Name%in%longvir,]
if(virus=="new") dupdat<-dupdat[!dupdat$Name%in%longvir,]

#techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA,vTN=NA,vTP=NA,vFN=NA,vFP=NA)
techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
#Only keep replicated samples
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$Samplename[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
ordup<-sdup
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<rnathr]<-0
sdup[,compnames][sdup[,compnames]>rnathr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
# techval[aaa,"vTN"]<-mean(unlist(ordup[,c(compnames[1],compnames[2])][sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0",]))
# techval[aaa,"vFP"]<-mean(unlist(ordup[,c(compnames[1],compnames[2])][sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1",]))
# techval[aaa,"vFN"]<-mean(unlist(ordup[,c(compnames[1],compnames[2])][sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0",]))
# techval[aaa,"vTP"]<-mean(unlist(ordup[,c(compnames[1],compnames[2])][sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1",]))
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$Samplename%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/RNAalign_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"/RNAalign_techval_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
rnaa<-techval

########################################
#
#RNA bracken (RNA b)
#
########################################

#Read RNA results for training set and test set
brtest<-fread(brackentest,data.table=F)
brtrain<-fread(brackentrain,data.table=F)
brtest$Total<-brtrain$Total<-brtest$taxonomy_id<-brtrain$taxonomy_id<-NULL
names(brtest)<-gsub("RPM_","",names(brtest))
names(brtrain)<-gsub("RPM_","",names(brtrain))
names(brtest)[names(brtest)=="name"]<-"Name"
names(brtrain)[names(brtrain)=="name"]<-"Name"
#Merge train and test set and only keep virus names and replicated samples
dupdat<-merge(brtest,brtrain,by="Name",all=T)
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$Bracken_name])]
if(virus=="core") dupdat<-dupdat[dupdat$Name%in%longvir,]
if(virus=="new") dupdat<-dupdat[!dupdat$Name%in%longvir,]

techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
#Only keep replicated samples
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$Bracken_name[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<brackthr]<-0
sdup[,compnames][sdup[,compnames]>=brackthr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$Bracken_name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/RNAbracken_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"/RNAbracken_techval_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
rnab<-techval


########################################
#
#small RNA align (small RNA a)
#
########################################

#Read RNA results for training set and test set
srtest<-fread(smallrnatest,data.table=F)
srtrain<-fread(smallrnatrain,data.table=F)
names(srtest)<-gsub("_FPKM","",names(srtest))
names(srtrain)<-gsub("_FPKM","",names(srtrain))
srtest$Total<-srtrain$Total<-NULL
#Merge train and test set and only keep virus names and replicated samples
dupdat<-merge(srtest,srtrain,by="Name",all=T)
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$smallalign_name])]
if(virus=="core") dupdat<-dupdat[dupdat$Name%in%longvir,]
if(virus=="new") dupdat<-dupdat[!dupdat$Name%in%longvir,]
techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
#Only keep replicated samples
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$smallalign_name[reps$Code==dup[aaa]]
#Special check for incomplete smallRNA datasets
if(sum(compnames%in%names(dupdat))<2) next
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<salignthr]<-0
sdup[,compnames][sdup[,compnames]>=salignthr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$smallalign_name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/smallRNAalign_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"/smallRNAalign_techval_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
smalla<-techval
########################################
#
#small RNA contig (small RNA c)
#
########################################
#Read RNA results for training set and test set
sctest<-fread(smallrnacontigtest,data.table=F)
sctrain<-fread(smallrnacontigtrain,data.table=F)
names(sctest)<-gsub("Percent_","",names(sctest))
names(sctrain)<-gsub("Percent_","",names(sctrain))
sctest$Total<-sctrain$Total<-sctest$Taxid<-sctrain$Total<-NULL
#Merge train and test set and only keep virus names and replicated samples
names(sctest)[names(sctest)=="Virus"|names(sctest)=="Species"]<-"Name"
names(sctrain)[names(sctrain)=="Virus"|names(sctrain)=="Species"]<-"Name"
dupdat<-merge(sctest,sctrain,by="Name",all=T)
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$smallcontig_name])]
if(virus=="core") dupdat<-dupdat[dupdat$Name%in%longvir,]
if(virus=="new") dupdat<-dupdat[!dupdat$Name%in%longvir,]
techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
#Only keep replicated samples
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$smallcontig_name[reps$Code==dup[aaa]]
#Special check for incomplete smallRNA datasets
if(sum(compnames%in%names(dupdat))<2) next
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<contigthr]<-0
sdup[,compnames][sdup[,compnames]>=contigthr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$smallcontig_name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/smallRNAcontig_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"/smallRNAcontig_techval_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
smallc<-techval



##########################################
#
#REAL TIME
#
##########################################


#Read RT-PCR results for training set and test set and compare them 
rtest<-fread(pcrtest,data.table=F)
names(rtest)<-gsub("\\\xa0","",names(rtest))  #Solve a stupid formatting problem
rtest$ID<-rtest$var<-NULL
row.names(rtest)<-rtest$sample
rtest$sample<-NULL
#RT files are transposed (because I am stupid) so I ahve to transform them 
rtest<-data.frame(t(rtest),check.names=F)
rtrain<-fread(pcrtrain,data.table=F)
rtrain$"viral_RNA/COX_mRNA"<-gsub("_FPKM","",rtrain$"viral_RNA/COX_mRNA")
rtrain[is.na(rtrain)]<-0
row.names(rtrain)<-rtrain$"viral_RNA/COX_mRNA"
rtrain$"viral_RNA/COX_mRNA"<-NULL
#RT files are transposed (because I am stupid) so I ahve to transform them 
rtrain<-data.frame(t(rtrain),check.names=F)
dupdat<-merge(rtest,rtrain,by="row.names")
setnames(dupdat,"Row.names","Name")
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$RT_Name])]
techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$RT_Name[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<RTthr]<-0
sdup[,compnames][sdup[,compnames]>RTthr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$RT_Name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/RT-PCR_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")

write.table(techval,paste(outdir,"RT-PCR_techval.txt",sep="/"),quote=F,row.names=F,sep="\t")
rtpcr<-techval





#Read ELISA results for training set and test set and compare them 
#Read ELISA results for test set
etest<-fread(elisatest,data.table=F)
names(etest)<-gsub("\\\xa0","",names(etest))  #Solve a stupid formatting problem
etest$ID<-etest$var<-NULL
row.names(etest)<-etest$sample
etest$sample<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etest<-data.frame(t(etest),check.names=F,stringsAsFactors=F)
etest[etest=="+"]<-1
etest[etest=="-"]<-0
etest[etest=="NaN"]<-0
#Read ELISA result for training set
etrain<-fread(elisatrain,data.table=F)
#Solve stupid case issue
names(etrain)[names(etrain)=="GFkV"]<-"GFKV"
row.names(etrain)<-etrain$"Sample name"
etrain$"Sample name"<-etrain$"Code"<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etrain<-data.frame(t(etrain),check.names=F)
dupdat<-merge(etest,etrain,by="row.names")
setnames(dupdat,"Row.names","Name")
dupdat<-dupdat[,c("Name",names(dupdat)[names(dupdat)%in%reps$ELISA_name])]
techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$ELISA_name[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<RTthr]<-0
sdup[,compnames][sdup[,compnames]>RTthr]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$ELISA_name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/ELISA_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"ELISA_techval.txt",sep="/"),quote=F,row.names=F,sep="\t")
elisa<-techval



#Create the either gold standard merging ELISA and RT-PCR 
#Read ELISA results for test set
etest<-fread(elisatest,data.table=F)
names(etest)<-gsub("\\\xa0","",names(etest))  #Solve a stupid formatting problem
etest$ID<-etest$var<-NULL
row.names(etest)<-etest$sample
etest$sample<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etest<-data.frame(t(etest),check.names=F,stringsAsFactors=F)
etest[etest=="+"]<-1
etest[etest=="-"]<-0
etest[etest=="NaN"]<-0
#Read ELISA result for training set
etrain<-fread(elisatrain,data.table=F)
#Solve stupid case issue
names(etrain)[names(etrain)=="GFkV"]<-"GFKV"
row.names(etrain)<-etrain$"Sample name"
etrain$"Sample name"<-etrain$"Code"<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etrain<-data.frame(t(etrain),check.names=F)
etrain[is.na(etrain)]<-0
#Read RT-PCR data
rtest<-fread(pcrtest,data.table=F)
names(rtest)<-gsub("\\\xa0","",names(rtest))  #Solve a stupid formatting problem
rtest$ID<-rtest$var<-NULL
row.names(rtest)<-rtest$sample
rtest$sample<-NULL
#RT files are transposed (because I am stupid) so I have to transform them 
rtest<-data.frame(t(rtest),check.names=F)
rtrain<-fread(pcrtrain,data.table=F)
rtrain$"viral_RNA/COX_mRNA"<-gsub("_FPKM","",rtrain$"viral_RNA/COX_mRNA")
row.names(rtrain)<-rtrain$"viral_RNA/COX_mRNA"
rtrain$"viral_RNA/COX_mRNA"<-NULL
rtrain[is.na(rtrain)]<-0
rtest[rtest<RTthr]<-0
rtest[rtest>=RTthr]<-1
rtrain[rtrain<RTthr]<-0
rtrain[rtrain>=RTthr]<-1
#RT files are transposed (because I am stupid) so I have to transform them 
rtrain<-data.frame(t(rtrain),check.names=F)
#Create merged datasets for elisa and RT-PCR
edupdat<-merge(etest,etrain,by="row.names")
setnames(edupdat,"Row.names","Name")
edupdat<-edupdat[,c("Name",names(edupdat)[names(edupdat)%in%reps$ELISA_name])]
edupdat<-merge(edupdat,data.frame(Name=shortvir),all=T)
edupdat[is.na(edupdat)]<-0
rdupdat<-merge(rtest,rtrain,by="row.names")
setnames(rdupdat,"Row.names","Name")
rdupdat<-rdupdat[,c("Name",names(rdupdat)[names(rdupdat)%in%reps$RT_Name])]
#Column order and row order are the same in rdupdat and edupdat (although colnames are slightly different in the two)
rdupdat<-merge(rdupdat,data.frame(Name=shortvir),all=T)
rdupdat[is.na(rdupdat)]<-0
alldupdat<-rdupdat
for(aaa in 2:ncol(alldupdat))
{
alldupdat[,aaa]<-alldupdat[,aaa]+as.numeric(edupdat[,aaa])
}
dupdat<-alldupdat

techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$RT_Name[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
sdup[is.na(sdup)]<-0
sdup[,compnames][sdup[,compnames]<1]<-0
sdup[,compnames][sdup[,compnames]>0]<-1
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0")
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1")
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0")
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1")
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$RT_Name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/either_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"either_techval.txt",sep="/"),quote=F,row.names=F,sep="\t")
goldst<-techval


#Create the concordant gold standard merging ELISA and RT-PCR 
#Read ELISA results for test set
etest<-fread(elisatest,data.table=F)
names(etest)<-gsub("\\\xa0","",names(etest))  #Solve a stupid formatting problem
etest$ID<-etest$var<-NULL
row.names(etest)<-etest$sample
etest$sample<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etest<-data.frame(t(etest),check.names=F,stringsAsFactors=F)
etest[etest=="+"]<-1
etest[etest=="-"]<-0
etest[etest=="NaN"]<-NA
#Read ELISA result for training set
etrain<-fread(elisatrain,data.table=F)
#Solve stupid case issue
names(etrain)[names(etrain)=="GFkV"]<-"GFKV"
row.names(etrain)<-etrain$"Sample name"
etrain$"Sample name"<-etrain$"Code"<-NULL
#ELISA files are transposed (because I am stupid) so I ahve to transform them 
etrain<-data.frame(t(etrain),check.names=F)
#Read RT-PCR data
rtest<-fread(pcrtest,data.table=F)
names(rtest)<-gsub("\\\xa0","",names(rtest))  #Solve a stupid formatting problem
rtest$ID<-rtest$var<-NULL
row.names(rtest)<-rtest$sample
rtest$sample<-NULL
#RT files are transposed (because I am stupid) so I have to transform them 
rtest<-data.frame(t(rtest),check.names=F)
rtrain<-fread(pcrtrain,data.table=F)
rtrain$"viral_RNA/COX_mRNA"<-gsub("_FPKM","",rtrain$"viral_RNA/COX_mRNA")
row.names(rtrain)<-rtrain$"viral_RNA/COX_mRNA"
rtrain$"viral_RNA/COX_mRNA"<-NULL
#This is the only case in which NA is actually meaning zero
rtrain[is.na(rtrain)]<-0
rtest[rtest<RTthr]<-0
rtest[rtest>=RTthr]<-1
rtrain[rtrain<RTthr]<-0
rtrain[rtrain>=RTthr]<-1
#RT files are transposed (because I am stupid) so I have to transform them 
rtrain<-data.frame(t(rtrain),check.names=F)
#Create the concordant set for training and test set
sharedvir<-Reduce(intersect,list(row.names(etrain),row.names(rtrain),row.names(etest),row.names(rtest)))
ctrain<-matrix(NA,nrow=length(sharedvir),ncol=ncol(rtrain),dimnames=list(sharedvir,colnames(rtrain)))
ctest<-matrix(NA,nrow=length(sharedvir),ncol=ncol(rtest),dimnames=list(sharedvir,colnames(rtest)))
#Loop on training set
for(aaa in row.names(ctrain))
{
	for(bbb in colnames(ctrain))
	{
	#If ELISA was not tested we cnnot have agreement and we skip!
	if(is.na(etrain[aaa,bbb])) next
	if(etrain[aaa,bbb]==rtrain[aaa,bbb]) ctrain[aaa,bbb]<-etrain[aaa,bbb]
	}
}

#Loop on test set
#I checked that samples are in the same order. I need to remove the "sulla falcon" sentence and I was lazy
colnames(etest)<-colnames(rtest)
for(aaa in row.names(ctest))
{
	for(bbb in colnames(ctest))
	{
	#If ELISA was not tested we cnnot have agreement and we skip!
	if(is.na(etest[aaa,bbb])) next
	if(etest[aaa,bbb]==rtest[aaa,bbb]) ctest[aaa,bbb]<-etest[aaa,bbb]
	}
}

#Create merged datasets for elisa and RT-PCR
edupdat<-merge(ctest,ctrain,by="row.names")
setnames(edupdat,"Row.names","Name")
#This is not a mistake, I did it because I already changed the name of ELISA to RT-PCR
edupdat<-edupdat[,c("Name",names(edupdat)[names(edupdat)%in%reps$RT_Name])]
dupdat<-edupdat

techval<-data.frame(ID=dup,TN=NA,TP=NA,FN=NA,FP=NA)
for(aaa in 1:length(dup))
{
#Create a small df with the two replicates
compnames<-reps$RT_Name[reps$Code==dup[aaa]]
sdup<-dupdat[,c("Name",compnames)]
sdup[,2:ncol(sdup)]<-apply(sdup[,2:ncol(sdup)],2,as.numeric)
#sdup[is.na(sdup)]<-0
techval[aaa,"TN"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="0",na.rm=T)
techval[aaa,"FP"]<-sum(sdup[,compnames[1]]=="0"&sdup[,compnames[2]]=="1",na.rm=T)
techval[aaa,"FN"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="0",na.rm=T)
techval[aaa,"TP"]<-sum(sdup[,compnames[1]]=="1"&sdup[,compnames[2]]=="1",na.rm=T)
}
#Yes, it sucks!!!
for(aaa in 2:length(names(dupdat)))
{
names(dupdat)[aaa]<-reps$Short[reps$RT_Name%in%names(dupdat)[aaa]]
}
write.table(dupdat,paste(outdir,"/concordant_techval_data_",virus,".txt",sep=""),quote=F,row.names=F,sep="\t")
write.table(techval,paste(outdir,"concordant_techval.txt",sep="/"),quote=F,row.names=F,sep="\t")
concordant<-techval

if(remove_changed)
{
browser()
}

stats<-matrix(NA,nrow=8,ncol=12,dimnames=list(c("RNAa","RNAb","smalla","smallc","elisa","rtpcr","either","concordant"),c("Number","Sens","Spec","NPV","PPV","Accuracy","F1","MCC","TP","TN","FP","FN")))
for(aaa in 1:nrow(stats))
{
if(aaa==1) mytab<-rnaa
if(aaa==2) mytab<-rnab
if(aaa==3) mytab<-smalla
if(aaa==4) mytab<-smallc
if(aaa==5) mytab<-elisa
if(aaa==6) mytab<-rtpcr
if(aaa==7) mytab<-goldst
if(aaa==8) mytab<-concordant
mytab[is.na(mytab)]<-0
stats[aaa,"TP"]<-TP<-sum(mytab["TP"])
stats[aaa,"TN"]<-TN<-sum(mytab["TN"])
stats[aaa,"FN"]<-FN<-sum(mytab["FN"])
stats[aaa,"FP"]<-FP<-sum(mytab["FP"])
stats[aaa,"Number"]<-sum(TP+TN+FN+FP)
stats[aaa,"Sens"]<-round(TP/(TP+FN),3)
stats[aaa,"Spec"]<-round(TN/(TN+FP),3)
stats[aaa,"NPV"]<-round(TN/(TN+FN),3)
stats[aaa,"PPV"]<-round(TP/(TP+FP),3)
stats[aaa,"Accuracy"]<-round((TP+TN)/sum(TP+TN+FN+FP),3)
stats[aaa,"F1"]<-round((2*TP)/(2*TP+FP+FN),3)
stats[aaa,"MCC"]<-round((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),3)
}

write.table(stats,paste(outdir,"stats_techval.txt",sep="/"),quote=F,sep="\t",col.names=NA)

}
compare_duplicate(alignrnatest=alignrnatest,alignrnatrain=alignrnatrain,brackentest=brackentest,brackentrain=brackentrain,
smallrnacontigtest=smallrnacontigtest,smallrnacontigtrain=smallrnacontigtrain,smallrnatest=smallrnatest,smallrnatrain=smallrnatrain,
pcrtest=pcrtest,pcrtrain=pcrtrain,elisatest=elisatest,elisatrain=elisatrain,namesfile=namesfile,
outdir=outdir,virus=virus,remove_changed=remove_changed)


