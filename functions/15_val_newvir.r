# Run with --help flag for help.
# Modified 26/01/2021 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-A", "--alignrnatest"), type="character", default="tables_val/val_RNAseq_VirDet.txt",
              help="RNA alignment file of test set [default= %default]", metavar="character"),
  make_option(c("-B", "--brackentest"), type="character", default="tables_val/val_RNA_bracken_VD_206.txt",
              help="Bracken summary file", metavar="character"),
  make_option(c("-C", "--smallrnacontigtest"), type="character", default="tables_val/smallRNA_contigs_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-S", "--smallrnatest"), type="character", default="tables_val/smallRNAseq_align.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-P", "--pcrtest"), type="character", default="Docs_val/Val_RT_Actin.csv",
              help="RT-PCR file", metavar="character"),
  make_option(c("-N", "--namesfile"), type="character", default="Docs_val/Sample_names.txt",
              help="Output table", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="tables_val/Val_new_virus.txt",
              help="Output table", metavar="character"),
  make_option(c("-F", "--remove.failed"), type="logical", default=TRUE,
              help="Remove the three problematic smallRNA files?", metavar="character"),
  make_option(c("-G", "--outdir"), type="character", default="tables_val/",
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

if (is.null(opt$brackentest)) {
  stop("WARNING: No brackentest specified with '-B' flag.")
} else {  cat ("Bracken summary file is ", opt$brackentest, "\n")
  brackentest <- opt$brackentest  
  }

if (is.null(opt$smallrnacontigtest)) {
  stop("WARNING: No smallrnacontigtest specified with '-C' flag.")
} else {  cat ("smallrnacontigtest file is ", opt$smallrnacontigtest, "\n")
  smallrnacontigtest <- opt$smallrnacontigtest  
  }

if (is.null(opt$smallrnatest)) {
  stop("WARNING: No smallrnatest specified with '-S' flag.")
} else {  cat ("smallrnatest file is ", opt$smallrnatest, "\n")
  smallrnatest <- opt$smallrnatest  
  }

if (is.null(opt$pcrtest)) {
  stop("WARNING: No pcrtest specified with '-P' flag.")
} else {  cat ("PCR test file is ", opt$pcrtest, "\n")
  pcrtest <- opt$pcrtest  
  }

  if (is.null(opt$outdir)) {
  stop("WARNING: No outdir specified with '-G' flag.")
} else {  cat ("outdir is ", opt$outdir, "\n")
  outdir <- opt$outdir  
  }

  if (is.null(opt$namesfile)) {
  stop("WARNING: No namesfile specified with '-N' flag.")
} else {  cat ("namesfile is ", opt$namesfile, "\n")
  namesfile <- opt$namesfile  
  }

  if (is.null(opt$remove.failed)) {
  stop("WARNING: No remove.failed specified with '-F' flag.")
} else {  cat ("remove.failed is ", opt$remove.failed, "\n")
  remove.failed <- opt$remove.failed  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


val_newvir<-function(alignrnatest,brackentest,smallrnacontigtest,smallrnatest,
						pcrtest,namesfile,
						outfile,outdir,remove.failed,
						RTthr=0.001,rnathr=19.28,brackthr=386,contigthr=0.97,salignthr=1572)
{
library("data.table")

#Another hardcoding... I really hope nobody ever sees this
longvir<-c("Grapevine yellow speckle viroid 1","Hop stunt viroid")
shortvir<-c("GySVd-1","HSVd")
#I select only samples analysed in the Phase 2
#These are the only ones for which I have the new viruses 
mynames<-fread(namesfile,data.table=F)
reps<-mynames[mynames$Phase>1,]
PCRdat<-fread(pcrtest,data.table=F)
names(PCRdat)<-gsub("\\\xa0","",names(PCRdat))  #Solve a stupid formatting problem
PCRdat$ID<-PCRdat$var<-NULL
rownames(PCRdat)<-PCRdat[,1]
PCRdat[,1]<-NULL
PCRdat<-t(PCRdat)
PCRdat<-PCRdat[sort(row.names(PCRdat)),]
#PCRdat<-round(PCRdat,3)
PCRdat<-PCRdat[shortvir,]

########################################
#
#RNA alignment (RNA a)
#
########################################

#Read RNA results for training set and test set
artest<-fread(alignrnatest,data.table=F)
artest$Total<-NULL
names(artest)<-gsub("_FPKM","",names(artest))
#Only keee the two viroids, change names
artest<-artest[artest$Name%in%longvir,]
#We want "GySVd-1" to be the first
artest<-artest[order(artest$Name),]
setnames(artest,reps$Samplename,reps$RT_Name)
artest$Name[artest$Name%in%longvir]<-shortvir[artest$Name%in%longvir]
row.names(artest)<-artest$Name
artest$Name<-NULL

RNAval<-artest
#Only keep replicated samples
#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(artest))
{
	for(bbb in colnames(artest))
	{
		if(PCRdat[aaa,bbb]>=RTthr & artest[aaa,bbb]<rnathr) RNAval[aaa,bbb]<-"FN"
		if(PCRdat[aaa,bbb]>=RTthr & artest[aaa,bbb]>=rnathr) RNAval[aaa,bbb]<-"TP"
		if(PCRdat[aaa,bbb]<RTthr & artest[aaa,bbb]<rnathr) RNAval[aaa,bbb]<-"TN"
		if(PCRdat[aaa,bbb]<RTthr & artest[aaa,bbb]>=rnathr) RNAval[aaa,bbb]<-"FP"
	}
}


########################################
#
#RNA bracken (RNA b)
#
########################################

brtest<-fread(brackentest,data.table=F)
brtest$Total<-NULL
names(brtest)<-gsub("RPM_","",names(brtest))
#Only keee the two viroids, change names
brtest<-brtest[brtest$name%in%longvir,]
setnames(brtest,reps$Bracken_name,reps$RT_Name)
brtest<-brtest[order(brtest$name),]
brtest$name[brtest$name%in%longvir]<-shortvir[brtest$name==longvir]
row.names(brtest)<-brtest$name
brtest$name<-brtest$taxonomy_id<-NULL
Bval<-brtest
#Only keep replicated samples
#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(brtest))
{
	for(bbb in colnames(brtest))
	{
		if(PCRdat[aaa,bbb]>=RTthr & brtest[aaa,bbb]<brackthr) Bval[aaa,bbb]<-"FN"
		if(PCRdat[aaa,bbb]>=RTthr & brtest[aaa,bbb]>=brackthr) Bval[aaa,bbb]<-"TP"
		if(PCRdat[aaa,bbb]<RTthr & brtest[aaa,bbb]<brackthr) Bval[aaa,bbb]<-"TN"
		if(PCRdat[aaa,bbb]<RTthr & brtest[aaa,bbb]>=brackthr) Bval[aaa,bbb]<-"FP"
	}
}



########################################
#
#small RNA align (small RNA a)
#
########################################

srtest<-fread(smallrnatest,data.table=F)
srtest$Total<-NULL
names(srtest)<-gsub("_FPKM","",names(srtest))
#Only keee the two viroids, change names
srtest<-srtest[srtest$Name%in%longvir,]
setnames(srtest,reps$smallalign_name,reps$RT_Name)
srtest$Name[srtest$Name%in%longvir]<-shortvir[srtest$Name%in%longvir]
row.names(srtest)<-srtest$Name
srtest$Name<-NULL
SAval<-srtest
#Only keep replicated samples
#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(srtest))
{
	for(bbb in colnames(srtest))
	{
		if(PCRdat[aaa,bbb]>=RTthr & srtest[aaa,bbb]<salignthr) SAval[aaa,bbb]<-"FN"
		if(PCRdat[aaa,bbb]>=RTthr & srtest[aaa,bbb]>=salignthr) SAval[aaa,bbb]<-"TP"
		if(PCRdat[aaa,bbb]<RTthr & srtest[aaa,bbb]<salignthr) SAval[aaa,bbb]<-"TN"
		if(PCRdat[aaa,bbb]<RTthr & srtest[aaa,bbb]>=salignthr) SAval[aaa,bbb]<-"FP"
	}
}

########################################
#
#small RNA contig (small RNA c)
#
########################################


sctest<-fread(smallrnacontigtest,data.table=F)
sctest$Total<-NULL
names(sctest)<-gsub("Percent_","",names(sctest))
sctest<-sctest[,grep("Number",names(sctest),invert=T)]
#Only keee the two viroids, change names
sctest<-sctest[sctest$Species%in%longvir,]
#Only select successful small RNA libraries
reps<-reps[!is.na(reps$smallcontig_name),]
setnames(sctest,reps$smallcontig_name,reps$RT_Name)
sctest$Species[sctest$Species%in%longvir]<-shortvir[sctest$Species%in%longvir]
row.names(sctest)<-sctest$Species
sctest$Species<-sctest$Taxid<-NULL
SCval<-sctest
#Only keep replicated samples
#Record TRUE and false positives in RNAseq virus detect
for(aaa in row.names(sctest))
{
	for(bbb in colnames(sctest))
	{
		if(PCRdat[aaa,bbb]>=RTthr & sctest[aaa,bbb]<contigthr) SCval[aaa,bbb]<-"FN"
		if(PCRdat[aaa,bbb]>=RTthr & sctest[aaa,bbb]>=contigthr) SCval[aaa,bbb]<-"TP"
		if(PCRdat[aaa,bbb]<RTthr & sctest[aaa,bbb]<contigthr) SCval[aaa,bbb]<-"TN"
		if(PCRdat[aaa,bbb]<RTthr & sctest[aaa,bbb]>=contigthr) SCval[aaa,bbb]<-"FP"
	}
}


#Remove the three low quality smallRNA samples
if(remove.failed) 
{
	cat("Removing failed smallRNA samples\n")
	SAval<-SAval[!names(SAval)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
	SCval<-SCval[!names(SCval)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
}





stats<-matrix(NA,nrow=4,ncol=12,dimnames=list(c("RNAa","RNAb","smalla","smallc"),c("Number","Sens","Spec","NPV","PPV","Accuracy","F1","MCC","TP","TN","FP","FN")))
for(aaa in 1:nrow(stats))
{
if(aaa==1) mytab<-RNAval
if(aaa==2) mytab<-Bval
if(aaa==3) mytab<-SAval
if(aaa==4) mytab<-SCval
mytab[is.na(mytab)]<-0
TP<-sum(mytab=="TP")
TN<-sum(mytab=="TN")
FN<-sum(mytab=="FN")
FP<-sum(mytab=="FP")
stats[aaa,"TP"]<-TP
stats[aaa,"TN"]<-TN
stats[aaa,"FP"]<-FP
stats[aaa,"FN"]<-FN
stats[aaa,"Number"]<-sum(TP+TN+FN+FP)
stats[aaa,"Sens"]<-round(TP/(TP+FN),3)
stats[aaa,"Spec"]<-round(TN/(TN+FP),3)
stats[aaa,"NPV"]<-round(TN/(TN+FN),3)
stats[aaa,"PPV"]<-round(TP/(TP+FP),3)
stats[aaa,"Accuracy"]<-round((TP+TN)/sum(TP+TN+FN+FP),3)
stats[aaa,"F1"]<-round((2*TP)/(2*TP+FP+FN),3)
stats[aaa,"MCC"]<-round((TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)),3)
}

nPCR<-data.frame(t(PCRdat),check.names=F)
names(nPCR)<-paste0("RT_",names(nPCR))
nartest<-data.frame(t(artest),check.names=F)
names(nartest)<-paste0("RNAa_",names(nartest))
nbrtest<-data.frame(t(brtest),check.names=F)
names(nbrtest)<-paste0("RNAb_",names(nbrtest))
nsrtest<-data.frame(t(srtest),check.names=F)
names(nsrtest)<-paste0("smalla_",names(nsrtest))
nsctest<-data.frame(t(sctest),check.names=F)
names(nsctest)<-paste0("smallc_",names(nsctest))
final_a<-merge(nPCR,nartest,by="row.names",all=T)
row.names(final_a)<-final_a$Row.names
final_a$Row.names<-NULL
final_b<-merge(final_a,nbrtest,by="row.names",all=T)
row.names(final_b)<-final_b$Row.names
final_b$Row.names<-NULL
final_c<-merge(final_b,nsrtest,by="row.names",all=T)
row.names(final_c)<-final_c$Row.names
final_c$Row.names<-NULL
final_d<-merge(final_c,nsctest,by="row.names",all=T)
row.names(final_d)<-final_d$Row.names
final_d$Row.names<-NULL

write.table(stats,paste(outdir,"stats_newvir.txt",sep="/"),quote=F,sep="\t",col.names=NA)
write.table(final_d,outfile,quote=F,sep="\t",col.names=NA)


}
val_newvir(alignrnatest=alignrnatest,brackentest=brackentest,
smallrnacontigtest=smallrnacontigtest,smallrnatest=smallrnatest,
pcrtest=pcrtest,namesfile=namesfile,remove.failed=remove.failed,
outfile=outfile,outdir=outdir)

