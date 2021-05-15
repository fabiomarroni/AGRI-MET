# Run with --help flag for help.
# Modified 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-B", "--brackenfile"), type="character", default="tables_val/Val_RNAb_concordant.txt",
              help="Bracken summary file", metavar="character"),
  make_option(c("-R", "--rnaseqfile"), type="character", default="tables_val/Val_RNAa_concordant.txt",
              help="Fasta index file", metavar="character"),
  make_option(c("-S", "--smallrnafile"), type="character", default="tables_val/Val_smalla_concordant.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-C", "--smallrnacontigfile"), type="character", default="tables_val/Val_smallc_concordant.txt",
              help="Small RNA assembly file", metavar="character"),
  make_option(c("-F", "--remove.failed"), type="logical", default=TRUE,
              help="Remove the three problematic smallRNA files?", metavar="character"),
  make_option(c("-V", "--virus"), type="character", default="core",
              help="Virus set to analyze. Possible values: 'core','novel','all'. See code for details [default= %default]", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="tables_val/final_stats_core_concordant.txt",
              help="Output table", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$brackenfile)) {
  stop("WARNING: No brackenfile specified with '-B' flag.")
} else {  cat ("Bracken summary file is ", opt$brackenfile, "\n")
  brackenfile <- opt$brackenfile  
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

  if (is.null(opt$remove.failed)) {
  stop("WARNING: No remove.failed specified with '-F' flag.")
} else {  cat ("remove.failed is ", opt$remove.failed, "\n")
  remove.failed <- opt$remove.failed  
  }

  if (is.null(opt$virus)) {
  stop("WARNING: No virus specified with '-V' flag.")
} else {  cat ("virus is ", opt$virus, "\n")
  virus <- opt$virus  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


  
 
explore_val<-function(rnaseqfile,smallrnafile,smallrnacontigfile,brackenfile,outfile,
						RTthr=0.001,rnathr=19.28,brackthr=386,contigthr=0.97,salignthr=1572,remove.failed,virus)
{
library("data.table")
library(scales)
library(RColorBrewer)
#We assume that the files have a name starting with "Val" and they report the result of validation, i.e. TP, FP,TN,FN
#We also assume that files with the same path and name but with Tab instead of Val are present and report the results of each method. 


#Another stupid hardcoding. Needed if we want to perform the validation only using viruses that have laready been tested in Phase 1
shortvir<-c("ArMV","GFLV","GFKV","GLRaV-1","GLRaV-2","GLRaV-3","GPGV","GRSPaV","GVA")

alignVal<-fread(rnaseqfile,data.table=F,header=T)
brackenVal<-fread(brackenfile,data.table=F,header=T)
salignVal<-fread(smallrnafile,data.table=F,header=T)
scontigVal<-fread(smallrnacontigfile,data.table=F,header=T)
alignTab<-fread(gsub("/Val_","/Tab_",rnaseqfile),data.table=F,header=T)
brackenTab<-fread(gsub("/Val_","/Tab_",brackenfile),data.table=F,header=T)
salignTab<-fread(gsub("/Val_","/Tab_",smallrnafile),data.table=F,header=T)
scontigTab<-fread(gsub("/Val_","/Tab_",smallrnacontigfile),data.table=F,header=T)
row.names(alignVal)<-alignVal[,1]
row.names(brackenVal)<-brackenVal[,1]
row.names(salignVal)<-salignVal[,1]
row.names(scontigVal)<-scontigVal[,1]
alignVal[,1]<-brackenVal[,1]<-salignVal[,1]<-scontigVal[,1]<-NULL
row.names(alignTab)<-alignTab[,1]
row.names(brackenTab)<-brackenTab[,1]
row.names(salignTab)<-salignTab[,1]
row.names(scontigTab)<-scontigTab[,1]
alignTab[,1]<-brackenTab[,1]<-salignTab[,1]<-scontigTab[,1]<-NULL
if(virus=="core")
{
alignVal<-alignVal[row.names(alignVal)%in%shortvir,]
brackenVal<-brackenVal[row.names(brackenVal)%in%shortvir,]
salignVal<-salignVal[row.names(salignVal)%in%shortvir,]
scontigVal<-scontigVal[row.names(scontigVal)%in%shortvir,]
alignTab<-alignTab[row.names(alignTab)%in%shortvir,]
brackenTab<-brackenTab[row.names(brackenTab)%in%shortvir,]
salignTab<-salignTab[row.names(salignTab)%in%shortvir,]
scontigTab<-scontigTab[row.names(scontigTab)%in%shortvir,]
}

#Remove the three low quality smallRNA samples
if(remove.failed) 
{
	cat("Removing failed smallRNA samples\n")
	salignTab<-salignTab[!names(salignTab)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
	salignVal<-salignVal[!names(salignVal)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
	scontigTab<-scontigTab[!names(scontigTab)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
	scontigVal<-scontigVal[!names(scontigVal)%in%c("CO-U1-1-1","DO-U2-2-1","MB-5-1")]
}

aV<-table(unlist(alignVal))
bV<-table(unlist(brackenVal))
scV<-table(unlist(scontigVal))
saV<-table(unlist(salignVal))
cat("Valid data point for smallRNAc", sum(scV),"\n")
cat("Valid data point for smallRNAa", sum(saV),"\n")


stats<-matrix(NA,nrow=4,ncol=12,dimnames=list(c("RNAa","RNAb","smalla","smallc"),c("Number","Sens","Spec","NPV","PPV","Accuracy","F1","MCC","TP","TN","FP","FN")))
for(aaa in 1:nrow(stats))
{
if(aaa==1) mytab<-aV
if(aaa==2) mytab<-bV
if(aaa==3) mytab<-saV
if(aaa==4) mytab<-scV
if(is.na(mytab["FN"])) mytab["FN"]<-0
if(is.na(mytab["FP"])) mytab["FP"]<-0
if(is.na(mytab["TN"])) mytab["TN"]<-0
if(is.na(mytab["TP"])) mytab["TP"]<-0
stats[aaa,"Number"]<-sum(mytab)
stats[aaa,"Sens"]<-round(mytab["TP"]/(mytab["TP"]+mytab["FN"]),3)
stats[aaa,"Spec"]<-round(mytab["TN"]/(mytab["TN"]+mytab["FP"]),3)
stats[aaa,"NPV"]<-round(mytab["TN"]/(mytab["TN"]+mytab["FN"]),3)
stats[aaa,"PPV"]<-round(mytab["TP"]/(mytab["TP"]+mytab["FP"]),3)
stats[aaa,"Accuracy"]<-round((mytab["TP"]+mytab["TN"])/sum(mytab),3)
stats[aaa,"F1"]<-round((2*mytab["TP"])/(2*mytab["TP"]+mytab["FP"]+mytab["FN"]),3)
stats[aaa,"MCC"]<-round((mytab["TP"]*mytab["TN"]-mytab["FP"]*mytab["FN"])/sqrt((mytab["TP"]+mytab["FP"])*(mytab["TP"]+mytab["FN"])*(mytab["TN"]+mytab["FP"])*(mytab["TN"]+mytab["FN"])),3)
stats[aaa,"TP"]<-mytab["TP"]
stats[aaa,"TN"]<-mytab["TN"]
stats[aaa,"FP"]<-mytab["FP"]
stats[aaa,"FN"]<-mytab["FN"]
}
write.table(stats,outfile,quote=F,sep="\t",col.names=NA)
#
}
explore_val(brackenfile=brackenfile,rnaseqfile=rnaseqfile,
smallrnafile=smallrnafile,smallrnacontigfile=smallrnacontigfile,outfile=outfile,remove.failed=remove.failed,
virus=virus)


#https://blog.revolutionanalytics.com/2016/08/roc-curves-in-two-lines-of-code.html
