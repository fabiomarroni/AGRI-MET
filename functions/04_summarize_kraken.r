# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="RNAseq/01a_kraken206",
              help="Input directory containing kraken reports", metavar="character"),
  make_option(c("-V", "--virtax"), type="character", default="RNAseq/01a_kraken207/Viruses.txt",
              help="File with viruses taxid", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="RNAseq/03_RNA_res/kraken_standard_206.txt", 
              help="output file name [default= %default]", metavar="character")
  # make_option(c("-T", "--treated_prefix"), type="character", default="experimental_",
              # help="Prefix identifying treated samples", metavar="character"),
  # make_option(c("-C", "--control_prefix"), type="character", default="control_",
              # help="Prefix identifying control samples", metavar="character"),
  # make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              # help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input folder specified with '-I' flag.")
} else {  cat ("Indir is ", opt$indir, "\n")
  indir <- opt$indir  
  }

if (is.null(opt$virtax)) {
  stop("WARNING: No virus taxonomy file specified with '-V' flag.")
} else {  cat ("Virus taxonomy file is ", opt$virtax, "\n")
  virtax <- opt$virtax  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No output file specified with '-I' flag.")
} else {  cat ("Output file is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

#Merge results obtained by virusdetect using blastn and blastx, and summarize statistics  
virus_kraken<-function(indir,outfile,virtax,level="S")
{
library("data.table")
setwd(indir)
fullfiles<-dir(pattern="report")
for(aaa in 1:length(fullfiles))
{
	mydata<-fread(fullfiles[aaa],data.table=F)
	mydata<-mydata[mydata$V4=="U"|mydata$V4==level,]
	mydata$RPM<-1000000*mydata$V2/sum(mydata$V2)
	mydata<-mydata[,c("V6","V5","RPM")]
	mydata<-mydata[mydata$V6!="unclassified",]
	setnames(mydata,"RPM",paste("RPM",gsub(".kraken.report.txt","",fullfiles[aaa]),sep="_"))
	if(aaa==1) findata<-mydata
	if(aaa>1) findata<-merge(findata,mydata,by=c("V6","V5"),all=TRUE)
}
findata[is.na(findata)]<-0
findata$Total<-rowSums(findata[,3:ncol(findata)])
findata<-findata[order(findata$Total,decreasing=T),]
findata<-findata[findata$Total>10,]
setnames(findata,c("V6","V5"),c("Name","Taxid"))
vtax<-fread(virtax,data.table=F)
virdata<-findata[findata$Taxid%in%as.character(unlist(vtax)),]
write.table(virdata,outfile,sep="\t",quote=F,row.names=F)
#Still need to read virus file RNAseq/01a_kraken207/Viruses.txt
#And remove all non-virus taxa
}
virus_kraken(indir=indir,virtax=virtax,outfile=outfile)