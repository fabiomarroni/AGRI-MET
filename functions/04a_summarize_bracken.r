# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="RNAseq/01a_kraken_VD206",
              help="Input directory containing kraken reports", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="tables/RNA_bracken_VD_206.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input folder specified with '-I' flag.")
} else {  cat ("Indir is ", opt$indir, "\n")
  indir <- opt$indir  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No output file specified with '-I' flag.")
} else {  cat ("Output file is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

#Merge results obtained by virusdetect using blastn and blastx, and summarize statistics  
virus_kraken<-function(indir,outfile,level="S")
{
library("data.table")
setwd(indir)
fullfiles<-dir(pattern="_S.bracken.txt")
krakfiles<-dir(pattern="kraken.report_bracken.txt")
for(aaa in 1:length(fullfiles))
{
	myunk<-fread(krakfiles[aaa],data.table=F)
	mykperc<-myunk$V1[myunk$V6=="root"]/100
	mydata<-fread(fullfiles[aaa],data.table=F)
	mydata$RPM<-mykperc*1000000*mydata$fraction_total_reads
	mydata<-mydata[,c("name","taxonomy_id","RPM")]
	setnames(mydata,"RPM",paste("RPM",gsub("_S.bracken.txt","",fullfiles[aaa]),sep="_"))
	if(aaa==1) findata<-mydata
	if(aaa>1) findata<-merge(findata,mydata,by=c("name","taxonomy_id"),all=TRUE)
}
findata[is.na(findata)]<-0
findata$Total<-rowSums(findata[,3:ncol(findata)])
findata<-findata[order(findata$Total,decreasing=T),]
#findata<-findata[findata$Total>40,]
write.table(findata,outfile,sep="\t",quote=F,row.names=F)
}
virus_kraken(indir=indir,outfile=outfile)