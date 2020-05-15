# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--infile"), type="character", default="2016-06/04_alignments_vv/cort_gva_2.tab",
              help="Input file", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No input file specified with '-I' flag.")
} else {  cat ("Infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$out)) {
  stop("WARNING: No output file specified with '-I' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  }

#Merge results obtained by virusdetect using blastn and blastx, and summarize statistics  
virus_TPM<-function(infile,outfile)
{
library("data.table")
fullfiles<-unlist(strsplit(infile,","))
for(aaa in 1:length(fullfiles))
{
if(length(grep("summary",fullfiles[aaa]))>0) next
mybase<-gsub(".tab","",basename(fullfiles[aaa]))
mydata<-fread(fullfiles[aaa],data.table=F)
mydata<-mydata[grep("VIT_",mydata[,1],invert=T),]
mydata<-mydata[mydata$TPM>0,]
mydata<-mydata[,c("Reference","Coverage","FPKM","TPM")]
mydata<-aggregate(mydata[,c("Coverage","FPKM","TPM")],by=list(mydata$Reference),FUN="sum")
setnames(mydata,c("Coverage","FPKM","TPM"),paste(mybase,c("Coverage","FPKM","TPM"),sep=""))
ifelse(exists("fulldata"),fulldata<-merge(fulldata,mydata,all=T),fulldata<-mydata)
}
write.table(fulldata,outfile,sep="\t",quote=F,row.names=F)
}
virus_TPM(infile=infile,outfile=outfile)