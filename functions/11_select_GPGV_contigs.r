# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="smallRNA/2018/02_VirusDetect/result_32_1.trimmed_1.fastq/32_1.trimmed_1.fastq.blastn.xls",
              help="Input file, containing pathway abundance per sample in matrix format", metavar="character"),
   make_option(c("-O", "--outfile"), type="character", default="GPGV/smallRNA_denovo/32_1.trimmed_1.fastq.blastn.fasta", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile)) {
  stop("WARNING: No input folder specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

extract.gpgv<-function(infile,
                    outfile)
{
library(seqinr)
library(data.table)
mydata<-fread(infile,data.table=F)
mydata$classification<-""
mydata$classification[grepl("(?i)(grapevine.*pinot.*gris.*virus)", mydata$Description)]<-"Grapevine Pinot gris virus"
mydata<-mydata[mydata$classification!="",]
mydata<-mydata[!duplicated(mydata$"#Contig_ID"),]
if(nrow(mydata)>0)
{
	write.fasta(as.list(mydata$Contig_Seq),mydata$"#Contig_ID",outfile,as.string=T)
}
}
extract.gpgv(infile=infile,outfile=outfile)