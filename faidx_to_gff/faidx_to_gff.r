# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-X", "--indexfile"), type="character", default="my_faidx.fai",
              help="Fasta index file", metavar="character"),
  make_option(c("-S", "--sourcestring"), type="character", default="ttsg",
              help="Text string indicating the source of the gene prediction", metavar="character"),
  make_option(c("-T", "--gtffile"), type="character", default="my_simple_gtf.gff", 
              help="Output gtf file [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indexfile)) {
  stop("WARNING: No index file specified with '-X' flag.")
} else {  cat ("indexfile is ", opt$indexfile, "\n")
  indexfile <- opt$indexfile  
  }

if (is.null(opt$sourcestring)) {
  stop("WARNING: No source string specified with '-S' flag.")
} else {  cat ("source string is ", opt$sourcestring, "\n")
  sourcestring <- opt$sourcestring  
  }


if (is.null(opt$gtffile)) {
  stop("WARNING: No gtffile file specified with '-T' flag.")
} else {  cat ("gtffile file is ", opt$gtffile, "\n")
  gtffile <- opt$gtffile  
  }

transcriptome_to_simple_gtf<-function(indexfile,gtffile,sourcestring)
{
library("data.table")
idx<-fread(indexfile,data.table=F)
idx<-idx[,1:2]
setnames(idx,"V2","V5")
idx$V2<-sourcestring
idx$V3<-"gene"
idx$V4<-1
idx$V6<-1000
idx$V7<-"+"
idx$V8<-"."
idx$V9<-paste("gene_id=",paste(idx$V1,".1",sep=""),"; Parent=",idx$V1,sep="")
idx<-idx[,sort(names(idx))]
idx$buttami<-seq(1:nrow(idx))
newidx<-idx
newidx$V3<-"exon"
idx<-rbind(idx,newidx)
idx<-idx[order(idx$buttami),]
idx$buttami<-NULL
write.table(idx,gtffile,quote=F,sep="\t",row.names=F,col.names=F)
}
transcriptome_to_simple_gtf(indexfile=indexfile,gtffile=gtffile,sourcestring=sourcestring)
