# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="clustalw2_out.aln",
              help="Clustalw2 multiple sequence alignment", metavar="character"),
   make_option(c("-O", "--outfile"), type="character", default="ciccio.fasta", 
              help="Trimmed fasta file [default= %default]", metavar="character")
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

#We remove long stretches of "-" at the beginning and end of at least one alignment. 
#This is likely due to the fact that at least one of the sequenced amplicons was shorter. 
#We don't want to use that "error" to infer phylogenies, do we?
remove.gap<-function(infile,
                    outfile)
{
library(seqinr)
myseq<-read.alignment(infile,format="clustal",forceToLower=F)
myseqmat<-as.matrix.alignment(myseq)
pino<-apply(myseqmat,2,table)
gino<-lapply(pino,"[","-")
#When no missing value (-) is present, the search for table name "-" will return an NA.
#The first NA we meet is the first position in which all the sequences are available 
#The last NA is the last position in which all the sequences are avilable
fullset<-which(is.na(gino))
beginpos<-min(fullset)
endpos<-max(fullset)
myseqmat<-myseqmat[,beginpos:endpos]
myseqstring<-apply(myseqmat,1,paste,sep="",collapse="")
write.fasta(as.list(myseqstring),row.names(myseqmat),outfile,as.string=T)
}
remove.gap(infile=infile,outfile=outfile)