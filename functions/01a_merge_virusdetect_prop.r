# Run with --help flag for help.
# Modified 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="",
              help="Fasta index file", metavar="character"),
  make_option(c("-P", "--pattern"), type="character", default="_summary.txt",
              help="Pattern", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="",
              help="Fasta index file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$indir)) {
  stop("WARNING: No input directory specified with '-G' flag.")
} else {  cat ("Input directory is ", opt$indir, "\n")
  indir <- opt$indir  
  }

if (is.null(opt$pattern)) {
  stop("WARNING: No pattern specified with '-P' flag.")
} else {  cat ("pattern is ", opt$pattern, "\n")
  pattern <- opt$pattern  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


merge_virus_contigs<-function(indir,pattern,outfile)
{
library("data.table")
setwd(indir)
myfiles<-dir(pattern=pattern)
#Merge all the smallRNA tables
name_codes<-rbind(paste(c("32_1","32_4","32_12"),pattern,sep=""),c("mont_neg_merged","cort_gva_merged","maln_lr3"))
mydat<-fread(myfiles[1],data.table=F)
newname<-name_codes[2,name_codes[1,]==myfiles[1]]
setnames(mydat,c("Number","Percent"),paste(c("Number","Percent"),newname,sep="_"))

for(aaa in 2:length(myfiles))
{
	tdat<-fread(myfiles[aaa],data.table=F)
    newname<-name_codes[2,name_codes[1,]==myfiles[aaa]]
    setnames(tdat,c("Number","Percent"),paste(c("Number","Percent"),newname,sep="_"))
	mydat<-merge(mydat,tdat,by=c("Species","Taxid"),all=T)
}

mydat[is.na(mydat)]<-0
write.table(mydat,outfile,quote=F,sep="\t",row.names=F)

}
merge_virus_contigs(indir=indir,pattern=pattern,outfile=outfile)
