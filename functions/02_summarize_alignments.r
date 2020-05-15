# Run with --help flag for help.
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="smallRNA/2018/04_alignments/32_1_alignment_stats.out",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="smallRNA/2018/04_alignments/32_1_alignment_compact.out", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$input)) {
  stop("WARNING: No input file specified with '-I' flag.")
} else {  cat ("Input file is ", opt$input, "\n")
  infile <- opt$input
  }

if (is.null(opt$out)) {
  stop("WARNING: No output file specified with '-O' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  #setwd(wd_location)  
  }

#Merge results obtained by virusdetect and summarize statistics
#I try to put some order in the mess in the virus names (e.g. some entries are uppercase and some lowercase and so on)
#I think the mess in the names is inherited by virusdetect from the ncbi database 
reformat_align<-function(infile=infile,outfile=outfile)
{
library("data.table")
mydat<-fread(infile,data.table=F)
#Regex a manetta!!!
# (?i) case insensitive
# .* separator between words can be any carachter repeated 0 to n times
# | or operator; separate two alternative orders of the words we want to search
mydat$classification<-"Other"
mydat$classification[grepl("(?i)(grapevine.*yellow.*speckle.*viroid|viroid.*yellow.*speckle.*grapevine)", mydat$V3)]<-"Grapevine yellow speckle viroid"
mydat$classification[grepl("(?i)(grapevine.*fleck.*virus)", mydat$V3)]<-"Grapevine fleck virus"
#mydat$classification[grepl("(?i)(grapevine.*asteroid.*mosaic-associated.*virus)", mydat$V3)]<-"Grapevine asteroid mosaic-associated virus"
mydat$classification[grepl("(?i)(rupestris.*vein.*feathering)", mydat$V3)]<-"Grapevine rupestris vein feathering virus"
mydat$classification[grepl("(?i)(grapevine.*fanleaf.*virus)", mydat$V3)]<-"Grapevine fanleaf virus"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 1)", mydat$V3)]<-"Grapevine leafroll-associated virus 1"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 2)", mydat$V3)]<-"Grapevine leafroll-associated virus 2"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 3)", mydat$V3)]<-"Grapevine leafroll-associated virus 3"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 4)", mydat$V3)]<-"Grapevine leafroll-associated virus 4"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 5)", mydat$V3)]<-"Grapevine leafroll-associated virus 5"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 6)", mydat$V3)]<-"Grapevine leafroll-associated virus 6"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 7)", mydat$V3)]<-"Grapevine leafroll-associated virus 7"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 8)", mydat$V3)]<-"Grapevine leafroll-associated virus 8"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 9)", mydat$V3)]<-"Grapevine leafroll-associated virus 9"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 10)", mydat$V3)]<-"Grapevine leafroll-associated virus 10"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 11)", mydat$V3)]<-"Grapevine leafroll-associated virus 11"
mydat$classification[grepl("(?i)(grapevine.*virus A )", mydat$V3)]<-"Grapevine virus A"
mydat$classification[grepl("(?i)(grapevine.*virus B )", mydat$V3)]<-"Grapevine virus B"
mydat$classification[grepl("(?i)(grapevine.*Syrah.*virus)", mydat$V3)]<-"Grapevine Syrah virus"
mydat$classification[grepl("(?i)(arabis.*mosaic.*virus)", mydat$V3)]<-"Arabis mosaic virus"
mydat$classification[grepl("(?i)(grapevine.*asteroid.*mosaic.*virus)", mydat$V3)]<-"Grapevine asteroid mosaic-associated virus"
mydat$classification[grepl("(?i)(hop.*stunt.*viroid)", mydat$V3)]<-"Hop Stunt viroid"
mydat$classification[grepl("(?i)(grapevine.*rootstock.*stem.*lesion.*associated)", mydat$V3)]<-"Grapevine rootstock stem lesion associated virus"
mydat$classification[grepl("(?i)(rupestris.*stem.*pitting.*associated)", mydat$V3)]<-"Grapevine rupestris stem pitting-associated virus"
mydat$classification[grepl("(?i)(grapevine.*pinot.*gris.*virus)", mydat$V3)]<-"Grapevine Pinot gris virus"
mydat$classification[grepl("(?i)(grapevine.*red.*globe.*virus)", mydat$V3)]<-"Grapevine Red Globe virus"
mytab<-aggregate(mydat$V2,by=list(mydat$classification),FUN="sum")
setnames(mytab,c("Classification","Reads"))
write.table(mytab,outfile,sep="\t",quote=F,row.names=F)
}

reformat_align(infile=infile,outfile=outfile)