# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gff"), type="character", default="02a_alignments_virusdb_vv/cort_gva_merged.gff3",
              help="Fasta index file", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="02a_alignments_virusdb_vv/cort_gva_merged_TPM_FPKM.txt",
              help="Fasta index file", metavar="character"),
  make_option(c("-T", "--taxfile"), type="character", default="/projects/novabreed/share/software/VirusDetect_v1.7/databases/vrl_genbank.info.gz", 
              help="Output file [default= %default]", metavar="character")
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

if (is.null(opt$gff)) {
  stop("WARNING: No gff file specified with '-G' flag.")
} else {  cat ("gff is ", opt$gff, "\n")
  gff <- opt$gff  
  }

if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

if (is.null(opt$taxfile)) {
  stop("WARNING: No taxfile specified with '-T' flag.")
} else {  cat ("taxfile is ", opt$taxfile, "\n")
  taxfile <- opt$taxfile  
  }

virus_TPM<-function(gff,taxfile,outfile)
{
library("data.table")
#Read gff results and taxonomy files
mydat<-fread(gff,skip=2,data.table=F)
mytax<-fread(paste("gunzip -c",taxfile),data.table=F)
#Remove results on vitis. We don't care about them anymore
mydat<-mydat[grep("chr",mydat$V1,invert=T),]
mydat<-mydat[grep("TPM",mydat$V9),]
mydat$FPKM<-as.numeric(gsub("\"","",(unlist(lapply(strsplit(unlist(lapply(strsplit(mydat$V9," FPKM "),"[",2)),";"),"[",1)))))
mydat$TPM<-as.numeric(gsub("\"","",(unlist(lapply(strsplit(unlist(lapply(strsplit(mydat$V9," TPM "),"[",2)),";"),"[",1)))))
mydat<-mydat[mydat$FPKM>0|mydat$TPM>0,]
mydat<-mydat[,c("V1","FPKM","TPM")]
myres<-merge(mytax,mydat,by="V1")
write.table(myres,outfile,quote=F,sep="\t",row.names=F)
}
virus_TPM(gff=gff,taxfile=taxfile,outfile=outfile)
