# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="smallRNA_val/02_webVirusDetect/result_55311ID1718COU121S201.fastq",
              help="Input directory", metavar="character"),
  make_option(c("-N", "--namefile"), type="character", default="", 
              help="File assigning species name to taxid [default= %default]", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="smallRNA_val/03_smallRNA_tables/55311ID1718COU121S201.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-S", "--summary"), type="character", default="smallRNA_val/03_smallRNA_tables/55311ID1718COU121S201_summary.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-T", "--taxfile"), type="character", default="", 
              help="Two column file with in column 1 seqid and in column 2 taxid [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$input)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Working directory is ", opt$input, "\n")
  input <- opt$input  
  }

if (is.null(opt$namefile)) {
  stop("WARNING: No namefile specified with '-N' flag.")
} else {  cat ("Namefile is ", opt$namefile, "\n")
  namefile <- opt$namefile  
  }

if (is.null(opt$summary)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output file is ", opt$summary, "\n")
  summary_file <- opt$summary  
  }

  if (is.null(opt$taxfile)) {
  stop("WARNING: No taxfile specified with '-T' flag.")
} else {  cat ("taxfile is ", opt$taxfile, "\n")
  taxfile <- opt$taxfile  
  }

  if (is.null(opt$out)) {
  stop("WARNING: No input directory specified with '-I' flag.")
} else {  cat ("Output file is ", opt$out, "\n")
  outfile <- opt$out  
  }

#Merge results obtained by virusdetect using blastn and blastx, and summarize statistics  
small_virus<-function(indir=input,outfile=outfile,namefile=namefile,taxfile=taxfile,summary_file=summary_file,loop=FALSE)
{
library("data.table")
myfiles<-dir(indir,full.names=T,pattern=".xls")
mysmall<-fread(myfiles[1],data.table=F)
if(length(myfiles)>1)
{
	for (aaa in 2:length(myfiles))
	{
	mysmall<-rbind(mysmall,fread(myfiles[aaa],data.table=F))
	}
}
mysmall<-mysmall[,!names(mysmall)%in%c("Contig_Seq","Hsp_strand")]
mytab<-table(mysmall$Description)
write.table(mytab,outfile,sep="\t",quote=F,row.names=F)
#Start second part of the function. We use the seqid to find the taxid and unambiguously associate the species to each contig 
#At present we are only associating the taxon. 
mysmall<-mysmall[,c("Hit_ID","Description")]
mytax<-fread(taxfile,data.table=F)
mytab<-data.frame(table(mysmall$Hit_ID),stringsAsFactors=F)
mytab$Var1<-as.character(mytab$Var1)

mytax$Seqid<-unlist(lapply(strsplit(mytax$V1,"\\."),"[",1))
smalltax<-merge(mytab,mytax,by.x="Var1",by.y="Seqid",all.x=T,all.y=F)
#rm(mytax)
setnames(smalltax,"V2","Taxid")
smalltax<-na.omit(smalltax)
smalltax<-aggregate(smalltax$Freq,by=list(smalltax$Taxid),FUN="sum")
setnames(smalltax,c("Taxid","Number"))

mynames<-fread(namefile,data.table=F)
#Remove tabs that remained in some columns
mynames$V2<-gsub("\t","",mynames$V2)
mynames$V3<-gsub("\t","",mynames$V3)
mynames$V4<-gsub("\t","",mynames$V4)
selnames<-mynames[mynames$V1%in%smalltax$Taxid,]
#Other tricks to deal with changes in the format of files or behaviour of fread.
#In the old version we had a V7 column with the taxonomy description, noe we have it in V4.
ifelse(sum(names(selnames)=="V7")>0,selnames<-selnames[selnames$V7=="scientific name",],selnames<-selnames[selnames$V4=="scientific name",]) 
#In the old version the species name was in column V3, now it is in column V4.
ifelse(length(table(selnames$V3))>length(table(selnames$V2)),selnames<-selnames[,c("V1","V3")],selnames<-selnames[,c("V1","V2")])
names(selnames)[2]<-"Species"
taxnames<-merge(smalltax,selnames,by.x="Taxid",by.y="V1",all.x=T,all.y=F)
taxnames$Percent<-round(100*taxnames$Number/sum(taxnames$Number),2)
write.table(taxnames,summary_file,sep="\t",quote=F,row.names=F)
}

small_virus(indir=input,outfile=outfile,taxfile=taxfile,summary_file=summary_file,namefile)
