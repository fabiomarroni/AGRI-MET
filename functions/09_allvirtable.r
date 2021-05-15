# Run with --help flag for help.
# Modified 05/02/2021 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile1"), type="character", default="tables/smallRNAseq_VirDet.txt",
              help="Phase 1 input file, containing viral abundance per sample in matrix format", metavar="character"),
  make_option(c("-i", "--infile2"), type="character", default="tables_val/smallRNAseq_align.txt",
              help="Phase 2 input file, containing viral abundance per sample in matrix format", metavar="character"),
  make_option(c("-V", "--virushostfile"), type="character", default="tables/Virus_Host.txt",
              help="Input file, containing association between virus and its host", metavar="character"),
  make_option(c("-T", "--threshold"), type="numeric", default=0.97,
              help="Threshold to declare presence", metavar="numeric"),
  make_option(c("-A", "--aboveonly"), type="logical", default=TRUE, 
              help="Report only viruses for which at least one sample was positive (i.e. not just >0 but above threshold)  [default= %default]", metavar="character"),
  make_option(c("-n", "--namecol"), type="character", default="smallalign_name", 
              help="Name of column to be used to read names for conversion  [default= %default]", metavar="character"),
  make_option(c("-N", "--namesfile"), type="character", default="Docs_val/Sample_names.txt",
              help="Name conversion table [default= %default]", metavar="character"),
   make_option(c("-O", "--outfile"), type="character", default="tables_val/test.xlsx", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$infile1)) {
  stop("WARNING: No input file from phase 1 specified with '-I' flag.")
} else {  cat ("infile 1 is ", opt$infile1, "\n")
  infile1 <- opt$infile1  
  }

if (is.null(opt$infile2)) {
  stop("WARNING: No input file from phase 2 specified with '-i' flag.")
} else {  cat ("infile 2 is ", opt$infile2, "\n")
  infile2 <- opt$infile2 
  }

if (is.null(opt$namecol)) {
  stop("WARNING: No namecol specified with '-n' flag.")
} else {  cat ("namecol is ", opt$namecol, "\n")
  namecol <- opt$namecol  
  }

if (is.null(opt$threshold)) {
  stop("WARNING: No threshold for declaring presence specified with '-T' flag.")
} else {  cat ("threshold is ", opt$threshold, "\n")
  threshold <- opt$threshold  
  }

  if (is.null(opt$virushostfile)) {
  stop("WARNING: No virus host file specified with '-V' flag.")
} else {  cat ("virushostfile is ", opt$virushostfile, "\n")
  virushostfile <- opt$virushostfile  
  }

  if (is.null(opt$namesfile)) {
  stop("WARNING: No namesfile specified with '-N' flag.")
} else {  cat ("namesfile is ", opt$namesfile, "\n")
  namesfile <- opt$namesfile  
  }

  if (is.null(opt$aboveonly)) {
  stop("WARNING: No aboveonly specified with '-A' flag.")
} else {  cat ("aboveonly is ", opt$aboveonly, "\n")
  aboveonly <- opt$aboveonly  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

do.table<-function(infile1,infile2,namesfile,aboveonly,
                    outfile,virushostfile,threshold,namecol)
{
library(data.table)
library(openxlsx)
indata1<-fread(infile1,data.table=F)
indata2<-fread(infile2,data.table=F)
mynames<-fread(namesfile,data.table=F)
#Clean the names. I remove all the possible prefixes and suffixes
names(indata1)[1]<-"Name"
if(length(grep("Total",names(indata1)))>0) indata1<-indata1[indata1$Total>0,]
indata1$Total<-indata1$taxonomy_id<-indata1$Taxid<-NULL
names(indata1)<-gsub("_FPKM","",names(indata1))
names(indata1)<-gsub("RPM_","",names(indata1))
indata1<-indata1[,grep("Number_",names(indata1),invert=T)]
names(indata1)<-gsub("Percent_","",names(indata1))
#Change to shortnames
names(indata1)[names(indata1)%in%mynames[,namecol]]<-mynames$Short[mynames[,namecol]%in%names(indata1)]
#Clean the names. I remove all the possible prefixes and suffixes
names(indata2)[1]<-"Name"
if(length(grep("Total",names(indata2)))>0) indata2<-indata2[indata2$Total>0,]
indata2$Total<-indata2$taxonomy_id<-indata2$Taxid<-NULL
names(indata2)<-gsub("_FPKM","",names(indata2))
names(indata2)<-gsub("RPM_","",names(indata2))
indata2<-indata2[,grep("Number_",names(indata2),invert=T)]
names(indata2)<-gsub("Percent_","",names(indata2))
#Change to shortnames
nname<-match(names(indata2)[2:ncol(indata2)],mynames[,namecol])
names(indata2)[2:ncol(indata2)]<-mynames$Short[nname]
#Create merged table
indata<-merge(indata1,indata2,by="Name",all=T)
#Remove the "Other" line and Vitis vinifera. We don't need it
indata<-indata[!indata$Name%in%c("Vitis vinifera","Other"),]
indata[is.na(indata)]<-0
indata[,2:ncol(indata)]<-round(indata[,2:ncol(indata)],2)
smeta<-fread(virushostfile,data.table=F)
indata<-merge(indata,smeta,by.x="Name",by.y="name",all.x=T)
row.names(indata)<-indata$Name
indata$Name<-NULL
#Trick to order in the way I like
indata$Number<-factor(indata$Category,levels=c("Grape","Other plants","No plants","PhiX"))
indata<-indata[order(indata$Number,row.names(indata)),]
indata$Number<-NULL
#If aboveonly is true we remove all viruses for which no sample is above threshold
if(aboveonly)
{
indata$pos<-apply(apply(indata[,names(indata)!="Category"],1,">",threshold),2,sum)
indata<-indata[indata$pos>0,]
}

#Format excel worksheet
wb <- createWorkbook()
addWorksheet(wb, "Viruses")
negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE", textDecoration = "Bold")
writeData(wb,"Viruses",indata,rowNames=TRUE)
conditionalFormatting(wb, "Viruses", cols=2:ncol(indata), rows=1:nrow(indata)+1, rule=paste0(">=",threshold), style = posStyle)
saveWorkbook(wb, outfile, TRUE)
}
do.table(infile1=infile1,infile2=infile2,namesfile=namesfile,aboveonly=aboveonly,outfile=outfile,threshold=threshold,virushostfile=virushostfile,namecol=namecol)