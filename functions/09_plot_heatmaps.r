# Run with --help flag for help.
# Modified 02/05/2018 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})


option_list = list(
  make_option(c("-I", "--infile"), type="character", default="tables/smallRNAseq_VirDet.txt",
              help="Input file, containing pathway abundance per sample in matrix format", metavar="character"),
  make_option(c("-V", "--virushostfile"), type="character", default="tables/Virus_Host.txt",
              help="Input file, containing association between virus and its host", metavar="character"),
  make_option(c("-T", "--threshold"), type="numeric", default=50,
              help="Threshold to declare presence", metavar="numeric"),
   make_option(c("-R", "--orignames"), type="character", default="", 
              help="Cultivar names as they are written in the input file [default= %default]", metavar="character"),
   make_option(c("-N", "--newnames"), type="character", default="", 
              help="Cultivar namesas they need to be plotted on the output file [default= %default]", metavar="character"),
   make_option(c("-O", "--outfile"), type="character", default="plots/test.pdf", 
              help="output file name [default= %default]", metavar="character")
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

if (is.null(opt$infile)) {
  stop("WARNING: No input folder specified with '-I' flag.")
} else {  cat ("infile is ", opt$infile, "\n")
  infile <- opt$infile  
  }

if (is.null(opt$orignames)) {
  stop("WARNING: No orignames for declaring presence specified with '-T' flag.")
} else {  cat ("orignames is ", opt$orignames, "\n")
  orignames <- opt$orignames  
  }

if (is.null(opt$newnames)) {
  stop("WARNING: No newnames for declaring presence specified with '-T' flag.")
} else {  cat ("newnames is ", opt$newnames, "\n")
  newnames <- opt$newnames  
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

if (is.null(opt$outfile)) {
  stop("WARNING: No outfile specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }

  
  #Perform clustering using data obtained from the function above (i.e. each taxonomy level in a different column)
#As a default we use the standard genus species taxonomy
#IMPORTANT: The higher taxonomy level should be the first, and the lower the second!!!
do.cluster<-function(infile,
                    outfile,virushostfile,threshold,orignames,newnames,
                    top=50,remove.extreme.high=F,nhigh=3,variable.name="Origin",add.col=T,use.new.names=T,
                    new.rotation=T,add.brackets=F)
{
library(data.table)
library(gplots)
library(RColorBrewer)
library(scales)
onames<-unlist(strsplit(orignames,","))
nnames<-unlist(strsplit(newnames,","))
if(length(onames)!=length(nnames)) stop("The length of the original names differs from the length of the new names")
mydat<-fread(infile,data.table=F)
mydat$taxonomy_id<-mydat$Total<-NULL
#In case we use the smallRNA results, the first columns is not "name", but "Species"
names(mydat)<-gsub("Species","name",names(mydat))
#Sometimes the "name" column is uppercase and we don't like it
names(mydat)<-gsub("Name","name",names(mydat))
#In addition we remove the colums including the number of contigs, and only use the percentage
mydat<-mydat[,grep("Number_",names(mydat),invert=T)]
mydat<-mydat[,grep("Number_",names(mydat),invert=T)]
#Some files also have a Taxid column that we do not want
mydat<-mydat[,!names(mydat)%in%"Taxid"]
row.names(mydat)<-mydat$name
mydat$name<-NULL
#Sort alphabetically 
mydat<-mydat[,sort(names(mydat))]
pino<-mydat
# length_palette<-8
# my_palette <- colorRampPalette(c("lightgoldenrodyellow","yellow","red"))(n=length_palette)
my_palette<-alpha(c("orangered","deepskyblue"),0.5)
my_palette<-c("tomato","deepskyblue")
my_palette<-c("white","deepskyblue")
#lightcoral
#Values equal to the threshold are considered positives.
my_breaks<-c(min(na.omit(as.vector(pino))),threshold-0.01,max(na.omit(as.vector(pino))))
pdf(outfile,width=6)
par(mar=c(0,0,0,0))
# pp<-data.frame(pino)
# names(pp)<-"name"
#Associate viruses to their host
smeta<-fread(virushostfile,data.table=F)
#browser()
ff<-merge(pino,smeta,by.x="row.names",by.y="name")
#Manually build the Set1 color palette, because I like it, but I had some problems when there were only two categories.
ff$mycol<-"#4DAF4A"
ff$mycol[ff$Category=="Grape"]<-"#E41A1C"
ff$mycol[ff$Category=="PhiX"]<-"#984EA3"
ff$mycol[ff$Category=="No plants"]<-"#377EB8"

ff$Order<-2
ff$Order[ff$Category=="Grape"]<-1
ff$Order[ff$Category=="PhiX"]<-4
ff$Order[ff$Category=="No plants"]<-3
ff<-ff[order(ff$Order),]
ff$Order<-NULL
forleg<-ff[!duplicated(ff$Category),]
legvar<-forleg$Category
legcol<-forleg$mycol
mycol<-ff$mycol
#Change grape names
setnames(ff,onames,nnames)

pino<-ff
row.names(pino)<-pino$Row.names
pino$Category<-pino$mycol<-pino$Row.names<-NULL
if(max(pino)<=10000) pino<-as.matrix(round(pino,2))
if(max(pino)>10000) pino<-as.matrix(round(pino,0))
#pino[pino==0]<-NA
heatmap.2(pino,  #Transposed matrix for better readability
            #Colv=FALSE, #No reordering on columns
			cellnote=pino, #Decide what to write in each cell
            notecol="black", #Color for characters in notes
            na.color=my_palette[1], #Color for plotting NA. I just did it because I don't like to have a lot of zeros written in the heatmap
            notecex=0.9, #Cex for notes
            col=my_palette,
            sepcolor=my_palette,
            sepwidth=c(0,0),
			Rowv=FALSE, #do not reorder rows
            Colv=FALSE, #do not reorder cols
            breaks=my_breaks,			
			RowSideColors=mycol,
            #density.info="none", #No density
            trace="none", #No trace
            dendrogram="none",
			#lmat=lmat,
            cexRow=0.8, #Resize row labels 
            cexCol=0.8, #Resize COl labels
			lhei=c(0.2,4),
            lwid=c(0.4,4),
			margins=c(9,16), #Change margins
            #margins=c(24,18), #Change margins
			key=F, #No color key
			#keysize=0.5 #Size of color key
            )
#legend(x=0.05,y=0.2,legend=legvar,col=legcol,lty=1,lwd=10,cex=0.8)
legend(x=0.08,y=0.12,legend=legvar,fill=legcol,cex=0.8,horiz=T,title="Viral host",bty="n")
#Plot color matrix
#colmat<-matrix(rep(rev(my_breaks),n=length(my_breaks)),nrow=length(my_breaks),ncol=length(my_breaks),byrow=T)
#Color legend
# image(x=seq(0.3,0.9,length.out=length(my_breaks)),y=seq(0.95,0.97,length.out=length(my_breaks)),z=t(colmat),add=T)
# text(x=0.3,y=0.9,"0",cex=0.8)
# text(x=0.9,y=0.9,round(max(pino),0),cex=0.8)
# text(x=0.55,y=0.9,bquote(log[10]~"OTU abundance"),cex=0.8)
dev.off()
}
do.cluster(infile=infile,outfile=outfile,threshold=threshold,virushostfile=virushostfile,orignames=orignames,newnames=newnames)