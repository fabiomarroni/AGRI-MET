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
   make_option(c("-O", "--outfile"), type="character", default="plots/test.html", 
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

#I read the new sample names as parameter, but for formatting the html table I needed to hardcode the new sample names!!!  
do.table<-function(infile,
                    outfile,virushostfile,threshold,orignames,newnames,
                    top=50,remove.extreme.high=F,nhigh=3,variable.name="Origin",add.col=T,use.new.names=T,
                    new.rotation=T,add.brackets=F)
{
library(data.table)
library(formattable)
library(knitr)
library(htmltools)
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

# my_palette<-alpha(c("orangered","deepskyblue"),0.5)
# my_palette<-c("tomato","deepskyblue")
# my_palette<-c("white","deepskyblue")
#Values equal to the threshold are considered positives.
smeta<-fread(virushostfile,data.table=F)
ff<-merge(pino,smeta,by.x="row.names",by.y="name")
#Manually build the Set1 color palette, because I like it, but I had some problems when there were only two categories.

ff$Order<-2
ff$Order[ff$Category=="Grape"]<-1
ff$Order[ff$Category=="PhiX"]<-4
ff$Order[ff$Category=="No plants"]<-3
ff<-ff[order(ff$Order),]
ff$Order<-NULL
forleg<-ff[!duplicated(ff$Category),]
legvar<-forleg$Category
legcol<-forleg$mycol
#Change grape names
setnames(ff,onames,nnames)

pino<-ff
pino[,nnames]<-round(pino[,nnames],2)
names(pino)[names(pino)=="Row.names"]<-"Virus"
# pino$CO<-round(pino$CO,0)
# pino$MN<-round(pino$MN,0)
# pino$MO<-round(pino$MO,0)
row.names(pino)<-NULL


#	mythreshold<-1572


my_font<-formatter("span",
  style = x ~ style("font-family" = "verdana",
	"font-size" = 10))
font_and_bar <- formatter("span",
  style = x ~ style(
    display = "inline-block",
    direction = "rtl",
    "border-radius" = "4px",
    "padding-right" = "2px",
    "background-color" = csscolor("lightgreen"),
    width = percent(proportion(x)),
    "font-family" = "verdana",
	"font-size" = 10))
font_and_grad <- formatter("span",
  style = x ~ style(
    display = "block",
    direction = "rtl",
    "border-radius" = "4px",
    "padding" = "2px",
    "background-color" = ifelse(x<threshold,"white",csscolor(gradient(as.numeric(x),min.color="yellow",max.color="darkorange"))),
    "font-family" = "verdana",
	"font-size" = 10))
my_cat<- formatter("span", 
    style = ~ style(
	display = "block", #fill the whole cell
	"font-family" = "verdana",
    "padding-right" = "8px",
	"font-size" = 10,
	"background-color" = ifelse(Category=="Grape",alpha("gray100",0.9),ifelse(Category=="Other plants","lightgreen","lightblue")))) #conditional color based on another column
#names(pino)<-my_font(names(pino))
myhtml<-format_table(pino, list(Virus=my_cat,CF=font_and_grad, U1=font_and_grad,U2=font_and_grad,MB=font_and_grad, 
MN=font_and_grad, MO=font_and_grad, RI=font_and_grad, Category=FALSE))

#Change header... This is the only way I found!
myhtml<-gsub("<th style=\"text-align:right","<th style=\"font-family: verdana; padding: 2px 8px ; font-size: 10; text-align:right",myhtml)
write(myhtml,outfile)



}
do.table(infile=infile,outfile=outfile,threshold=threshold,virushostfile=virushostfile,orignames=orignames,newnames=newnames)