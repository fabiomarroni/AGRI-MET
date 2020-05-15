# Run with --help flag for help.
# Modified 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--indir"), type="character", default="02a_alignments_virusdb_vv",
              help="Fasta index file", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="tables/RNAseq_VirDet.txt",
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

if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


merge_virus_TPM<-function(indir,outfile)
{
library("data.table")
setwd(indir)
#Merge all the TPM and FPKM tables
myfiles<-dir(pattern="_TPM_FPKM.txt")
mysamples<-gsub("_TPM_FPKM.txt","",myfiles)
mydat<-fread(myfiles[1],data.table=F)
setnames(mydat,c("TPM","FPKM"),c(paste(mysamples[1],"TPM",sep="_"),paste(mysamples[1],"FPKM",sep="_")))
for(aaa in 2:length(myfiles))
{
	tdat<-fread(myfiles[aaa],data.table=F)
	setnames(tdat,c("TPM","FPKM"),c(paste(mysamples[aaa],"TPM",sep="_"),paste(mysamples[aaa],"FPKM",sep="_")))
	mydat<-merge(mydat,tdat,by=c("V1","V2","V3","V4","V5","V6","V7"),all=T)
}
#Remove TPM. Right now I am interested in FPKM.
mydat<-mydat[,grep("TPM",names(mydat),invert=T)]
mydat[is.na(mydat)]<-0
mydat$total<-rowSums(mydat[,grep("FPKM",names(mydat))])

#Regex a manetta!!!
# (?i) case insensitive
# .* separator between words can be any carachter repeated 0 to n times
# | or operator; separate two alternative orders of the words we want to search
mydat$classification<-"Other"
mydat$classification[grepl("(?i)(grapevine.*yellow.*speckle.*viroid|viroid.*yellow.*speckle.*grapevine)", mydat$V4)]<-"Grapevine yellow speckle viroid 1"
mydat$classification[grepl("(?i)(grapevine.*fleck.*virus)", mydat$V4)]<-"Grapevine fleck virus"
mydat$classification[grepl("(?i)(grapevine.*asteroid.*mosaic-associated.*virus)", mydat$V4)]<-"Grapevine asteroid mosaic-associated virus"
mydat$classification[grepl("(?i)(rupestris.*vein.*feathering)", mydat$V4)]<-"Grapevine rupestris vein feathering virus"
mydat$classification[grepl("(?i)(grapevine.*fanleaf.*virus)", mydat$V4)]<-"Grapevine fanleaf virus"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus.*Carn)", mydat$V4)]<-"Grapevine leafroll-associated virus Carn"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 1)", mydat$V4)]<-"Grapevine leafroll-associated virus 1"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 2)", mydat$V4)]<-"Grapevine leafroll-associated virus 2"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 3)", mydat$V4)]<-"Grapevine leafroll-associated virus 3"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 4)", mydat$V4)]<-"Grapevine leafroll-associated virus 4"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 5)", mydat$V4)]<-"Grapevine leafroll-associated virus 5"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 6)", mydat$V4)]<-"Grapevine leafroll-associated virus 6"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 7)", mydat$V4)]<-"Grapevine leafroll-associated virus 7"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 8)", mydat$V4)]<-"Grapevine leafroll-associated virus 8"
mydat$classification[grepl("(?i)(grapevine.*leafroll.*associated.*virus 9)", mydat$V4)]<-"Grapevine leafroll-associated virus 9"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 10)", mydat$V4)]<-"Grapevine leafroll-associated virus 10"
mydat$classification[grepl("(?i)(grapevine.*leafroll-associated.*virus 11)", mydat$V4)]<-"Grapevine leafroll-associated virus 11"
mydat$classification[grepl("(?i)(grapevine.*virus A )", mydat$V4)]<-"Grapevine virus A"
mydat$classification[grepl("(?i)(grapevine.*virus B )", mydat$V4)]<-"Grapevine virus B"
mydat$classification[grepl("(?i)(grapevine.*virus D )", mydat$V4)]<-"Grapevine virus D"
mydat$classification[grepl("(?i)(grapevine.*virus E )", mydat$V4)]<-"Grapevine virus E"
mydat$classification[grepl("(?i)(grapevine.*virus F )", mydat$V4)]<-"Grapevine virus F"
mydat$classification[grepl("(?i)(grapevine.*virus Q )", mydat$V4)]<-"Grapevine virus Q"
mydat$classification[grepl("(?i)(grapevine.*Syrah.*virus.*1)", mydat$V4)]<-"Grapevine Syrah virus 1"
mydat$classification[grepl("(?i)(arabis.*mosaic.*virus)", mydat$V4)]<-"Arabis mosaic virus"
mydat$classification[grepl("(?i)(grapevine.*asteroid.*mosaic.*virus)", mydat$V4)]<-"Grapevine asteroid mosaic-associated virus"
mydat$classification[grepl("(?i)(hop.*stunt.*viroid)", mydat$V4)]<-"Hop stunt viroid"
mydat$classification[grepl("(?i)(grapevine.*rootstock.*stem.*lesion.*associated)", mydat$V4)]<-"Grapevine rootstock stem lesion associated virus"
mydat$classification[grepl("(?i)(rupestris.*stem.*pitting.*associated)", mydat$V4)]<-"Grapevine rupestris stem pitting-associated virus"
mydat$classification[grepl("(?i)(grapevine.*pinot.*gris.*virus)", mydat$V4)]<-"Grapevine Pinot gris virus"
mydat$classification[grepl("(?i)(grapevine.*red.*globe.*virus)", mydat$V4)]<-"Grapevine Red Globe virus"
mydat$classification[grepl("(?i)(grapevine.*satellite.*virus)", mydat$V4)]<-"Grapevine satellite virus"
mydat$classification[grepl("(?i)(grapevine.*deformation.*virus)", mydat$V4)]<-"Grapevine deformation virus"
mydat$classification[grepl("(?i)(grapevine.*associated.*narnavirus.*1)", mydat$V4)]<-"Grapevine associated narnavirus-1"


#Insert viruses not belonging to grapevine
mydat$classification[grepl("(?i)(Andean.*potato.*latent.*virus)", mydat$V4)]<-"Andean potato latent virus"
mydat$classification[grepl("(?i)(Erysimum.*latent.*virus)", mydat$V4)]<-"Erysimum latent virus"
mydat$classification[grepl("(?i)(UNVERIFIED.*Cactus.*virus.*X)", mydat$V4)]<-"Unverified cactus virus X"
mydat$classification[grepl("(?i)(Turnip.*yellow.*mosaic.*virus)", mydat$V4)]<-"Turnip yellow mosaic virus"
mydat$classification[grepl("(?i)(Okra.*mosaic.*virus)", mydat$V4)]<-"Okra mosaic virus"
mydat$classification[grepl("(?i)(Maize.*rayado.*fino.*virus)", mydat$V4)]<-"Maize rayado fino virus"
mydat$classification[grepl("(?i)(Blackberry.*virus.*S)", mydat$V4)]<-"Blackberry virus S"
mydat$classification[grepl("(?i)(Physalis.*mottle.*tymovirus)", mydat$V4)]<-"Physalis mottle tymovirus"
mydat$classification[grepl("(?i)(Poinsettia.*mosaic.*virus)", mydat$V4)]<-"Poinsettia mosaic virus"
mydat$classification[grepl("(?i)(Apple.*mosaic.*virus)", mydat$V4)]<-"Apple mosaic virus"
mydat$classification[grepl("(?i)(Apple.*stem.*pitting.*virus)", mydat$V4)]<-"Apple stem pitting virus"
mydat$classification[grepl("(?i)(Little.*cherry.*virus.*1)", mydat$V4)]<-"Little cherry virus 1"
mydat$classification[grepl("(?i)(Twisted.*stalk.*leaf.*streak.*virus)", mydat$V4)]<-"Twisted-stalk leaf streak virus"
mydat$classification[grepl("(?i)(Fig.*fleck.*associated.*virus)", mydat$V4)]<-"Fig fleck-associated virus"
mydat$classification[grepl("(?i)(Mertensia.*leaf.*curl.*virus)", mydat$V4)]<-"Mertensia leaf curl virus"
mydat$classification[grepl("(?i)(Olive.*latent.*virus.*3)", mydat$V4)]<-"Olive latent virus 3"
mydat$classification[grepl("(?i)(Eggplant.*mosaic.*virus.*3)", mydat$V4)]<-"Eggplant mosaic virus"
mydat$classification[grepl("(?i)(Cassia.*yellow.*mosaic.*associated)", mydat$V4)]<-"Cassia yellow mosaic-associated virus"
mydat$classification[grepl("(?i)(Citrus.*sudden.*death.*associated.*virus)", mydat$V4)]<-"Citrus sudden death-associated virus"
mydat$classification[grepl("(?i)(Citrus.*bent.*leaf.*viroid)", mydat$V4)]<-"Citrus bent leaf viroid"
mydat$classification[grepl("(?i)(Citrus.*psorosis.*virus)", mydat$V4)]<-"Citrus psorosis virus"
mydat$classification[grepl("(?i)(Diascia.*yellow.*mottle.*virus)", mydat$V4)]<-"Diascia yellow mottle virus"
mydat$classification[grepl("(?i)(Dulcamara.*mottle.*virus)", mydat$V4)]<-"Dulcamara mottle virus"
mydat$classification[grepl("(?i)(Bamboo.*mosaic.*virus)", mydat$V4)]<-"Bamboo mosaic virus"
mydat$classification[grepl("(?i)(Soybean.*associated.*ourmiavirus.*2)", mydat$V4)]<-"Soybean-associated ourmiavirus 2"
mydat$classification[grepl("(?i)(Soybean.*associated.*endornavirus.*1)", mydat$V4)]<-"Soybean-associated endornavirus 1"
mydat$classification[grepl("(?i)(*Oat.*blue.*dwarf.*virus)", mydat$V4)]<-"Oat blue dwarf virus"
mydat$classification[grepl("(?i)(*potato.*virus.*T)", mydat$V4)]<-"Potato virus T"
mydat$classification[grepl("(?i)(Discula.*destructiva.*virus.*1)", mydat$V4)]<-"Discula destructiva virus 1"
mydat$classification[grepl("(?i)(Maize.*stripe.*virus)", mydat$V4)]<-"Maize stripe virus"
mydat$classification[grepl("(?i)(Hop.*trefoil.*cryptic.*virus.*1)", mydat$V4)]<-"Hop trefoil cryptic virus 2"
mydat$classification[grepl("(?i)(Heterobasidion.*partitivirus.*17)", mydat$V4)]<-"Heterobasidion partitivirus 17"
mydat$classification[grepl("(?i)(Botrytis.*cinerea.*mitovirus.*2)", mydat$V4)]<-"Botrytis cinerea mitovirus 2"

findat<-aggregate(mydat[,grep("FPKM",names(mydat))],by=list(mydat$classification),FUN=sum)
findat$Total<-rowSums(findat[,grep("FPKM",names(findat))])
#Only report results for which FPKM/TPM is greater than 5 across samples
#This step has been omitted, and removed from M&M in the paper
#findat<-findat[findat$Total>=5,]
setnames(findat,"Group.1","Name")
#Change names of small-RNA data
if(length(grep("32_1",names(findat)))>0)
{
setnames(findat,c("32_1.trimmed_FPKM","32_12.trimmed_FPKM","32_4.trimmed_FPKM"),c("mont_neg_merged","maln_lr3","cort_gva_merged"))
}
write.table(findat,outfile,quote=F,sep="\t",row.names=F)


}
merge_virus_TPM(indir=indir,outfile=outfile)
