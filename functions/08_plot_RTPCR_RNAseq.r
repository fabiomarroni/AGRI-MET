# Run with --help flag for help.
# Modified 02/05/2019 by Fabio Marroni
suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-B", "--brackenfile"), type="character", default="tables/RNA_bracken_nt_206.txt",
              help="Bracken summary file", metavar="character"),
  make_option(c("-I", "--includesmall"), type="logical", default="FALSE",
              help="Should results be shown only for samples that include smallRNA seq data?", metavar="logical"),
  make_option(c("-P", "--pcrfile"), type="character", default="tables/RT-PCR.txt",
              help="RT-PCR file", metavar="character"),
  make_option(c("-E", "--elisafile"), type="character", default="tables/Elisa_res.txt",
              help="ELISA file", metavar="character"),
  make_option(c("-R", "--rnaseqfile"), type="character", default="tables/RNAseq_VirDet.txt",
              help="Fasta index file", metavar="character"),
  make_option(c("-S", "--smallrnafile"), type="character", default="tables/smallRNAseq_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-C", "--smallrnacontigfile"), type="character", default="tables/smallRNA_contigs_VirDet.txt",
              help="Small RNA file", metavar="character"),
  make_option(c("-O", "--outfile"), type="character", default="tables/RNA_smallRNA_RTPCR_gold_std.txt",
              help="Fasta index file", metavar="character"),
  make_option(c("-G", "--graphdir"), type="character", default="plots/",
              help="Fasta index file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#print("USAGE: $ 01_small_RNA.r -I input dir/ -O summarized output file")

if (is.null(opt$brackenfile)) {
  stop("WARNING: No brackenfile specified with '-B' flag.")
} else {  cat ("Bracken summary file is ", opt$brackenfile, "\n")
  brackenfile <- opt$brackenfile  
  }

if (is.null(opt$includesmall)) {
  stop("WARNING: No includesmall parameter specified with '-I' flag.")
} else {  cat ("includesmall is ", opt$includesmall, "\n")
  includesmall <- opt$includesmall  
  }

if (is.null(opt$pcrfile)) {
  stop("WARNING: No pcrfile specified with '-P' flag.")
} else {  cat ("PCR file is ", opt$pcrfile, "\n")
  pcrfile <- opt$pcrfile  
  }

  if (is.null(opt$elisafile)) {
  stop("WARNING: No elisafile specified with '-E' flag.")
} else {  cat ("ELISA file is ", opt$elisafile, "\n")
  elisafile <- opt$elisafile  
  }

if (is.null(opt$rnaseqfile)) {
  stop("WARNING: No rnaseqfile specified with '-R' flag.")
} else {  cat ("RNA seq file is ", opt$rnaseqfile, "\n")
  rnaseqfile <- opt$rnaseqfile  
  }

if (is.null(opt$smallrnafile)) {
  stop("WARNING: No smallrnafile specified with '-S' flag.")
} else {  cat ("smallrnafile file is ", opt$smallrnafile, "\n")
  smallrnafile <- opt$smallrnafile  
  }

if (is.null(opt$smallrnacontigfile)) {
  stop("WARNING: No smallrnacontigfile specified with '-C' flag.")
} else {  cat ("smallrnacontigfile file is ", opt$smallrnacontigfile, "\n")
  smallrnacontigfile <- opt$smallrnacontigfile  
  }

  if (is.null(opt$graphdir)) {
  stop("WARNING: No graphdir specified with '-G' flag.")
} else {  cat ("graphdir is ", opt$graphdir, "\n")
  graphdir <- opt$graphdir  
  }

  if (is.null(opt$outfile)) {
  stop("WARNING: No out file specified with '-O' flag.")
} else {  cat ("outfile is ", opt$outfile, "\n")
  outfile <- opt$outfile  
  }


  
simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

#This is the rocit version of plot modified by me to make it more customizable
plot.rocit <- function(x, col = c("#2F4F4F", "#BEBEBE"),
                       legend = TRUE, legendpos = "bottomright",
                       YIndex = TRUE, values = TRUE, add.x.lab=FALSE, add.y.lab=FALSE, xtitle=NULL, ytitle=NULL, ... = NULL,add.title=F)
{

  col <- rep(col, 2)
  y1 <- x$TPR
  x1 <- x$FPR
  plot(x1, y1, type = "l",
       ylab = "Sensitivity (TPR)",
       xlab = "1-Specificity (FPR)",
       col = col[1], lwd = 2, axes = FALSE )
  ifelse(add.x.lab,axis(side=1,labels=seq(0,1,0.2),at=seq(0,1,0.2)),axis(side=1,labels=FALSE))
  ifelse(add.y.lab,axis(side=2,labels=seq(0,1,0.2),at=seq(0,1,0.2)),axis(side=2,labels=FALSE))
  if(add.title) title(ylab = "Sensitivity (TPR)",xlab = "1-Specificity (FPR)",outer=TRUE)
  grid(col = "gray60")
  abline(0, 1, lwd = 2, col = col[2], lty = 2)
  diff <- y1 - x1
  maxIndex <- which.max(diff)
  xYI <- x1[maxIndex]
  yYI <- y1[maxIndex]
  cYI <- x$Cutoff[maxIndex]
  if(YIndex){
    points(x = xYI, y = yYI, pch = 16, cex = 1)
    points(x = xYI, y = yYI, pch = 3, cex = 3)
    text(x = (xYI + 0.2), y = (yYI - 0.1),
         "Optimal (Youden Index) point", font = 3)
  }


  methodName <- x$method
  methodName <- paste0(toupper(substr(methodName, 1, 1)),
                       substr(methodName, 2, nchar(methodName)))
  maintext <- paste(methodName, "ROC plot")

  if(legend){
    legend(legendpos,
           c(paste(methodName, "ROC curve"), "Chance line"),
           lty = c(1, 2), lwd = 1.5, col = col,cex=0.9)
  }

  if(values){
    return(invisible(list(method = x$method,
                          AUC = x$AUC,
                          Cutoff = x$Cutoff,
                          `TPR/Sensitivity/Recall` = x$TPR,
                          `FPR/1-Specificity` = x$FPR,
                          `optimal Youden Index point` =
                            c(value = max(diff),
                              FPR = xYI,
                              TPR = yYI,
                              cutoff = cYI))))


  }
}
  
RTPCR_RNAseq<-function(pcrfile,elisafile,rnaseqfile,smallrnafile,smallrnacontigfile,brackenfile,outfile,graphdir,includesmall,col.list=c("black","grey"))
{
library("data.table")
library(scales)
library(RColorBrewer)
library(ROCit)
#We assume that the files have the same numbe of columns and with the same name
#Read results of RTPCR 
PCRdat<-fread(pcrfile,data.table=F,fill=T)
rownames(PCRdat)<-PCRdat[,1]
PCRdat[,1]<-NULL
PCRdat<-PCRdat[,sort(names(PCRdat))]
row.names(PCRdat)<-gsub("_FPKM","",row.names(PCRdat))
PCRdat[is.na(PCRdat)]<-0
#Read ELISA results
ELISAdat<-fread(elisafile,data.table=F)
ELISAdat$Code<-NULL
row.names(ELISAdat)<-ELISAdat$"Sample name"
ELISAdat$"Sample name"<-NULL

#Read results of RNA alignment 
#Change names of Viruses: only viruses with short name corresponding to those present in RT or ELISA files are included in analysis.
RNAdat<-fread(rnaseqfile,data.table=F)
RNAdat$short<-NA
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
RNAdat$short[RNAdat$Name=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
RNAdat$short[RNAdat$Name=="Grapevine rupestris stem pitting-associated virus"]<-"RPSaV"
RNAdat$short[RNAdat$Name=="Grapevine Pinot gris virus"]<-"GPGV"
RNAdat$short[RNAdat$Name=="Grapevine virus A"]<-"GVA"
RNAdat$short[RNAdat$Name=="Grapevine fanleaf virus"]<-"GFLV"
RNAdat$short[RNAdat$Name=="Grapevine fleck virus"]<-"GFkV"
RNAdat$short[RNAdat$Name=="Arabis mosaic virus"]<-"ArMV"

RNAdat<-RNAdat[!is.na(RNAdat$short),]
rownames(RNAdat)<-RNAdat$short
RNAdat$Name<-RNAdat$Total<-RNAdat$short<-NULL
RNAdat<-data.frame(t(RNAdat),check.names=F)
RNAdat<-RNAdat[,sort(names(RNAdat))]
row.names(RNAdat)<-gsub("_FPKM","",row.names(RNAdat))

#Read results of RNA assignment using bracken
#Change names of Viruses: only viruses with short name corresponding to those present in RT or ELISA files are included in analysis.
Bdat<-fread(brackenfile,data.table=F)
Bdat$short<-NA
Bdat$short[Bdat$name=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
Bdat$short[Bdat$name=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
Bdat$short[Bdat$name=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
Bdat$short[Bdat$name=="Grapevine rupestris stem pitting-associated virus"]<-"RPSaV"
Bdat$short[Bdat$name=="Grapevine Pinot gris virus"]<-"GPGV"
Bdat$short[Bdat$name=="Grapevine virus A"]<-"GVA"
Bdat$short[Bdat$name=="Grapevine fanleaf virus"]<-"GFLV"
Bdat$short[Bdat$name=="Grapevine fleck virus"]<-"GFkV"
Bdat$short[Bdat$name=="Arabis mosaic virus"]<-"ArMV"

Bdat<-Bdat[!is.na(Bdat$short),]
row.names(Bdat)<-Bdat$short
Bdat$Total<-Bdat$name<-Bdat$taxonomy_id<-Bdat$short<-NULL
Bdat<-data.frame(t(Bdat),check.names=F)
Bdat<-Bdat[,sort(names(Bdat))]
row.names(Bdat)<-gsub("RPM_","",row.names(Bdat))
#Read results of smallRNAseq
smalldat<-fread(smallrnafile,data.table=F)
if(sum(names(smalldat)%in%"Name")>0) setnames(smalldat,"Name","Species")
smalldat$Taxid<-NULL
smalldat$short<-NA
smalldat$short[smalldat$Species=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
smalldat$short[smalldat$Species=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
smalldat$short[smalldat$Species=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
smalldat$short[smalldat$Species=="Grapevine rupestris stem pitting-associated virus"]<-"RPSaV"
smalldat$short[smalldat$Species=="Grapevine Pinot gris virus"]<-"GPGV"
smalldat$short[smalldat$Species=="Grapevine virus A"]<-"GVA"
smalldat$short[smalldat$Species=="Grapevine fanleaf virus"]<-"GFLV"
smalldat$short[smalldat$Species=="Grapevine fleck virus"]<-"GFkV"
smalldat$short[smalldat$Species=="Arabis mosaic virus"]<-"ArMV"
smalldat<-smalldat[!is.na(smalldat$short),]
row.names(smalldat)<-smalldat$short
smalldat$Species<-smalldat$short<-NULL
smalldat<-data.frame(t(smalldat),check.names=F,stringsAsFactors=F)

#Read results of smallRNAseq based on contig %

consmalldat<-fread(smallrnacontigfile,data.table=F)
if(length(grep("Number_",names(consmalldat)))>0)
{
consmalldat<-consmalldat[,grep("Number_",names(consmalldat),invert=T)]
names(consmalldat)<-gsub("Percent_","",names(consmalldat))
if(sum(names(consmalldat)%in%"Name")>0) setnames(consmalldat,"Name","Species")
consmalldat$Taxid<-NULL
} else {
setnames(consmalldat,c("Virus","U1","MN","MO"),c("Species","cort_gva_merged","maln_lr3","mont_neg_merged"))
}
consmalldat$short<-NA
consmalldat$short[consmalldat$Species=="Grapevine leafroll-associated virus 1"]<-"GLRaV-1"
consmalldat$short[consmalldat$Species=="Grapevine leafroll-associated virus 2"]<-"GLRaV-2"
consmalldat$short[consmalldat$Species=="Grapevine leafroll-associated virus 3"]<-"GLRaV-3"
consmalldat$short[consmalldat$Species=="Grapevine rupestris stem pitting-associated virus"]<-"RPSaV"
consmalldat$short[consmalldat$Species=="Grapevine Pinot gris virus"]<-"GPGV"
consmalldat$short[consmalldat$Species=="Grapevine virus A"]<-"GVA"
consmalldat$short[consmalldat$Species=="Grapevine fanleaf virus"]<-"GFLV"
consmalldat$short[consmalldat$Species=="Grapevine fleck virus"]<-"GFkV"
consmalldat$short[consmalldat$Species=="Arabis mosaic virus"]<-"ArMV"
#Horrible patch to adjust for the fact that in the contig smallRNA analysis we never find ArMV which is tested by ELISA, and thu should be counted as negative.
tomatch<-data.frame(short=c("GLRaV-1","GLRaV-2","GLRaV-3","RPSaV","GPGV","GVA","GFLV","GFkV","ArMV"))
pp<-merge(consmalldat,tomatch,by="short",all.y=T)                  
consmalldat<-pp
consmalldat<-consmalldat[!is.na(consmalldat$short),]
row.names(consmalldat)<-consmalldat$short
consmalldat$Species<-consmalldat$short<-NULL
consmalldat[is.na(consmalldat)]<-0
consmalldat<-data.frame(t(consmalldat),check.names=F,stringsAsFactors=F)


#Create the final data frame. We use Bdat (just in case) becasue is the one with the most analysed samples
finres<-data.frame(expand.grid(row.names(Bdat),colnames(Bdat)),stringsAsFactors=F)
setnames(finres,c("Var1","Var2"),c("Var","Virus"))
finres$RTPCR<-finres$ELISA<-NA
finres$aRNA<-finres$bRNA<-finres$smallRNA<-finres$consmallRNA<-0
finres$Var<-as.character(finres$Var)
finres$Virus<-as.character(finres$Virus)
for(aaa in 1:nrow(finres))
{
	if(length(grep(finres$Virus[aaa],names(PCRdat)))>0) finres$RTPCR[aaa]<- PCRdat[finres$Var[aaa],finres$Virus[aaa]]
	finres$aRNA[aaa]<- RNAdat[finres$Var[aaa],finres$Virus[aaa]]
	finres$bRNA[aaa]<- Bdat[finres$Var[aaa],finres$Virus[aaa]]
	if(length(grep(finres$Virus[aaa],names(smalldat)))>0) finres$smallRNA[aaa]<- smalldat[finres$Var[aaa],finres$Virus[aaa]]
	if(length(grep(finres$Virus[aaa],names(consmalldat)))>0) finres$consmallRNA[aaa]<- consmalldat[finres$Var[aaa],finres$Virus[aaa]]
    if(length(grep(finres$Virus[aaa],names(ELISAdat)))>0) finres$ELISA[aaa]<-ELISAdat[finres$Var[aaa],finres$Virus[aaa]]
}
#Correct for the fact that some cultivar were not analyzed with smallrna
finres$smallRNA[finres$Var=="ries_ara_merged"|finres$Var=="malb_lr1"|finres$Var=="dolc_rl31"|finres$Var=="cabf_neg"]<-NA
finres$labels_RT<-NA
finres$labels_RT[finres$RTPCR>0]<-TRUE
finres$labels_RT[finres$RTPCR==0]<-FALSE
finres$labels_ELISA<-NA
finres$labels_ELISA[finres$ELISA>0]<-TRUE
finres$labels_ELISA[finres$ELISA==0]<-FALSE
finres$labels<-NA
finres$labels[finres$labels_RT|finres$labels_ELISA]<-TRUE
finres$labels[!finres$labels_RT&is.na(finres$labels_ELISA)]<-FALSE
finres$labels[is.na(finres$labels_RT)&!(finres$labels_ELISA)]<-FALSE
finres$labels[!(finres$labels_RT)&!(finres$labels_ELISA)]<-FALSE

finres$aRNA<-round(finres$aRNA,2)
write.table(finres,outfile,sep="\t",quote=F,row.names=F)
#Remove Cabernet franc from ROC analyses, because cabernet didn't have RT-PCR or ELISA test.
finres<-finres[finres$Var!="cabf_neg",]

#finres<-finres[complete.cases(finres),]
mycor<-cor(finres[,4:8])
write.table(mycor,gsub(".txt","_cor.txt",outfile),sep="\t",quote=F)
if(includesmall) finres<-finres[!is.na(finres$smallRNA),]
#finres<-finres[!is.na(finres$both),]
mypdf<-ifelse(includesmall,"small_ROC.pdf","ROC.pdf")

pdf(paste(graphdir,mypdf,sep=""))
ROCit_obj_a_ELISA <- rocit(score=finres$aRNA,class=finres$labels_ELISA)
ROCit_obj_b_ELISA <- rocit(score=finres$bRNA,class=finres$labels_ELISA)
ROCit_obj_small_ELISA <- rocit(score=finres$smallRNA,class=finres$labels_ELISA)
ROCit_obj_consmall_ELISA <- rocit(score=finres$consmallRNA,class=finres$labels_ELISA)
if(includesmall) par(mfrow=c(3,4),oma = c(5,4,0,0) + 0.1,mar = c(0.8,0.8,0.8,0.8) + 0.1)
if(!includesmall) par(mfrow=c(3,2),oma = c(5,4,0,0) + 0.1,mar = c(0.8,0.8,0.8,0.8) + 0.1)
plot.rocit(ROCit_obj_a_ELISA,YIndex=F,axes=F,col=col.list,xtitle="RNAa",ytitle="ELISA",add.y.lab=TRUE,add.title=TRUE,legend=TRUE)
mtext("ELISA",side=2,cex=0.6,line=2)
plot.rocit(ROCit_obj_b_ELISA,YIndex=F,col=col.list,legend=TRUE)
if(includesmall) 
{
	plot.rocit(ROCit_obj_small_ELISA,YIndex=F,col=col.list,legend=TRUE)
	plot.rocit(ROCit_obj_consmall_ELISA,YIndex=F,col=col.list,legend=TRUE)
}	
	
ROCit_obj_a_RT <- rocit(score=finres$aRNA,class=finres$labels_RT)
ROCit_obj_b_RT <- rocit(score=finres$bRNA,class=finres$labels_RT)
ROCit_obj_small_RT <- rocit(score=finres$smallRNA,class=finres$labels_RT)
ROCit_obj_consmall_RT <- rocit(score=finres$consmallRNA,class=finres$labels_RT)

plot.rocit(ROCit_obj_a_RT,YIndex=F,add.y.lab=TRUE,col=col.list,legend=TRUE)
mtext("RT-PCR",side=2,cex=0.6,line=2)
plot.rocit(ROCit_obj_b_RT,YIndex=F,col=col.list,legend=TRUE)
if(includesmall) 
{
	plot.rocit(ROCit_obj_small_RT,YIndex=F,col=col.list,legend=TRUE)
	plot.rocit(ROCit_obj_consmall_RT,YIndex=F,col=col.list,legend=TRUE)
}
ROCit_obj_a_either <- rocit(score=finres$aRNA,class=finres$labels)
ROCit_obj_b_either <- rocit(score=finres$bRNA,class=finres$labels)
ROCit_obj_small_either <- rocit(score=finres$smallRNA,class=finres$labels)
ROCit_obj_consmall_either <- rocit(score=finres$consmallRNA,class=finres$labels)
plot.rocit(ROCit_obj_a_either,YIndex=F,add.x.lab=TRUE,add.y.lab=TRUE,col=col.list,legend=TRUE)
mtext("Both",side=2,cex=0.6,line=2)
mtext(bquote("total" ~  RNA^a),side=1,cex=0.6,line=2)
plot.rocit(ROCit_obj_b_either,YIndex=F,add.x.lab=TRUE,col=col.list,legend=TRUE)
mtext(bquote("total" ~ RNA^b),side=1,cex=0.6,line=2)
if(includesmall) 
{
	plot.rocit(ROCit_obj_small_either,YIndex=F,add.x.lab=TRUE,col=col.list,legend=TRUE)
	mtext(bquote("small" ~ RNA^a),side=1,cex=0.6,line=2)
	plot.rocit(ROCit_obj_consmall_either,YIndex=F,add.x.lab=TRUE,col=col.list,legend=TRUE)
	mtext(bquote("small" ~ RNA^c),side=1,cex=0.6,line=2)
}
dev.off()
#Never plotted the "consistent" results
finres$cons<-FALSE
finres$cons[finres$ELISA>0&finres$RTPCR>0]<-TRUE
finres$cons[finres$ELISA<1&finres$RTPCR<0.001]<-TRUE
finres<-finres[finres$cons,]

ROCit_obj_a_cons <- rocit(score=finres$aRNA,class=finres$labels)
ROCit_obj_b_cons <- rocit(score=finres$bRNA,class=finres$labels)
ROCit_obj_small_cons <- rocit(score=finres$smallRNA,class=finres$labels)
ROCit_obj_consmall_cons <- rocit(score=finres$consmallRNA,class=finres$labels)


#Summarize ROC results
opt_youden_a_ELISA<-ROCit_obj_a_ELISA$Cutoff[which.max(ROCit_obj_a_ELISA$TPR-ROCit_obj_a_ELISA$FPR)]
opt_TPR_a_ELISA<-ROCit_obj_a_ELISA$TPR[which.max(ROCit_obj_a_ELISA$TPR-ROCit_obj_a_ELISA$FPR)]
opt_FPR_a_ELISA<-ROCit_obj_a_ELISA$FPR[which.max(ROCit_obj_a_ELISA$TPR-ROCit_obj_a_ELISA$FPR)]
opt_youden_a_RT<-ROCit_obj_a_RT$Cutoff[which.max(ROCit_obj_a_RT$TPR-ROCit_obj_a_RT$FPR)]
opt_TPR_a_RT<-ROCit_obj_a_RT$TPR[which.max(ROCit_obj_a_RT$TPR-ROCit_obj_a_RT$FPR)]
opt_FPR_a_RT<-ROCit_obj_a_RT$FPR[which.max(ROCit_obj_a_RT$TPR-ROCit_obj_a_RT$FPR)]
opt_youden_a_either<-ROCit_obj_a_either$Cutoff[which.max(ROCit_obj_a_either$TPR-ROCit_obj_a_either$FPR)]
opt_TPR_a_either<-ROCit_obj_a_either$TPR[which.max(ROCit_obj_a_either$TPR-ROCit_obj_a_either$FPR)]
opt_FPR_a_either<-ROCit_obj_a_either$FPR[which.max(ROCit_obj_a_either$TPR-ROCit_obj_a_either$FPR)]
opt_youden_a_cons<-ROCit_obj_a_cons$Cutoff[which.max(ROCit_obj_a_cons$TPR-ROCit_obj_a_cons$FPR)]
opt_TPR_a_cons<-ROCit_obj_a_cons$TPR[which.max(ROCit_obj_a_cons$TPR-ROCit_obj_a_cons$FPR)]
opt_FPR_a_cons<-ROCit_obj_a_cons$FPR[which.max(ROCit_obj_a_cons$TPR-ROCit_obj_a_cons$FPR)]
opt_youden_b_ELISA<-ROCit_obj_b_ELISA$Cutoff[which.max(ROCit_obj_b_ELISA$TPR-ROCit_obj_b_ELISA$FPR)]
opt_TPR_b_ELISA<-ROCit_obj_b_ELISA$TPR[which.max(ROCit_obj_b_ELISA$TPR-ROCit_obj_b_ELISA$FPR)]
opt_FPR_b_ELISA<-ROCit_obj_b_ELISA$FPR[which.max(ROCit_obj_b_ELISA$TPR-ROCit_obj_b_ELISA$FPR)]
opt_youden_b_RT<-ROCit_obj_b_RT$Cutoff[which.max(ROCit_obj_b_RT$TPR-ROCit_obj_b_RT$FPR)]
opt_TPR_b_RT<-ROCit_obj_b_RT$TPR[which.max(ROCit_obj_b_RT$TPR-ROCit_obj_b_RT$FPR)]
opt_FPR_b_RT<-ROCit_obj_b_RT$FPR[which.max(ROCit_obj_b_RT$TPR-ROCit_obj_b_RT$FPR)]
opt_youden_b_either<-ROCit_obj_b_either$Cutoff[which.max(ROCit_obj_b_either$TPR-ROCit_obj_b_either$FPR)]
opt_TPR_b_either<-ROCit_obj_b_either$TPR[which.max(ROCit_obj_b_either$TPR-ROCit_obj_b_either$FPR)]
opt_FPR_b_either<-ROCit_obj_b_either$FPR[which.max(ROCit_obj_b_either$TPR-ROCit_obj_b_either$FPR)]
opt_youden_b_cons<-ROCit_obj_b_cons$Cutoff[which.max(ROCit_obj_b_cons$TPR-ROCit_obj_b_cons$FPR)]
opt_TPR_b_cons<-ROCit_obj_b_cons$TPR[which.max(ROCit_obj_b_cons$TPR-ROCit_obj_b_cons$FPR)]
opt_FPR_b_cons<-ROCit_obj_b_cons$FPR[which.max(ROCit_obj_b_cons$TPR-ROCit_obj_b_cons$FPR)]

opt_youden_small_ELISA<-ROCit_obj_small_ELISA$Cutoff[which.max(ROCit_obj_small_ELISA$TPR-ROCit_obj_small_ELISA$FPR)]
opt_TPR_small_ELISA<-ROCit_obj_small_ELISA$TPR[which.max(ROCit_obj_small_ELISA$TPR-ROCit_obj_small_ELISA$FPR)]
opt_FPR_small_ELISA<-ROCit_obj_small_ELISA$FPR[which.max(ROCit_obj_small_ELISA$TPR-ROCit_obj_small_ELISA$FPR)]
opt_youden_small_RT<-ROCit_obj_small_RT$Cutoff[which.max(ROCit_obj_small_RT$TPR-ROCit_obj_small_RT$FPR)]
opt_TPR_small_RT<-ROCit_obj_small_RT$TPR[which.max(ROCit_obj_small_RT$TPR-ROCit_obj_small_RT$FPR)]
opt_FPR_small_RT<-ROCit_obj_small_RT$FPR[which.max(ROCit_obj_small_RT$TPR-ROCit_obj_small_RT$FPR)]
opt_youden_small_either<-ROCit_obj_small_either$Cutoff[which.max(ROCit_obj_small_either$TPR-ROCit_obj_small_either$FPR)]
opt_TPR_small_either<-ROCit_obj_small_either$TPR[which.max(ROCit_obj_small_either$TPR-ROCit_obj_small_either$FPR)]
opt_FPR_small_either<-ROCit_obj_small_either$FPR[which.max(ROCit_obj_small_either$TPR-ROCit_obj_small_either$FPR)]
opt_youden_small_cons<-ROCit_obj_small_cons$Cutoff[which.max(ROCit_obj_small_cons$TPR-ROCit_obj_small_cons$FPR)]
opt_TPR_small_cons<-ROCit_obj_small_cons$TPR[which.max(ROCit_obj_small_cons$TPR-ROCit_obj_small_cons$FPR)]
opt_FPR_small_cons<-ROCit_obj_small_cons$FPR[which.max(ROCit_obj_small_cons$TPR-ROCit_obj_small_cons$FPR)]

opt_youden_consmall_ELISA<-ROCit_obj_consmall_ELISA$Cutoff[which.max(ROCit_obj_consmall_ELISA$TPR-ROCit_obj_consmall_ELISA$FPR)]
opt_TPR_consmall_ELISA<-ROCit_obj_consmall_ELISA$TPR[which.max(ROCit_obj_consmall_ELISA$TPR-ROCit_obj_consmall_ELISA$FPR)]
opt_FPR_consmall_ELISA<-ROCit_obj_consmall_ELISA$FPR[which.max(ROCit_obj_consmall_ELISA$TPR-ROCit_obj_consmall_ELISA$FPR)]
opt_youden_consmall_RT<-ROCit_obj_consmall_RT$Cutoff[which.max(ROCit_obj_consmall_RT$TPR-ROCit_obj_consmall_RT$FPR)]
opt_TPR_consmall_RT<-ROCit_obj_consmall_RT$TPR[which.max(ROCit_obj_consmall_RT$TPR-ROCit_obj_consmall_RT$FPR)]
opt_FPR_consmall_RT<-ROCit_obj_consmall_RT$FPR[which.max(ROCit_obj_consmall_RT$TPR-ROCit_obj_consmall_RT$FPR)]
opt_youden_consmall_either<-ROCit_obj_consmall_either$Cutoff[which.max(ROCit_obj_consmall_either$TPR-ROCit_obj_consmall_either$FPR)]
opt_TPR_consmall_either<-ROCit_obj_consmall_either$TPR[which.max(ROCit_obj_consmall_either$TPR-ROCit_obj_consmall_either$FPR)]
opt_FPR_consmall_either<-ROCit_obj_consmall_either$FPR[which.max(ROCit_obj_consmall_either$TPR-ROCit_obj_consmall_either$FPR)]
opt_youden_consmall_cons<-ROCit_obj_consmall_cons$Cutoff[which.max(ROCit_obj_consmall_cons$TPR-ROCit_obj_consmall_cons$FPR)]
opt_TPR_consmall_cons<-ROCit_obj_consmall_cons$TPR[which.max(ROCit_obj_consmall_cons$TPR-ROCit_obj_consmall_cons$FPR)]
opt_FPR_consmall_cons<-ROCit_obj_consmall_cons$FPR[which.max(ROCit_obj_consmall_cons$TPR-ROCit_obj_consmall_cons$FPR)]

rres<-matrix(c(opt_youden_a_ELISA,ROCit_obj_a_ELISA$AUC,opt_TPR_a_ELISA,opt_FPR_a_ELISA,length(ROCit_obj_a_ELISA$pos_D),length(ROCit_obj_a_ELISA$neg_D),
               opt_youden_b_ELISA,ROCit_obj_b_ELISA$AUC, opt_TPR_b_ELISA,opt_FPR_b_ELISA,length(ROCit_obj_b_ELISA$pos_D),length(ROCit_obj_b_ELISA$neg_D),
               opt_youden_consmall_ELISA,ROCit_obj_consmall_ELISA$AUC,opt_TPR_consmall_ELISA,opt_FPR_consmall_ELISA,length(ROCit_obj_consmall_ELISA$pos_D),length(ROCit_obj_consmall_ELISA$neg_D),
               opt_youden_small_ELISA,ROCit_obj_small_ELISA$AUC,opt_TPR_small_ELISA,opt_FPR_small_ELISA,length(ROCit_obj_small_ELISA$pos_D),length(ROCit_obj_small_ELISA$neg_D),
               opt_youden_a_RT,ROCit_obj_a_RT$AUC,opt_TPR_a_RT,opt_FPR_a_RT,length(ROCit_obj_a_RT$pos_D),length(ROCit_obj_a_RT$neg_D),
               opt_youden_b_RT,ROCit_obj_b_RT$AUC,opt_TPR_b_RT,opt_FPR_b_RT,length(ROCit_obj_b_RT$pos_D),length(ROCit_obj_b_RT$neg_D),
               opt_youden_consmall_RT,ROCit_obj_consmall_RT$AUC,opt_TPR_consmall_RT,opt_FPR_consmall_RT,length(ROCit_obj_consmall_RT$pos_D),length(ROCit_obj_consmall_RT$neg_D),
               opt_youden_small_RT,ROCit_obj_small_RT$AUC,opt_TPR_small_RT,opt_FPR_small_RT,length(ROCit_obj_small_RT$pos_D),length(ROCit_obj_small_RT$neg_D),
               opt_youden_a_either,ROCit_obj_a_either$AUC,opt_TPR_a_either,opt_FPR_a_either,length(ROCit_obj_a_either$pos_D),length(ROCit_obj_a_either$neg_D),
               opt_youden_b_either,ROCit_obj_b_either$AUC,opt_TPR_b_either,opt_FPR_b_either,length(ROCit_obj_b_either$pos_D),length(ROCit_obj_b_either$neg_D),
               opt_youden_consmall_either,ROCit_obj_consmall_either$AUC,opt_TPR_consmall_either,opt_FPR_consmall_either,length(ROCit_obj_consmall_either$pos_D),length(ROCit_obj_consmall_either$neg_D),
               opt_youden_small_either,ROCit_obj_small_either$AUC,opt_TPR_small_either,opt_FPR_small_either,length(ROCit_obj_small_either$pos_D),length(ROCit_obj_small_either$neg_D),
			   opt_youden_a_cons,ROCit_obj_a_cons$AUC,opt_TPR_a_cons,opt_FPR_a_cons,length(ROCit_obj_a_cons$pos_D),length(ROCit_obj_a_cons$neg_D),
               opt_youden_b_cons,ROCit_obj_b_cons$AUC,opt_TPR_b_cons,opt_FPR_b_cons,length(ROCit_obj_b_cons$pos_D),length(ROCit_obj_b_cons$neg_D),
               opt_youden_consmall_cons,ROCit_obj_consmall_cons$AUC,opt_TPR_consmall_cons,opt_FPR_consmall_cons,length(ROCit_obj_consmall_cons$pos_D),length(ROCit_obj_consmall_cons$neg_D),
               opt_youden_small_cons,ROCit_obj_small_cons$AUC,opt_TPR_small_cons,opt_FPR_small_cons,length(ROCit_obj_small_cons$pos_D),length(ROCit_obj_small_cons$neg_D)),
      ncol=6,byrow=T,dimnames=list(c("aRNA_ELISA","bRNA_ELISA","contig_small_ELISA","small_ELISA","aRNA_RT","bRNA_RT","contig_small_RT","small_RT","aRNA_either","bRNA_either","contig_small_either","small_either","aRNA_consistent","bRNA_consistent","contig_small_consistent","small_consistent"),c("Youden","AUC","TPR","FPR","Positives","Negatives")))
mypattern<-ifelse(includesmall,"small_RNA_RTPCR_ROC","RNA_RTPCR_ROC")
write.table(rres,gsub("RNA_smallRNA_RTPCR_gold_std",mypattern,outfile),sep="\t",quote=F)


}
RTPCR_RNAseq(brackenfile=brackenfile,pcrfile=pcrfile,elisafile=elisafile,rnaseqfile=rnaseqfile,smallrnafile=smallrnafile,smallrnacontigfile=smallrnacontigfile,outfile=outfile,graphdir=graphdir,includesmall=includesmall)


#https://blog.revolutionanalytics.com/2016/08/roc-curves-in-two-lines-of-code.html
