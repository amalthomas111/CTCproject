#!/usr/bin/Rscript
# Author: A. Thomas
# Rscript to create footprinting plots 
# Input required: coverage R objects for the bams in savedobj folder
# To Run: Rscript script <motifname> <motif index shown by MotifDb list> 
# <bam1 object name> <bam2 object name> <de regions between bam1 and 
# bam2 (optional)>
# E.g Rscript plot_motiffootprint.R MYC 3 brain_fullcoverage.Robj 
# nonbrain_fullcoverage.Robj  brain_nonbrain_DE.bed
# If motif index is not known just run Rscript script <motifname>. 
# Then from the search results of the motif by motifdb package choose the index.
# Providing zero as the index script tries to choose vertebrate motif
# from JASPAR, Jolma, HOCOMOCO or SwissRegulon
# DE bed file needs to have a header with atleast three columns
# chr\tstart\tend.
# To creat coverage object for the bamfiles run create_bamcoverage.R
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0)
{
  stop(paste("Usage:Rscript script <motifname> <motifindex(0=auto,manual no.)>"
             ," <bam1 coverage object> <bam2 coverage object> ",
             "<de_region btwn bam1 and bam2(optional)>"), call.=FALSE)
}

if(!file.exists(paste0("plotFootprints_edit.R"))){
 stop("plotFootprints_edit.R not found")
}
suppressPackageStartupMessages({
library(ATACseqQC)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
library(motifStack)
library(BiocParallel)
library(GenomicAlignments)
library(ChIPpeakAnno)
source("plotFootprints_edit.R")
})

motif = args[1]
q = query(MotifDb,motif)
q = as.list(q)
cat("\n\n---------------------------------\n\n")
names(q)
cat("\n\n---------------------------------\n\n")

if (length(args)<4)
{
  stop(paste("Usage:Rscript script <motifname> <motifindex(0=auto,manual no.)>"
             ," <bam1 coverage object> <bam2 coverage object> ",
             "<de_region btwn bam1 and bam2(optional)>"), call.=FALSE)
}
cat("\nMotif Selection\n")
if(args[2]==0){
  index = 0
  if(any(grepl("Vertebrata-JASPAR_CORE",names(q)))){
    index = grep("Vertebrata-JASPAR_CORE",names(q))[1]
    cat("\nVertebrata-JASPAR_CORE found\n")
  }else if(any(grepl("Hsapiens-JASPAR_CORE",names(q)))){
    index = grep("Hsapiens-JASPAR_CORE",names(q))[1]
    cat("\nHsapiens-JASPAR_CORE found\n")
  }else if(any(grepl("Hsapiens-jaspar2016",names(q)))){
    index = grep("Hsapiens-jaspar2016",names(q))[1]
    cat("\nHsapiens-jaspar2016 found\n")
  }else if(any(grepl("Hsapiens-JASPAR_2014",names(q)))){
    index = grep("JASPAR_2014",names(q))[1]
    cat("\n JASPAR_2014 found\n")
  }else if(any(grepl("Hsapiens-jolma2013",names(q)))){
    index = grep("Hsapiens-jolma2013",names(q))[1]
    cat("\nHsapiens-jolma2013 found\n")
  }else if(any(grepl("Hsapiens-HOCOMOCOv10",names(q)))){
    index = grep("Hsapiens-HOCOMOCOv10",names(q))[1]
    cat("\nHsapiens-HOCOMOCOv10 found\n")
  }else if(any(grepl("Hsapiens-HOCOMOCOv10",names(q)))){
    index = grep("Hsapiens-HOCOMOCOv10",names(q))[1]
    cat("\nHsapiens-HOCOMOCOv10 found\n")
  }else if(any(grepl("Hsapiens-SwissRegulon",names(q)))){
    index = grep("Hsapiens-SwissRegulon",names(q))[1]
    cat("\nHsapiens-SwissRegulon found\n")
  }else
    stop("\nNo Vertebrata core/human motif found!!\n",call.=FALSE)
  }else{
  index=args[2]
  cat("\nManual index=",index," selected\n")
}

cat("\nSelected Motif\n",names(q)[as.integer(index)])

motifname=""
if(basename(names(q)[as.integer(index)]) == "Homer"){
 motifname=strsplit(names(q)[as.integer(index)],"/")[[1]][2]
}else{
 motifname=names(q)[as.integer(index)]
}

cat("\n\nStart:")
print(Sys.time())
cat("\n")

genome <- Hsapiens
seqlev = paste0("chr",c(c(1:22),"X","Y"))
upstream = 100
downstream = 100
name1 = gsub("_fullcoverage.Robj","",args[3])
name2 = gsub("_fullcoverage.Robj","",args[4])
outputname = paste0(name1,"_",name2)
cat("\nOutputname:",outputname,"\n")

#### Reading first object #####
if(file.exists(paste0("savedobj/",args[3]))){
	load(paste0("savedobj/",args[3]))
	cat("\nbam1 coverage object found\n")
}else {
 stop(paste("bam1 coverage object not found in savedobj folder",
 "run create_bamcoverage.R for the bam file"))
}
cvgSum.lib1 = bamIn[names(bamIn) %in% seqlev]

#### Reading second object #####
if(file.exists(paste0("savedobj/",args[4]))){
	load(paste0("savedobj/",args[4]))
	cat("\nbam2 coverage object found\n")
}else {
 stop(paste("bam2 coverage object not found in savedobj folder",
 "run create_bamcoverage.R for the bam file"))
}
cvgSum.lib2 = bamIn[names(bamIn) %in% seqlev]


mainDir=getwd()
subDir="plots"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="savedobj"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="savedobj/motif"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="DE_Stat"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)


pfm = q[[as.integer(index)]]
wid <- ncol(pfm)

if(file.exists(paste0("savedobj/motif/",motifname,".Robj"))){
	load(paste0("savedobj/motif/",motifname,".Robj"))
	cat("\nMotif object already found\n")
}else {
  cat("\nCalculating potential binding sites of the motif\n")
	mt <- matchPWM(motifStack::pfm2pwm(pfm), genome, min.score = "95%", 
                   with.score=TRUE,
                   exclude=names(genome)[!names(genome) %in% seqlev])
	mt <- mt[seqnames(mt) %in% seqlev]
	save(mt,file= paste0("savedobj/motif/",motifname,".Robj"))
	lines = paste0(motif,"\t",names(q)[as.integer(index)],"\t",length(mt),"\n")
	write(lines,file = "savedobj/motif/no_of_motifs_auto.txt",append = TRUE)
}

mt.s <- split(mt, seqnames(mt))
seqlev <- intersect(unique(names(cvgSum.lib1),names(cvgSum.lib2)), names(mt.s))
mt.s <- mt.s[seqlev]

cvgSum.lib1 <- cvgSum.lib1[seqlev]
cvgSum.lib2 <- cvgSum.lib2[seqlev]
cvglist = list()
cvglist[["+"]] = cvgSum.lib1
cvglist[["-"]] = cvgSum.lib2
cvgSum <- cvglist[["+"]] + cvglist[["-"]]
mt.v <- Views(cvgSum, mt.s)
mt.s <- mt.s[viewSums(mt.v)>0]
mt <- unlist(mt.s)

sigs <- featureAlignedSignal(cvglists=cvglist,
                             feature.gr=reCenterPeaks(mt, width=1),
                             upstream=upstream+floor(wid/2),
                             downstream=downstream+ceiling(wid/2),
                             n.tile=upstream+downstream+wid)

s=colMeans(do.call(cbind, sigs))

pdf(file = paste0("plots/",outputname,"_",motifname,"_genomewide.pdf"))
plotFootprints(s,Mlen=wid, motif=new("pfm", mat = as.matrix(pfm), name= motif),
               lab1 = name1, lab2 = name2,xlab= names(q)[as.integer(index)])
dev.off()

print(Sys.time())
cat("\nGenome wide Done:\n")

if (length(args)!=5)
{
  stop("No DE region given. Exiting!!\n", call.=FALSE)
}
defile = args[5]
path="DE_regions/"
df = read.table(paste0(path,defile),sep = "\t",header = T)
head(df)

deregion = GRanges(seqnames=df$chr,
                   ranges=IRanges(start=df$start+1,end=df$end),
                   strand=rep("*",length(df$chr)))
#head(deregion)
length(deregion)


mt.de = mt[mt %over% deregion]
chroms = dput(as.character(unique(seqnames(deregion))))
seqlev.de = chroms
cat("No of motifs in de regions:",length(mt.de))
lines = paste0(outputname,"\t",names(q)[as.integer(index)],"\t",
               args[5],"\t",length(mt.de))
write(lines,file = "DE_Stat/noofmotifs_in_DE.tsv",append = TRUE)
if(length(mt.de) >= 20){
  #seqinfo(mt.de) <- Seqinfo(seqlev.de, seqlengths = seqlengths(mt.de))
  mt.de.s <- split(mt.de, seqnames(mt.de))
  seqlev <- intersect(names(cvgSum), names(mt.de.s))
  mt.de.s <- mt.de.s[seqlev]
  mt.de.v <- Views(cvgSum, mt.de.s)
  mt.de.s <- mt.de.s[viewSums(mt.de.v)>0]
  mt.de <- unlist(mt.de.s)
  sigs.de <- featureAlignedSignal(cvglists=cvglist,
                                  feature.gr=reCenterPeaks(mt.de, width=1),
                                  upstream=upstream+floor(wid/2),
                                  downstream=downstream+ceiling(wid/2),
                                  n.tile=upstream+downstream+wid)

  s.de=colMeans(do.call(cbind, sigs.de))
  pdf(file = paste0("plots/",outputname,"_",motifname,"_deregion.pdf"))
  plotFootprints(s.de,Mlen=wid, motif=new("pfm", mat = as.matrix(pfm),
                                          name= motif),
                 lab1 = name1, lab2 = name2, xlab= names(q)[as.integer(index)])
  dev.off()
}else {
  cat("\n # motif found in DE region is less than 20\n")
}

print(Sys.time())
cat("\nMotif DE region Done:\n")
