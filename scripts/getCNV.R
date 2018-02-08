args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript getCNV.R <VARSCAN.copynumber.called> <outputname>", call.=FALSE)
}
suppressPackageStartupMessages({
  library(DNAcopy)
})

cnv.input = args[1]
outputname = args[2]
cn = read.table(cnv.input,header = TRUE,sep = "\t")
chroms = c(paste0("chr",c(1:22,"M","X","Y")))
cn  = cn[cn$chrom %in% chroms,]
head(cn)

mainDir=getwd()
subDir="plots"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="savedobj"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
subDir="CNVsegments"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

CNA.object = CNA(genomdat = cn$adjusted_log_ratio,chrom = cn$chrom, maploc = cn$chr_start, data.type = "logratio")
smoothed.CNA.object <- smooth.CNA(CNA.object)
sdundo.CNA.object <- segment(smoothed.CNA.object, undo.splits="sdundo",undo.SD=3,verbose=1)
psegs = segments.p(sdundo.CNA.object)
write.table(psegs[,2:10], file=paste0("CNVsegments/",outputname,"_CNV.cb.tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
save(sdundo.CNA.object,file = paste0("savedobj/",outputname,"_cnaobj.Robj"))

pdf(file = paste0("plots/",outputname,"_CNV.pdf"))
DNAcopy:::plot.DNAcopy(sdundo.CNA.object,plot.type="c", cbys.layout = c(1,1), cbys.nchrom = 1)
dev.off()
