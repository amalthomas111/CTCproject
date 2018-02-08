args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3)
{
  stop("Usage: Rscript runexomeCNV.R <tumor.sample_interval_summary> <normal.sample_interval_summary> <outputname.", call.=FALSE)
}
suppressPackageStartupMessages({
  library(ExomeCNV)
})

# output files:
# exomeCNVoutput/outputname.cnv.png
# exomeCNVoutput/outputname.cnv.txt
# exomeCNVoutput/outputname.exon.lrr.txt
# exomeCNVoutput/outputname.segment.copynumber.txt
# exomeCNVoutput/outputname.segment.lrr.txt

mainDir=getwd()
subDir="exomeCNVoutput/"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)

chr.list = paste0("chr",c(1:22,"X","Y"))

tumor = read.coverage.gatk(args[1])
normal = read.coverage.gatk(args[2])
outputname = args[3]

tumor$probe = gsub("chr","",tumor$probe)
normal$probe = gsub("chr","", normal$probe)
tumor$chr = gsub("chrchr","chr",tumor$chr)
normal$chr = gsub("chrchr","chr",normal$chr)

demo.logR = calculate.logR(normal, tumor,normal.chrs = chr.list)
demo.eCNV = c()
for (i in 1:length(chr.list)) {
  idx = (normal$chr == chr.list[i])
  ecnv = classify.eCNV(normal=normal[idx,], tumor=tumor[idx,], logR=demo.logR[idx], min.spec=0.9999, min.sens=0.9999, option="spec", c=0.5, l=70)
  demo.eCNV = rbind(demo.eCNV, ecnv)
}

demo.cnv = multi.CNV.analyze(normal, tumor, logR=demo.logR, all.cnv.ls=list(demo.eCNV), coverage.cutoff=5, min.spec=0.99, min.sens=0.99, option="auc", c=0.5)



myplot = function (all.ecnv, pch = "*", lim.quantile = 0.99, style = "idx", 
                   bg.cnv = NULL, line.plot = FALSE) {
  if (is.null(bg.cnv)) {
    lim.logR = quantile(all.ecnv$logR[all.ecnv$logR != Inf & 
                                        all.ecnv$logR != -Inf], c(1 - lim.quantile, lim.quantile), 
                        na.rm = TRUE)
  }
  else {
    lim.logR = quantile(bg.cnv$logR[bg.cnv$logR != Inf & 
                                      bg.cnv$logR != -Inf], c(1 - lim.quantile, lim.quantile), 
                        na.rm = TRUE)
  }
  if (lim.logR[1] == -Inf) 
    lim.logR[1] = min(all.ecnv$logR[all.ecnv$logR != -Inf], 
                      na.rm = TRUE)
  if (lim.logR[2] == Inf) 
    lim.logR[2] = max(all.ecnv$logR[all.ecnv$logR != Inf], 
                      na.rm = TRUE)
  max.cn = max(all.ecnv$copy.number, na.rm = TRUE)
  reds = rainbow(2^(max.cn - 2), start = 3/4, end = 0)[2^(max.cn:3 - 
                                                            2)]
  colors = c("green4", "gold", reds)
  chr.list = unique(as.character(all.ecnv$chr))
  #dimx = floor(sqrt(length(chr.list)))
  #dimy = ceiling(length(chr.list)/dimx)
  #par(mfrow = c(dimx, dimy))
  
  for (chr in chr.list) {
    ExomeCNV:::do.plot.one.eCNV(all.ecnv[all.ecnv$chr == chr, ], colors = colors, 
                                ylim = lim.logR, main = chr, pch = pch, style = style, 
                                bg.cnv = bg.cnv[bg.cnv$chr == chr, ], line.plot = line.plot)
  }
  
}

pdf(file=paste0("exomeCNVoutput/",outputname,".pdf"))
myplot(demo.cnv, lim.quantile=0.99, style="bp", bg.cnv=demo.eCNV, line.plot=T)
dev.off()

write.output(demo.eCNV, demo.cnv, paste0("exomeCNVoutput/",outputname))