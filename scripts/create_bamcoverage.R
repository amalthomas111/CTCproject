#!/usr/bin/Rscript
# Rscript to create coverage object for bam file
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1)
{
  stop("Usage: Rscript script bamfile", call.=FALSE)
}

suppressPackageStartupMessages({
library(GenomicAlignments)
})


chroms = paste0("chr", c(1:22, "X","Y"))
seqlev = chroms


filename = args[1]
name = gsub(".sorted.bam","",filename)
bamIn <- mapply(readGAlignments,filename,filename, SIMPLIFY = FALSE)
bamIn <- lapply(bamIn, as, Class = "GRanges")
if(class(bamIn)!="GRangesList") bamIn <- GRangesList(bamIn)
bamIn <- unlist(bamIn)
bamIn <- promoters(bamIn, upstream=0, downstream=1)
#save(bamIn, file = paste0("savedobj/",name,".Robj"))

#normalize for seq depth. Change the scale if required
scale = 10^8/length(bamIn)
cat("\n",name,":\n"scale:",scale,"\n")
bamIn = coverage(bamIn, weight = scale)
subDir="savedobj"
dir.create(file.path(mainDir,subDir), showWarnings = FALSE)
save(bamIn, file = paste0("savedobj/",name,"_fullcoverage.Robj"))
