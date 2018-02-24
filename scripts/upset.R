#libraries
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2)
{
  stop("Usage: Rscript chipcounts.R up_distribution outputname", call.=FALSE)
}

suppressPackageStartupMessages({
library(UpSetR)
  })

inputfile = args[1]
outputname = args[2]

df = read.table(inputfile,sep="\t",header = T,row.names = 1)

pdf(file=paste0(outputname,".pdf"))
upset(df, sets = dput(as.character(colnames(df))), sets.bar.color = "#56B4E9",
order.by = "freq", empty.intersections = "on")
dev.off()

