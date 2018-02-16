#!/usr/bin/Rscript
#Author: A.Thomas
#Script to get promoter location(+/-2.5Kb from TSS) for hg19 genes
#using ENSEMBL biomart

suppressPackageStartupMessages({
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
})

width=2500
name=as.character(width/1000)
outputname=paste0("hg19_promoters_",name,"kb")
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
#GRCh=37)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                     GRCh=37)
#head(listFilters(ensembl))
#head(listAttributes(ensembl))
#head(listAttributes(ensembl))

genes <- getBM(attributes=c('ensembl_gene_id','chromosome_name',
                            'start_position','end_position','strand',
                            'transcript_start'),  mart = ensembl)

head(genes)
dim(genes)
genes.uniq = genes[!duplicated(genes$ensembl_gene_id),]
dim(genes.uniq)
#keep only those genes present in our GTF file
gtf_genes = read.table('hg19.genenames')
genes.uniq.selected = genes.uniq[genes.uniq$ensembl_gene_id  %in% 
                                 unlist(gtf_genes),]
dim(genes.uniq.selected)

gr1_TSS = with(genes.uniq.selected,
               GRanges(seqnames = paste0("chr",chromosome_name),
                    ranges = IRanges(start=transcript_start,
                    end=transcript_start),strand=strand,id=ensembl_gene_id))
head(gr1_TSS)
gr1_promoter = flank(gr1_TSS,width,both=TRUE)
start(gr1_promoter) = ifelse(start(gr1_promoter) < 1,1,start(gr1_promoter))
head(gr1_promoter)
print(outputname)
export(gr1_promoter,paste0(outputname,".bed"),format="bed")
system(paste0("sort -k1,1 -k2,2n ",outputname,".bed > ",
              outputname,".sorted.bed"))
system(paste0("cut -f1-3 ",outputname,".bed| sort -k1,1 -k2,2n |
              bedtools merge -i - > ",outputname,".sorted.merged.bed"))
system(paste0("awk -v OFS='\\t' '{print \"region_\"NR,$1,$2,$3,\".\"}' ",
              outputname,".sorted.merged.bed > ",outputname,".merged.saf"))
