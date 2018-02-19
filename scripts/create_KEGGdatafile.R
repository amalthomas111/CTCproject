#Author: A.Thomas
#Script to create Kegg pathway data files
suppressPackageStartupMessages({
  library(dplyr)
  library(goseq)
  library(KEGGREST)
  library(biomaRt)
})


human <- useEnsembl('ensembl', dataset = 'hsapiens_gene_ensembl', GRCh=37)
attributes <- c('ensembl_gene_id',
                'entrezgene')

ensembl.entrez <- getBM(attributes,
                                  values="*", 
                                  mart = human,
                                  bmHeader=FALSE)
ensembl.entrez = ensembl.entrez[complete.cases(ensembl.entrez),]
ensembl2entrez <- as.character(ensembl.entrez$entrezgene)
names(ensembl2entrez) <- ensembl.entrez$ensembl_gene_id
ensembl2entrez <- as.list(ensembl2entrez)

get_kegg <- function(l)
{
  return( keggLink("pathway", l))
}

entrez.kegg <- paste('hsa', unique(ensembl2entrez), sep=":")
entrez.kegg.splitted <- split( entrez.kegg, 
                              ceiling(seq_along( entrez.kegg)/20))

entrez.kegg.pathways <- lapply(entrez.kegg.splitted, get_kegg) 
entrez.kegg.pathways.red <- Reduce(c, entrez.kegg.pathways)
grepKEGG <- function(geneid, pathway_map){
  if (!is.na(geneid))
  {
    return (unique(unlist(pathway_map[names(pathway_map) == 
                          paste('hsa', geneid, sep=':')])))
  }
}
gene2cat <- lapply(ensembl2entrez,grepKEGG, entrez.kegg.pathways.red)
dir.create(file.path(getwd(),"data"), showWarnings = FALSE)
save(entrez.kegg.pathways,file=paste0("data/Kegg.pathway.red_",
                                      Sys.Date(),".Robj"))
save(gene2cat,file=paste0("data/Kegg.gene2cat_",Sys.Date(),".Robj"))
