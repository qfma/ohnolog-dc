# Information taken from https://www.biostars.org/p/83901/
# First, import the GTF-file that you have also used as input for htseq-count
# library(GenomicFeatures)
# txdb <- makeTranscriptDbFromGFF("yourFile.gtf",format="gtf")
# # then collect the exons per gene id
# exons.list.per.gene <- exonsBy(txdb,by="gene")
# # then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
# exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})

library(biomaRt)
library(GenomicFeatures)
txdb <- makeTranscriptDbFromBiomart(dataset="ggallus_gene_ensembl",
									biomart="ENSEMBL_MART_ENSEMBL",
									host='feb2014.archive.ensembl.org') # Use Ensembl75
txdb
exons.list.per.gene <- exonsBy(txdb,by="gene")
exonic.gene.sizes <- lapply(exons.list.per.gene,function(x){sum(width(reduce(x)))})
# Convert list into dataframe and reverse rows/ columns
exonic_length <- t(as.data.frame(exonic.gene.sizes))
#IMPORTANT: NEED TO REMOVE HEADER V1!
write.table(exonic_length, file="2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length",quote=F)