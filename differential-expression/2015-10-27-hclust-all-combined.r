library("edgeR")


# =============================================
# Normalize gonad tissue
# =============================================
data <- read.csv("../count-extraction/2015-10-29-MF-ALL-TISSUES-COUNTS-COMBINED.expr", header=T, row.names=1, stringsAsFactors=FALSE)
dim(data)
names(data)
conditions <- factor(c(rep("M", 16),rep("F", 15)))
conditions
expr <- DGEList(counts=data, group=conditions)
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr)
# Compare Male - Female expression
# . Note
# that the first group listed in the pair is the baseline for the comparison—so if exactTest
# the pair is c("A","B") then the comparison is B - A, so genes with positive
# log-fold change are up-regulated in group B compared with group A (and vice
# versa for genes with negative log-fold change
de <- exactTest(expr, c("M","F"))
#plotMDS(expr)
# Save all tags in this table. nrow(expr) makes sure that all genes are
# included. We need to select the $table, because topTags also contains other
# elements like comparison, test
detag_table <- topTags(de, n=nrow(expr))$table
write.table(detag_table, file="2015-10-29-MF-ALL-TISSUES-detags-edgeR.csv",quote=F, sep=",")

# Export the cpm values for all genes
# prior.count=0.125 is the default for exactTest() and is needed here to make values comparable.
expr_norm <- cpm(expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125)
write.table(expr_norm, file="2015-10-29-MF-ALL-TISSUES-CPM-NORM-log2-edgeR.csv",quote=F, sep=",")

# RPKM calculation
# Gene length based on combined exonic length
#IMPORTANT: NEED TO REMOVE HEADER V1!
gene_length <- read.table("2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length")
expressed_genes <- rownames(expr)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
dim(gene_length)
# Reorder gene length vector to match order in expr object
gene_length_ordered <- gene_length[match(rownames(expr),gene_length$V1),]
# Check if reorder has worked ( should return TRUE)
all(gene_length_ordered$V1 == rownames(expr))
expr_norm <- rpkm(expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125, gene.length=gene_length_ordered$V2)
write.table(expr_norm, file="2015-10-29-MF-ALL-TISSUES-RPKM-NORM-log2-edgeR.csv",quote=F, sep=",")

# =============================================
# Hiearchical clusters of all samples
# =============================================
data <- read.csv("2015-10-29-MF-ALL-TISSUES-CPM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
plot(hclust(dist(t(data))))
png(file="2015-10-29-MF-ALL-TISSUES-CPM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()
# =============================================
# Hiearchical clusters of all samples
# =============================================
data <- read.csv("2015-10-29-MF-ALL-TISSUES-RPKM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
plot(hclust(dist(t(data))))
png(file="2015-10-29-MF-ALL-TISSUES-RPKM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()