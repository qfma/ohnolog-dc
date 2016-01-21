unlink(".RData")
#REMOVE everything from the workspace
rm(list=ls(all=TRUE))
library("edgeR")

# =============================================
# Normalize gonad tissue
# =============================================
gonad_data <- read.csv("../count-extraction/2015-10-27-MF-GONAD-COUNTS-COMBINED.expr", header=T, row.names=1, stringsAsFactors=FALSE)
dim(gonad_data)
names(gonad_data)
# cor(data)
class(gonad_data[,1])
class(gonad_data[,2])
class(gonad_data[,3])
class(gonad_data[,4])
class(gonad_data[,5])
class(gonad_data[,6])
class(gonad_data[,7])
class(gonad_data[,8])
gonad_conditions <- factor(c("M","M","M","M","F","F","F","F"))
gonad_expr <- DGEList(counts=gonad_data, group=gonad_conditions)
#Default is TMM
# manual: The default method for computing these scale factors uses a 
# trimmed mean of Mvalues (TMM) between each pair of samples [21]
gonad_expr <- calcNormFactors(gonad_expr)
# Calculate common and tagwise dispersion at the same time (recommended)
gonad_expr <- estimateDisp(gonad_expr)
# Print expression samples
gonad_expr$samples
# Compare Male - Female expression
# . Note
# that the first group listed in the pair is the baseline for the comparison—so if exactTest
# the pair is c("A","B") then the comparison is B - A, so genes with positive
# log-fold change are up-regulated in group B compared with group A (and vice
# versa for genes with negative log-fold change
gonad_de <- exactTest(gonad_expr, pair=c("M","F"))
#plotMDS(expr)
# Save all tags in this table. nrow(expr) makes sure that all genes are
# included. We need to select the $table, because topTags also contains other
# elements like comparison, test
gonad_detag_table <- topTags(gonad_de, n=nrow(gonad_expr))$table
write.table(gonad_detag_table, file="2015-10-27-MF-GONAD-detags-edgeR.csv",quote=F, sep=",")

# Export the cpm values for all genes
# prior.count=0.125 is the default for exactTest() and is needed here to make values comparable.
gonad_expr_norm <- cpm(gonad_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125)
write.table(gonad_expr_norm, file="2015-10-27-MF-GONAD-CPM-NORM-log2-edgeR.csv",quote=F, sep=",")

# RPKM calculation
# Gene length based on combined exonic length
#IMPORTANT: NEED TO REMOVE HEADER V1!
gene_length <- read.table("2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length")
expressed_genes <- rownames(gonad_expr)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
dim(gene_length)
# Reorder gene length vector to match order in expr object
gene_length_ordered <- gene_length[match(rownames(gonad_expr),gene_length$V1),]
# Check if reorder has worked ( should return TRUE)
all(gene_length_ordered$V1 == rownames(gonad_expr))
gonad_expr_norm_rpkm <- rpkm(gonad_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125, gene.length=gene_length_ordered$V2)
write.table(gonad_expr_norm_rpkm, file="2015-10-27-MF-GONAD-RPKM-NORM-log2-edgeR.csv",quote=F, sep=",")

# =============================================
# Normalize spleen tissue
# =============================================
spleen_data <- read.csv("../count-extraction/2015-10-27-MF-SPLEEN-COUNTS-COMBINED.expr", header=T, row.names=1, stringsAsFactors=FALSE)
dim(spleen_data)
names(spleen_data)
# cor(data)
class(spleen_data[,1])
class(spleen_data[,2])
class(spleen_data[,3])
class(spleen_data[,4])
class(spleen_data[,5])
class(spleen_data[,6])
class(spleen_data[,7])
class(spleen_data[,8])
spleen_conditions <- factor(c("M","M","M","M","F","F","F","F"))
spleen_expr <- DGEList(counts=spleen_data, group=spleen_conditions)
#Default is TMM
spleen_expr <- calcNormFactors(spleen_expr)
# Calculate common and tagwise dispersion at the same time (recommended)
spleen_expr <- estimateDisp(spleen_expr)
# Compare Male - Female expression
# Print expression samples
spleen_expr$samples
spleen_de <- exactTest(spleen_expr, pair=c("M","F"))
#plotMDS(expr)
# Save all tags in this table. nrow(expr) makes sure that all genes are
# included. We need to select the $table, because topTags also contains other
# elements like comparison, test
spleen_detag_table <- topTags(spleen_de, n=nrow(spleen_expr))$table
write.table(spleen_detag_table, file="2015-10-27-MF-SPLEEN-detags-edgeR.csv",quote=F, sep=",")

# Export the cpm values for all genes
# prior.count=0.125 is the default for exactTest() and is needed here to make values comparable.
spleen_expr_norm <- cpm(spleen_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125)
write.table(spleen_expr_norm, file="2015-10-27-MF-SPLEEN-CPM-NORM-log2-edgeR.csv",quote=F, sep=",")

# RPKM calculation
# Gene length based on combined exonic length
#IMPORTANT: NEED TO REMOVE HEADER V1!
gene_length <- read.table("2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length")
expressed_genes <- rownames(spleen_expr)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
dim(gene_length)
# Reorder gene length vector to match order in expr object
gene_length_ordered <- gene_length[match(rownames(spleen_expr),gene_length$V1),]
# Check if reorder has worked ( should return TRUE)
all(gene_length_ordered$V1 == rownames(spleen_expr))
spleen_expr_norm_rpkm <- rpkm(spleen_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125, gene.length=gene_length_ordered$V2)
write.table(spleen_expr_norm_rpkm, file="2015-10-27-MF-SPLEEN-RPKM-NORM-log2-edgeR.csv",quote=F, sep=",")

# =============================================
# Normalize heart tissue
# =============================================
heart_data <- read.csv("../count-extraction/2015-10-27-MF-HEART-COUNTS-COMBINED.expr", header=T, row.names=1, stringsAsFactors=FALSE)
dim(heart_data)
names(heart_data)
# cor(data)
class(heart_data[,1])
class(heart_data[,2])
class(heart_data[,3])
class(heart_data[,4])
class(heart_data[,5])
class(heart_data[,6])
class(heart_data[,7])
class(heart_data[,8])
heart_conditions <- factor(c("M","M","M","M","F","F","F","F"))
heart_expr <- DGEList(counts=heart_data, group=heart_conditions)
#Default is TMM
heart_expr <- calcNormFactors(heart_expr)
# Calculate common and tagwise dispersion at the same time (recommended)
heart_expr <- estimateDisp(heart_expr)
# Print expression samples
heart_expr$samples
# Compare Male - Female expression
# . Note
# that the first group listed in the pair is the baseline for the comparison—so if58 exactTest
# the pair is c("A","B") then the comparison is B - A, so genes with positive
# log-fold change are up-regulated in group B compared with group A (and vice
# versa for genes with negative log-fold change).
heart_de <- exactTest(heart_expr, pair=c("M","F"))
#plotMDS(expr)
# Save all tags in this table. nrow(expr) makes sure that all genes are
# included. We need to select the $table, because topTags also contains other
# elements like comparison, test
heart_detag_table <- topTags(heart_de, n=nrow(heart_expr))$table
write.table(heart_detag_table, file="2015-10-27-MF-HEART-detags-edgeR.csv",quote=F, sep=",")

# Export the cpm values for all genes
# prior.count=0.125 is the default for exactTest() and is needed here to make values comparable.
heart_expr_norm <- cpm(heart_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125)
write.table(heart_expr_norm, file="2015-10-27-MF-HEART-CPM-NORM-log2-edgeR.csv",quote=F, sep=",")

# RPKM calculation
# Gene length based on combined exonic length
#IMPORTANT: NEED TO REMOVE HEADER V1!
gene_length <- read.table("2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length")
expressed_genes <- rownames(heart_expr)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
dim(gene_length)
# Reorder gene length vector to match order in expr object
gene_length_ordered <- gene_length[match(rownames(heart_expr),gene_length$V1),]
# Check if reorder has worked ( should return TRUE)
all(gene_length_ordered$V1 == rownames(heart_expr))
heart_expr_norm_rpkm <- rpkm(heart_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125, gene.length=gene_length_ordered$V2)
write.table(heart_expr_norm_rpkm, file="2015-10-27-MF-HEART-RPKM-NORM-log2-edgeR.csv",quote=F, sep=",")
# =============================================
# Normalize liver tissue, excluded wonky female liver sample
# =============================================
liver_data <- read.csv("../count-extraction/2015-10-27-MF-LIVER-COUNTS-COMBINED.expr", header=T, row.names=1, stringsAsFactors=FALSE)
dim(liver_data)
names(liver_data)
# cor(data)
class(liver_data[,1])
class(liver_data[,2])
class(liver_data[,3])
class(liver_data[,4])
class(liver_data[,5])
class(liver_data[,6])
class(liver_data[,7])
liver_conditions <- factor(c("M","M","M","M","F","F","F"))
liver_expr <- DGEList(counts=liver_data, group=liver_conditions)
#Default is TMM
liver_expr <- calcNormFactors(liver_expr)
# Calculate common and tagwise dispersion at the same time (recommended)
liver_expr <- estimateDisp(liver_expr)
# Print expression samples
liver_expr$samples
# Compare Male - Female expression
liver_de <- exactTest(liver_expr, pair=c("M","F"))
#plotMDS(expr)
# Save all tags in this table. nrow(expr) makes sure that all genes are
# included. We need to select the $table, because topTags also contains other
# elements like comparison, test
liver_detag_table <- topTags(liver_de, n=nrow(liver_expr))$table
write.table(liver_detag_table, file="2015-10-27-MF-LIVER-detags-edgeR.csv",quote=F, sep=",")

# Export the cpm values for all genes
# prior.count=0.125 is the default for exactTest() and is needed here to make values comparable.
liver_expr_norm <- cpm(liver_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125)
write.table(liver_expr_norm, file="2015-10-27-MF-LIVER-CPM-NORM-log2-edgeR.csv",quote=F, sep=",")

# Gene length based on combined exonic length
#IMPORTANT: NEED TO REMOVE HEADER V1!
gene_length <- read.table("2015-10-29-Gallus_gallus.Galgal4.75.exonic.gene.length")
expressed_genes <- rownames(liver_expr)
gene_length <- subset(gene_length, V1 %in% expressed_genes)
dim(gene_length)
# Reorder gene length vector to match order in expr object
gene_length_ordered <- gene_length[match(rownames(liver_expr),gene_length$V1),]
class(gene_length_ordered$V2)
# Check if reorder has worked ( should return TRUE)
all(gene_length_ordered$V1 == rownames(liver_expr))
liver_expr_norm_rpkm <- rpkm(liver_expr, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.125, gene.length=gene_length_ordered$V2)
write.table(liver_expr_norm_rpkm, file="2015-10-27-MF-LIVER-RPKM-NORM-log2-edgeR.csv",quote=F, sep=",")
