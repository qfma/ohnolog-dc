# =============================================
# Hiearchical clusters of spleen samples
# =============================================
data <- read.csv("2015-10-27-MF-SPLEEN-CPM-NORM-log2-edgeR.csv", header=T, row.names=1)
# data <- read.csv("2015-10-27-MF-SPLEEN-RPKM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
cor(data)
plot(hclust(dist(t(data))))
png(file="2015-10-27-MF-SPLEEN-CPM-NORM-log2-edgeR-hclust.png")
# png(file="2015-10-27-MF-SPLEEN-RPKM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()

# =============================================
# Hiearchical clusters of gonad samples
# =============================================
data <- read.csv("2015-10-27-MF-GONAD-CPM-NORM-log2-edgeR.csv", header=T, row.names=1)
# data <- read.csv("2015-10-27-MF-GONAD-RPKM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
head(data)
cor(data)
plot(hclust(dist(t(data))))
png(file="2015-10-27-MF-GONAD-CPM-NORM-log2-edgeR-hclust.png")
# png(file="2015-10-27-MF-GONAD-RPKM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()

# =============================================
# Hiearchical clusters of heart samples
# =============================================
data <- read.csv("2015-10-27-MF-HEART-CPM-NORM-log2-edgeR.csv", header=T, row.names=1)
# data <- read.csv("2015-10-27-MF-HEART-RPKM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
cor(data)
plot(hclust(dist(t(data))))
png(file="2015-10-27-MF-HEART-CPM-NORM-log2-edgeR-hclust.png")
# png(file="2015-10-27-MF-HEART-RPKM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()
# =============================================
# Hiearchical clusters of liver samples
# =============================================
data <- read.csv("2015-10-27-MF-LIVER-CPM-NORM-log2-edgeR.csv", header=T, row.names=1)
# data <- read.csv("2015-10-27-MF-LIVER-RPKM-NORM-log2-edgeR.csv", header=T, row.names=1)
dim(data)
names(data)
cor(data)
plot(hclust(dist(t(data))))
png(file="2015-10-27-MF-LIVER-CPM-NORM-log2-edgeR-hclust.png")
# png(file="2015-10-27-MF-LIVER-RPKM-NORM-log2-edgeR-hclust.png")
plot(hclust(dist(t(data))))
dev.off()