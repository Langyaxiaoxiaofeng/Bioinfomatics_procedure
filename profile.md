
## R包安装
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("minfi", "ChAMPdata", "Illumina450ProbeVariants.db", "sva",
    "IlluminaHumanMethylation450kmanifest", "limma", "RPMM", "DNAcopy", "preprocessCore",
    "impute", "marray", "wateRmelon", "goseq", "plyr", "GenomicRanges", "RefFreeEWAS",
    "qvalue", "isva", "doParallel", "bumphunter", "quadprog", "shiny", "shinythemes",
    "plotly", "RColorBrewer", "DMRcate", "dendextend", "IlluminaHumanMethylationEPICmanifest",
    "FEM", "matrixStats", "missMethyl", "combinat"))
```

## 读取数据
```R
# 矩阵数据
data <- read.table("../../data/data152519/GSE152519_series_matrix.txt.gz", header = T, comment.char = "!") #comment.char = "!" 是为了跳过"!"开头的行
# 接下来整理行名和列名

# 使用getGEO读取或下载数据
gset <- getGEO(GEO = "GSE152519", filename = "destdir/data.txt.gz", destdir = "destdir", GSEMatrix =TRUE, getGPL=FALSE)
# 有时候不能加下面这两行
# if (length(gset) > 1) idx <- grep("GPL21145", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 提取表达矩阵
ex <- exprs(gset)
```
## log2转换（看了半天没看明白）
```R
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { 
  ex[which(ex <= 0)] <- NaN
  ex <- log2(ex) 
}
```

## 看看数据质量
```R
# box-and-whisker plot
pdf(filename = "1.boxplot.pdf", width=3+ncol(gset)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("GSE152519", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

# expression value distribution plot
pdf(file = "1.value_distribution.pdf")
par(mar=c(4,4,2,1))
title <- paste ("GSE152519", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)
dev.off()

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
pdf(file = "1.mean-variance.pdf")
plotSA(lmFit(ex), main="Mean variance trend, GSE152519")
dev.off()

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
library("maptools")  # point labels without overlaps
pdf(file = "1.umap.pdf")
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
dev.off()
```