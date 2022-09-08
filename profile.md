## 读取数据
#### h5文件
```R
data_sample <- Read10X_h5("sample.h5) #读取数据
seurat_sample <- CreateSeuratObject(counts = data_sample, project = "sample", min.cells = 3, min.feature = 200) #创建seurat对象
```
#### 10X数据
###### 10X数据解压后有barcodes.tsv、genes.tsv、matrix.mtx三个文件
```R
pbmc.data <- Read10X(data.dir = "dir_path") #读取数据
seurat_sample <- CreateSeuratObject(counts = pbmc.data, project = "sample", min.cells = 3, min.features = 200) #创建seurat对象
```

## 合并seurat对象
```R
seurat_merge <- merge(seurat_sample1, seurat_sample2, add.cell.ids = c("sample1", "sample2"))
```

## 数据过滤
```R
seurat_sample[["percent.mt"]] <- PercentageFeatureSet(object = seurat_sample, pattern = "^MT-") #计算线粒体基因含量
seurat_filter=subset(x = seurat_sample, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10) #设置过滤条件对数据进行过滤
```

## 画图1
```R
VlnPlot(object = seurat_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) #小提琴图
plot1 <- FeatureScatter(object = seurat_filter, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=1.5)
plot2 <- FeatureScatter(object = seurat_filter, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=1.5)
CombinePlots(plots = list(plot1, plot2)) #画散点图
```

## 数据进一步处理
```R
# NormalizeData可消除不同测序深度的影响(应该是分位数标准化)
suerat_data <- NormalizeData(object = seurat_filter, normalization.method = "LogNormalize", scale.factor = 10000) #数据标准化
seurat_data <- FindVariableFeatures(object = seurat_data, selection.method = "vst", nfeatures = 1500) #找高变基因
# 查看前10的高变基因
top10 <- head(x = VariableFeatures(object = seurat_data), 10)
plot1 <- VariableFeaturePlot(object = seurat_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
```

## 数据分析
#### PCA
```R
# ScaleData对基因表达量的数值进行z-score转换，使其符合norm(mean = 0, sd = 1)
seurat_data = ScaleData(seurat_data)          
seurat_data = RunPCA(object = seurat_data, npcs = 20, pc.genes=VariableFeatures(object = seurat_data)) #PCA
VizDimLoadings(object = seurat_data, dims = 1:4, reduction = "pca",nfeatures = 20) #可视化对每个主成分影响比较大的基因集
DimPlot(object = seurat_data, reduction = "pca") #查看主成分坐标下的散点图
DimHeatmap(object = seurat_data, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2) #查看主成分坐标下的热图
# 筛选信息量大的主成分
# 随机置换数据的一部分子集（默认1%）再运行PCA，构建了一个’null distribution’的特征分数，重复这一步。最终会识别出低p-value特征的显著PCs
seurat_data <- JackStraw(object = seurat_data, num.replicate = 100) 
seurat_data <- ScoreJackStraw(object = seurat_data, dims = 1:15)
JackStrawPlot(object = seurat_data, dims = 1:15)
ElbowPlot(seurat_data) #（肘方法）展示PCs的信息量
pcSelect = 16 #选择合适数量的主成分
```

#### TSNE
```R
seurat_data <- FindNeighbors(object = seurat_data, reduction = "pca", dims = 1:pcSelect)       
seurat_data <- FindClusters(object = seurat_data, resolution = 0.5)
write.table(seurat_data$seurat_clusters,file="Cluster.txt",quote=F,sep="\t",col.names=F)        
seurat_data <- RunTSNE(object = seurat_data, dims = 1:pcSelect, check_duplicates = FALSE) #大的数据集可能发生自然重复，所以check_duplicates = FALSE
write.table(Embeddings(object = seurat_data[["tsne"]]),file="tsneaxis.txt",quote=F,sep="\t",col.names=T)
TSNEPlot(object = seurat_data, pt.size = 0.5, label = TRUE) 
```

#### UMAP
```R
seurat_data <- RunUMAP(object = seurat_data, dims = 1:pcSelect)
UMAPPlot(object = seurat_data, pt.size = 0.5, label = TRUE)
write.table(Embeddings(object = seurat_data[["umap"]]),file="umapaxis.txt",quote=F,sep="\t",col.names=T)
```

## 细胞注释
#### 找marker基因
```R
# 找每个cluster的marker基因
seurat.markers <- FindAllMarkers(object = seurat_data,
                               only.pos = FALSE, #Only return positive markers
                               min.pct = 0.25, #只测试在两个cluster中检测到至少在min.pct%中表达的基因。通过不测试那些很少表达的基因来加快速度
                               logfc.threshold = logFCfilter)
sig.markers <- seurat.markers[(abs(as.numeric(as.vector(seurat.markers$avg_log2FC))) > logFCfilter &
							   as.numeric(as.vector(seurat.markers$p_val_adj)) < adjPvalFilter),]
write.table(sig.markers, file="clusterMarkers.txt", sep="\t", row.names=F, quote=F)
top10.markers <- seurat.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object = seurat_data, features = top10.markers$gene) + NoLegend()
VlnPlot(object = seurat_data, features = row.names(sig.markers)[1:2])
```

#### 展示感兴趣的基因
```R
for(showGene in showGenes){
  print(showGene)
  p <- FeaturePlot(object = seurat_data, reduction = "tsne", features = showGene,
  				   cols = c("grey", "red"), pt.size = 0.5, order = TRUE)
  print(p)
}

FeaturePlot(object = seurat_data, reduction = "tsne", features = showGenes, 
			cols = c("#000081","#6b2588","#d59c2a","#af1a1a"),
			pt.size = 0.5,
            order = TRUE)

DotPlot(object = seurat_data, features = showGenes)

avg_exp <- AverageExpression(seurat_data, verbose=T) #verbose是日志显示选项
avg_exp <- avg_exp$RNA
avg_exp <- t(avg_exp[showGenes,])
pheatmap(avg_exp, cluster_rows=F, cluster_cols=F, scale="column", fontsize=40)
```
#### 根据之前的展示对细胞进行手动注释
###### 1、直接对seurat_data注释
```R
new.cluster.cell <- c("cell.type1", "cell.type2", "cell.type3")
seurat_data@meta.data$celltype <- seurat_data@meta.data$seurat_clusters
levels(seurat_data@meta.data$celltype) <- new.cluster.cell
```
###### 2、用RenameIdents注释
```R
new.cluster.cell <- c("cell.type1", "cell.type2", "cell.type3")
names(new.cluster.cell) <- levels(seurat_data)
seurat_data <- RenameIdents(seurat_data, new.cluster.cell)
```

## 拟时分析
```R

```
