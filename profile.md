## 读取数据
### h5文件
```R
data_sample <- Read10X_h5("sample.h5) #读取数据
seurat_sample <- CreateSeuratObject(counts = data_sample, project = "sample", min.cells = 3, min.feature = 200) #创建seurat对象
```
### 10X数据
##### 10X数据解压后有barcodes.tsv、genes.tsv、matrix.mtx三个文件
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
data <- NormalizeData(object = seurat_filter, normalization.method = "LogNormalize", scale.factor = 10000) #数据标准化
data <- FindVariableFeatures(object = data, selection.method = "vst", nfeatures = 1500) #找高变基因
```




