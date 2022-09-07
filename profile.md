## 读取数据
### h5文件
```R
data_sample1 <- Read10X_h5("sample1.h5) #读取数据
seurat_sample1 <- CreateSeuratObject(data_sample1, project = "data_sample") #创建seurat对象
```

