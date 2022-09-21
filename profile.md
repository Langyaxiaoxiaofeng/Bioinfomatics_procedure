
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

##