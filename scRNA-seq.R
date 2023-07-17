library(Seurat)
data_dir <- '/Sshare/home/dengyihan/scRNA-seq/pbmc_1k_v3/outs/filtered_feature_bc_matrix'

# 建立Seurat对象------------------------------------------------
exp <- Read10X(data.dir = data_dir)
exp <- CreateSeuratObject(exp, project = "scRNA", min.cells = 5, min.features = 200)

# 预处理-------------------------------------------------
exp <- PercentageFeatureSet(exp, pattern = "^MT-", col.name = "percent.mt")
exp <- subset(exp, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 20)
exp <- subset(exp, subset = nCount_RNA <= 20000)
exp <- NormalizeData(exp, normalization.method = "LogNormalize", scale.factor = 10000)
exp <- FindVariableFeatures(exp, selection.method = "vst", nfeatures = 2000)
exp <- ScaleData(exp, vars.to.regress = "percent.mt")
exp <- RunPCA(exp, features = VariableFeatures(object = exp))
ElbowPlot(exp)

# 聚类-----------------------------------
exp <- FindNeighbors(exp, dims = 1:20)
exp <- FindClusters(exp, resolution = 0.8)
exp <- RunUMAP(exp, dims = 1:20)
exp <- RunTSNE(exp, dims = 1:20)

# 可视化---------------------------------------------
DimPlot(exp, label = T, label.size = 10)
FeaturePlot(exp, features = c('CD8A', 'CD4'))
VlnPlot(exp, features = c('CD8A', 'CD4'), pt.size = 0)