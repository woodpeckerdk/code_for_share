#加载包
rm(list = ls())
library(tidyverse)
library(Seurat)
library(patchwork)

# 加载单细胞数据集,创建Seurat对象
sc.data <- Read10X(data.dir = "/home/woodpecker/R/DATA/bladder_cancer/ZFS")
sc_obj <- CreateSeuratObject(counts = sc.data, project = "ZFS", min.cells = 3, min.features = 200)
sc_obj
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
# 将QC指标可视化为小提琴图
VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#FeatureScatter通常用于可视化特征-特征关系，
#但可用于对象计算的任何内容，即对象元数据中的列、PC分数等。
plot1 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("QC1.pdf")
#根据以上3个QC指标筛选合适的细胞和基因，LogNormalizeData规范化为data,ScaleData归一化为scale存储在[[“RNA”]]中，#高度可变特征的识别（特征选择），并将可变基因可视化。
sc_obj <- subset(sc_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 25)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sc_obj), 10)  # 识别前10个高可变基因
plot1 <- VariableFeaturePlot(sc_obj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave(plot = plot1 +plot2,"高可变基因.pdf")
all.genes <- rownames(sc_obj)
sc_obj <- ScaleData(sc_obj, features = all.genes)


#根据高可变基因特征进行线性降维,
sc_obj <- RunPCA(sc_obj, features = VariableFeatures(object = sc_obj))
# 通过几种不同的方式检查和可视化PCA结果
print(sc_obj[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(sc_obj, dims = 1:2, reduction = "pca")  #点图可视化PCA结果
DimPlot(sc_obj, reduction = "pca") + NoLegend()       #DimPlot可视化
DimHeatmap(sc_obj, dims = 1:9, cells = 500, balanced = TRUE)#热图可视化PCA结果（PC1-9）
ElbowPlot(sc_obj)  #确定数据集的维度

#降维聚类，
sc_obj <- FindNeighbors(sc_obj, dims = 1:10)  #细化边缘权重
sc_obj <- FindClusters(sc_obj, resolution = 0.5)  #迭代地将细胞分组在一起，目的是优化标准模块化函数。
sc_obj <- RunUMAP(sc_obj, dims = 1:10)
DimPlot(sc_obj, reduction = "umap",label = T)
 
#找到单个cluster 2的marker、找两个cluster之间的marker，找所有cluster之间的marker
cluster2.markers <- FindMarkers(sc_obj, ident.1 = 2) 
cluster5.markers <- FindMarkers(sc_obj, ident.1 = 5, ident.2 = c(0, 3))
sc_obj_markers <- FindAllMarkers(sc_obj, only.pos = TRUE)
#将找到的marker进行筛选过滤
sc_obj_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)
cluster0_markers <- FindMarkers(sc_obj, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
#小提琴图、featureplot可视化特征表达
VlnPlot(sc_obj, features = c("MS4A1", "CD79A"))
VlnPlot(sc_obj, features = c("MS4A1","CD79A"), layer = "counts", log = TRUE)#指定layer中的counts/data/scale
FeaturePlot(sc_obj, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A"))
#top10基因热图可视化
sc_obj_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(sc_obj, features = top10$gene) + NoLegend()
#定义feature,气泡图可视化feature的表达情况
genes_to_check <- c('ENPP3','FCER1A','FCGR1A','FCGR2A','FCGR2B')
DotPlot(sc_obj,features = genes_to_check,
             assay = "RNA",group.by = 'seurat_clusters') + coord_flip() +ggtitle(label = T)
#参考：
genes_to_check <- c('PTPRC','CD3D','CD3E','CD4','CD8A', ##T cells
                    'ENPP3','FCER1A','FCGR1A','FCGR2A','FCGR2B',  ##Mast cells
                    'ACTA2','PDGFRA','PDGFRB','THY1', #Fibroblast
                    'CD19','CD79A','MS4A1', ## B cells
                    'TAGLN2','CD5', ## Actived B cells
                    'CD27','CD38','LY9','LATR1','ICAM1','KIT',  ## plasma
                    'NKG7','GNLY', #NK cells
                    'CXCR5','CCR7',  ## Memory T cell
                    'CD6','IL7R','IL2RA','IKZF2', ##Treg
                    'CD33','ENTPD1', ##MDSC
                    'S100A8','S100A9','S100A12', ## Neutrophil
                    'CD68','CD163','MRC1','MSR1','CXCL10','CCL18', ## Macrophage
                    'PECAM1','VWF','MCAM','CD34','ESM1', ## Endothelial
                    'ALDH1A1','KRT18','PROM1', )## Cancer stem cell
########  第三部分 细胞注释  ########
celltype <- data.frame(clusterID = 0:14,
                      celltype = "unknown")
celltype[celltype$clusterID %in% c(0,1,2,3),2] = "T cells"
celltype[celltype$clusterID %in% c(4),2] = "NK"
celltype[celltype$clusterID %in% c(5,6,7),2] ="M1"
celltype[celltype$clusterID %in% c(8,9),2] = "Plasma"
celltype[celltype$clusterID %in% c(10),2] = "actived B cell"
celltype[celltype$clusterID %in% c(11,12,13,14),2] = "Monocyte"
celltype
table(celltype$celltype)
sce_in <- sc_obj
sce_in@meta.data$celltype <-  "NA"
##注释
for(i in 1:nrow(celltype)){
  sce_in@meta.data[which(sce_in@meta.data$seurat_clusters == celltype$clusterID[i]),"celltype"] <- celltype$celltype[i]}
table(sce_in@meta.data$celltype)
sce <- sce_in
p.umap.cell <- DimPlot(sce,reduction = "umap",group.by = "celltype",label = T)
p.umap.cell
ggsave(plot = p.umap.cell,filename = "umap_celltype.pdf",width = 10,height = 10)
DefaultAssay <- "RNA"  #DefaultAssay设置为RNA意味着接下来的分析将基于原始值
saveRDS(sce,file = "after_annotation.rds")
