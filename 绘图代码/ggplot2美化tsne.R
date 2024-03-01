#使用conda env: r4.3.0，r4.3.2
rm(list = ls())
library(SingleR)
library(Seurat)
library(tidyverse)
library(ggrepel)
#加载seurat对象
data <- readRDS("RData_and_rds/after_annotation.rds")
data
# seurat 包中`DimPlot`函数一行代码绘制umap图
DimPlot(data, group.by = c("celltype"),reduction = "tsne")
ggsave("orign_tsne.pdf")
# 二 ggplot2绘制umap图
# umap图所需的数据就是每个cell的坐标以及cluster或者celltype信息，然后绘制点图
tsne = data@reductions$tsne@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = data@meta.data$celltype) # 注释后的label信息 ，改为cell_type
head(tsne)
# 2.2 ggplot2 绘制umap图
# 调整color，颜色列表来自于https://www.jianshu.com/p/67d2decf5517
allcolour=c("#267365","#BC28C7","#F29F05","#F28705","#F23030","#F266AB")
p <- ggplot(tsne,aes(x= tSNE_1 , y = tSNE_2 ,color = cell_type)) + 
geom_point(size = 0.8 , alpha =1 )  +  scale_color_manual(values = allcolour)
ggsave("ggplot_tsne.pdf",width = 10,height = 7)
# 3.1 调整umap图 - theme
# 主题的调整比较简单，去掉网格线，坐标轴和背景色即可
p2 <- p  +
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))
ggsave("ggplot_tsne_p2.pdf",width = 10,height = 7)
# 3.2 调整umap图 - legend
# legeng部分去掉legend.title后，调整标签大小，标签点的大小以及 标签之间的距离
p3 <- p2 +         
        theme(
          legend.title = element_blank(), #去掉legend.title 
          legend.key=element_rect(fill='white'), #
        legend.text = element_text(size=20), #设置legend标签的大小
        legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 
ggsave("ggplot_tsne_p3.pdf",width = 10,height = 7)


