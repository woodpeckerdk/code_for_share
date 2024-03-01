#需要r4.3.0
library(Seurat)
library(tidyverse)
library(scCustomize) # 需要Seurat版本4.3.0
library(viridis)
library(RColorBrewer)
library(gridExtra)
data <- readRDS("RData_and_rds/Myeloid_after_annotation.rds")
data <- data[,data$celltype !="other"]
data$celltype["S100P+ Macro"] <- "MAST4+ Macro"
data$celltype <- str_replace(data$celltype,"S100P\\+ Macro", "MAST4+ Macro")
unique(data$celltype)
saveRDS(data,'Myeloid_after_annotation_modif.rds')
p <- DotPlot(data,features = genes_to_check,
             assay = "RNA",group.by = 'celltype')  +ggtitle(label = T)
ggsave("orign_dotplot.pdf",width = 20,height = 7)
#整理的一些marker可以添加到这里
genes_to_check <- c('IL1B','NLRP3','CCL3',
                    'FOLR2','SLC40A1',
                      'LIPA',
                      'MKI67','TOP2A','UBE2C','TK1',
                      'FCGR2B', 'RGCC', 'IL4I1',"ALDH2","CSF2RA",
                      "SLC14A1","IRF6","MAST4")
# 美化
p1 <- DotPlot(data, features = genes_to_check ,
        assay='RNA' ,group.by = 'celltype') + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ #自定义轴标签
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#E3B527','#E33C22')) #颜色
p1
ggsave("beauty_dotplot.pdf",width = 10,height = 5)
