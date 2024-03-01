#运行的conda env: r4.3.2
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(ggplot2)#用于ggsave保存pdf
library(ggnewscale)
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
test_markers <- read.csv( "/home/woodpecker/R/R教程/GO_KEGG_GSEA分析/data_file/test_markers.csv") 
info <- test_markers[,c(1,3)]
colnames(info) <- c("SYMBOL","Log2FoldChange") #统一重命名列
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
gene <- bitr(info$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)#gene ID转换
#2.3 GO分析:
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T)
write.csv(GO,file = "result_GO.csv")
#2.4 KEGG分析:
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
write.csv(KEGG,file = "result_KEGG.csv")
# 2.5 GSEA分析:
# 由于GSEA需要差异倍数的信息即Log2FoldChange，我们先要对gene转换后的ID和读入信息进行合并。
info_merge <- merge(info,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$Log2FoldChange
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE) #降序排列
GSEA_KEGG <- gseKEGG(GSEA_input, 
                     organism = KEGG_database, 
                     pvalueCutoff = 0.05)#GSEA富集分析
write.csv(GSEA_KEGG,file = "result_GSEA_KEGG.csv")
#3. GO_KEGG_GSEA基础可视化:
#3.1 GO/KEGG富集柱状图+点状图:
barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
ggsave("GO_bar.pdf")
barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
ggsave("KEGG_bar.pdf")
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
ggsave("GO_dot.pdf")
dotplot(KEGG)
ggsave("KEGG_dot.pdf")
#GSEA可视化
library(ggridges)
ridgeplot(GSEA_KEGG) 
ggsave("ridgeplot.pdf")
gseaplot2(GSEA_KEGG,1)
ggsave("gseaplot.pdf")
gseaplot2(GSEA_KEGG,1:25)#30是根据ridgeplot中有30个富集通路得到的
ggsave("gseaplot25.pdf")