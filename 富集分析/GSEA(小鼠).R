rm(list = ls())
library(tidyverse)
library(enrichplot)
library(clusterProfiler)
library(ggnewscale)
library(org.Mm.eg.db)

different_geneset <- read.csv( "file.csv") 
gene <- different_geneset[,c(1,2)]
colnames(gene)[1] <-"SYMBOL"
GO_database <- 'org.Mm.eg.db'
KEGG_database <- 'mmu' #KEGG分析指定物种
ENTREZID_id <- bitr(gene$SYMBOL,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)#gene ID转换

info_merge <- merge(ENTREZID_id,gene,by='SYMBOL')#合并转换后的基因ID和Log2FoldChange
GSEA_input <- info_merge$log2FC
names(GSEA_input) = info_merge$ENTREZID
GSEA_input = sort(GSEA_input, decreasing = TRUE) #降序排列
GSEA_KEGG <- gseKEGG(GSEA_input, 
                     organism = KEGG_database, 
                     pvalueCutoff = 0.05)#GSEA富集分析
write.csv(GSEA_KEGG,file = "result_GSEA_KEGG.csv")

GSEA_GO <- gseGO(GSEA_input, 
                     OrgDb = 'org.Mm.eg.db',
                     pvalueCutoff = 0.05)#GSEA富集分析
write.csv(GSEA_GO,file = "result_GSEA_GO.csv")