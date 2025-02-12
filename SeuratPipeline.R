library(Seurat)

nodData= readRDS("./data/cleanSeurat.RDS")

rawCounts<- nodData@assays$RNA@counts
cleanMetaData= nodData@meta.data
cleanMetaData= cleanMetaData[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "sample type", "attackedornot")]

nodData= CreateSeuratObject(counts = rawCounts, meta.data = cleanMetaData)

saveRDS(nodData,"./data/rawSeurat.RDS")



