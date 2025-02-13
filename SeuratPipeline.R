library(Seurat)

nodData= readRDS("./data/cleanSeurat.RDS")

rawCounts<- nodData@assays$RNA@counts
cleanMetaData= nodData@meta.data
cleanMetaData= cleanMetaData[, c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "sample type", "attackedornot", "Cell type")]

nodData= CreateSeuratObject(counts = rawCounts, meta.data = cleanMetaData)

saveRDS(nodData,"./data/rawSeurat.RDS")


nodData=SCTransform(nodData)
nodData <- RunPCA(nodData, verbose = FALSE)
nodData <- RunUMAP(nodData, dims = 1:30, verbose = FALSE)
saveRDS(nodData,"./data/UMAP.RDS")

nodData= readRDS("./data/UMAP.RDS")

metadata= nodData@meta.data

DimPlot(nodData)