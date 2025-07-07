#Library package
library(Seurat)
library(patchwork)
library(dplyr)
library(ggsci)
library(cowplot)
library(ggplot2)
library(DoubletFinder)
#input 
pbmc <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200)
pbmc
pbmc@meta.data
pbmc@assays$RNA
table(pbmc.combined$orig.ident)
all.genes <- rownames(pbmc.combined)
pbmc.combined <- ScaleData(pbmc.combined, features = all.genes)
FindVariableFeatures(pbmc.combined)
#filtration
pbmc <- PercentageFeatureSet(pbmc.combined, pattern = "^MT-", col.name = "percent.mt")
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)
pbmc <- RunPCA(pbmc, verbose = FALSE)
DimPlot(pbmc, label = TRUE)
DimPlot(pbmcplot, label = TRUE, cols = custom_colors)
#SubSertoli cells
subSertoli<-subset(pbmcplot,idents=c("Sertoli cells"))
subSertoli <- ScaleData(subSertoli, verbose = FALSE)
subSertoli <- FindNeighbors(subSertoli, reduction = "pca", dims = 1:15)
subSertoli <- FindClusters(subSertoli, resolution = 0.5)
subSertoli <- RunUMAP(subSertoli, reduction = "pca", dims = 1:15)
DimPlot(subSertoli, label = TRUE)
pd <- new('AnnotatedDataFrame', data = pbmc3@meta.data)
fData <- data.frame(gene_short_name = row.names(subSertoli), row.names = row.names(subSertoli))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle <- newCellDataSet(data,
                          phenoData = pd,
                          featureData = fd,
                          lowerDetectionLimit = 0.5,
                          expressionFamily = negbinomial.size())

monocle <- estimateSizeFactors(monocle)
monocle <- estimateDispersions(monocle)
monocle <- detectGenes(monocle,min_expr = 0.1)
print(head(fData(monocle)))
#Functional Enrichment Analysis
kegg <- data.frame(data$Description,data$count ,data$pval)
colnames(kegg) <- c("Description","count","pval")
p <- ggplot(data=kegg,aes(x=Description,y=count,fill=pval))
p1 <- p + geom_bar(stat="identity") + coord_flip()
p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))
p3 <- p2 + ylim(0,60) +xlim(0,60)+ scale_fill_gradient(low="#FDA7DF",high="#d1d8e0")
p4 <- p3 + scale_x_discrete(limits=rev(kegg[,1])) +labs(x="",y="",title="KEGG")
print(p4)
dev.off()



