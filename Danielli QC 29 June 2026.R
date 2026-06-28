#Quality control of Danielli intergrated dataset.
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)

RMS <- readRDS("RMS_atlas_final_20240130.rds")

#gives you number of genes and cells.
dim(RMS)

#shows metadata columns.
head(RMS@meta.data)

RMS$allcells <- "ALL CELLS"

qc <- VlnPlot(RMS, features = "nFeature_RNA", group.by = "allcells", pt.size = 0)
qc
rm(qc)
RMS.QC@assays$
RMS.QC <- RMS
DefaultAssay(RMS.QC) <- "ADT"
RMS.QC <- PercentageFeatureSet(RMS, pattern = "^RP[SL]", col.name = "percent.rb")
head(RMS.QC@meta.data)

RMS.QC <- JoinLayers(RMS.QC)

FeatureScatter(RMS.QC, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "allcells")

FeatureScatter(RMS.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "allcells")

#Violin plots for the four different QC parameters.
VlnPlot(RMS.QC, features = "nFeature_RNA", group.by = "allcells", pt.size = 0)
VlnPlot(RMS.QC, features = "nFeature_RNA", group.by = "allcells")

VlnPlot(RMS.QC, features = "nCount_RNA", group.by = "allcells", pt.size = 0)
VlnPlot(RMS.QC, features = "nCount_RNA", group.by = "allcells")

VlnPlot(RMS.QC, features = "percent.mt", group.by = "allcells", pt.size = 0)
VlnPlot(RMS.QC, features = "percent.mt", group.by = "allcells")

VlnPlot(RMS.QC, features = "percent.rb", group.by = "allcells", pt.size = 0)
VlnPlot(RMS.QC, features = "percent.rb", group.by = "allcells")

#Violin plots by origin.

VlnPlot(RMS.QC, features = "nCount_RNA", group.by = "origin", pt.size = 0)
VlnPlot(RMS.QC, features = "nCount_RNA", group.by = "origin")

VlnPlot(RMS.QC, features = "nFeature_RNA", group.by = "origin", pt.size = 0)
VlnPlot(RMS.QC, features = "nFeature_RNA", group.by = "origin")

VlnPlot(RMS.QC, features = "percent.mt", group.by = "origin", pt.size = 0)
VlnPlot(RMS.QC, features = "percent.mt", group.by = "origin")

VlnPlot(RMS.QC, features = "percent.rb", group.by = "origin", pt.size = 0)
VlnPlot(RMS.QC, features = "percent.rb", group.by = "origin")


doublets <- read.table("data/update/scrublet_calls.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
RMS.QC <- AddMetaData(RMS.QC,doublets)
head(srat[[]])

remotes::install_github("andreaskapou/SeuratPipe@release/v1.0.0")
SeuratPipe::install_scrublet()
library(reticulate)
counts.RMS <- GetAssayData(RMS.QC, layer = "counts")


# Scrublet needs cells as rows and genes as columns
py_run_string("import scrublet as scr")
py_run_string("import scipy.sparse")

# Send the transposed matrix to Python
py$counts_T <- Matrix::t(counts_matrix)

py_run_string("
scrub = scr.Scrublet(counts_T)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
")

#Add percentage ribosomal to RMS.QC.

RMS_RNA <- LayerData(object = RMS, assay = "RNA", layer = "data")
ribo_genes <- grep(pattern = "^RP[LS]", x = rownames(RMS_RNA), value = TRUE)
percent_ribo <- (Matrix::colSums(RMS_RNA[ribo_genes, ]) / Matrix::colSums(RMS_RNA)) * 100
RMS.QC <- AddMetaData(object = RMS_RNA, metadata = percent_ribo, col.name = "percent.rb")
RMS.QC[["percent.rb"]] <- percent_ribo

head(RMS.QC@meta.data)
FeatureScatter(RMS.QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "allcells")
FeatureScatter(RMS.QC, feature1 = "percent.rb", feature2 = "percent.mt", group.by = "allcells")

#doublet scores using RNA layer.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("scran", "scater"))
BiocManager::install("scDblFinder")
library(scran)
library(scDblFinder)

RMS.QC.sce <- as.SingleCellExperiment(RMS.QC)
RMS.QC.sce <- scDblFinder(RMS.QC.sce, samples = "id")

#So RNA data is most likely normalized and not raw so cannot calculate doubelt counts!!!. 
is_integer(RMS_RNA)

RMS_RNA[1:5, 1:5]

summary(as.vector(RMS_RNA[1:100, 1:100]))
