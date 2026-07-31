#Quality control of Danielli intergrated dataset split by study.
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(scCustomize)

RMS <- readRDS("RMS_atlas_final_20240130.rds")

#gives you number of genes and cells for the whole integrated dataset.
DefaultAssay(RMS) <- "RNA"
dim(RMS) #63174 genes and 107523 cells overall. 

#See metadata and already has n_count, n_feature and percent.mt.
head(RMS@meta.data)

#Shows the different labels in the origin column which gives names of studies. 
unique(RMS$origin)

#Violin plot number of features
VlnPlot(RMS, features = "nFeature_RNA", group.by = "origin", pt.size = 0, combine = TRUE) + NoLegend() +
  ggplot2::labs(x= "Dataset", y= "Number of Features", title = NULL)

#violin plot of percentage of mitochondrial genes. 
VlnPlot(RMS, features = "percent.mt", group.by = "origin", pt.size = 0, combine = TRUE) + NoLegend() +
  ggplot2::labs(x= "Dataset", y= "Percent mitochondrial genes (%)", title = NULL)

unique(RMS@assays)

DefaultAssay(RMS) <- "integrated"

#Due to the violin plot of percent mito having some cells reaching over the 20% filter,
#I wanted to check how many of these cells there were and whether I would need to reintegrate the data.
#as was only a very small percentage it is just noise and nothing to worry about. 
table(RMS$percent.mt > 20)
mean(RMS$percent.mt > 20) *100 #0.0177% 

RMS.fil <- subset(RMS, subset = (origin == "Patel et al." & percent.mt <= 15) | (origin != "Patel et al." & percent.mt <= 20))

# Check cluster proportions before vs after
# Before
p1 <- DimPlot(RMS, label = TRUE) + ggtitle("Before filtering")

# After
p2 <- DimPlot(RMS.fil, label = TRUE) + ggtitle("After filtering (>20% MT removed)")

p1 + p2

# Did any clusters disappear or shrink dramatically?
table(Idents(RMS)) - table(Idents(RMS.fil))

#Violin plot of mito genes after filtering.
VlnPlot(RMS.fil, features = "percent.mt", group.by = "origin", pt.size = 0, combine = TRUE) + NoLegend() +
  ggplot2::labs(x= "Dataset", y= "Percent mitochondrial genes (%)", title = NULL)
