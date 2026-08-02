library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(scCustomize)
library(plyr)
library(viridis)

RMS <- readRDS("RMS_atlas_final_20240130.rds")
rm(RMS)

RMS.fil <- subset(RMS, subset = (origin == "Patel et al." & percent.mt <= 15) | (origin != "Patel et al." & percent.mt <= 20))

RMS.wei.pat.fil <- subset(RMS.fil, subset = (origin== "Wei et al." & nFeature_RNA >=1000) |
                            (origin == "Patel et al." & nFeature_RNA >= 400) |
                            (origin == "Weng et al.") |
                            (origin == "Danielli et al."))

dim(RMS.wei.pat.fil) #106842 cells now after filtering. 

RMS.wei.pat.fil$newfusion <- mapvalues(x= RMS$fusion, from= c("PAX3::FOXO1", "PAX7::FOXO1"),
                                       to= c("FP-RMS", "FP-RMS"))

RMS.wei.pat.fil.noMYOD1 <- subset(RMS.wei.pat.fil, subset = newfusion != "MYOD1")

unique(RMS.wei.pat.fil.noMYOD1$newfusion)

DefaultAssay(RMS.wei.pat.fil.noMYOD1) <- "RNA"

NKligs.inhib <- c("CD274","PDCD1LG2","CLEC2D","LGALS9","ENTPD1", "NT5E","PVR", "PVRL2")

plots.inhib <- lapply(NKligs.inhib, function(g) {FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = g,
                                                             min.cutoff = 0, max.cutoff = "q95") })
wrap_plots(plots.inhib) + plot_layout(guides = "collect")

rm(plots.inhib)

RMS.wei.pat.fil.noMYOD1@active.assay

FeaturePlot(RMS.wei.pat.active.assayFeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "PVRL2", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "LGALS9", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "CD274", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "PDCD1LG2", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "CLEC2D", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "ENTPD1", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "NT5E", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "PVR", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "PVRL2", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "CD48", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "CD58", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "NCR3LG1", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "ICAM1", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "MICA", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "MICB", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "ULBP1", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "ULBP2", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "ULBP3", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "RAET1G", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "RAET1E", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "RAET1L", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

#Why is it that on the RMS object I can't find HLA-E or A/B/C but on this filtered version I can? Not even sure if I should have filtered it?
FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "HLA-E", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "HHLA-2", cols = c("white", "navy"), min.cutoff = 0, max.cutoff =  "q98", pt.size = 1)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = NKligs.inhib)

table(RMS.wei.pat.fil.noMYOD1$origin)

#FNFP UMAP
Idents(RMS.wei.pat.fil.noMYOD1) <- RMS.wei.pat.fil.noMYOD1$newfusion
DimPlot(RMS.wei.pat.fil.noMYOD1, reduction = "umap_rpca", pt.size = 1, cols = c("#FF8247", 
                                                                                
                                                                                "#7A67EE"))

#FN PAX3 and 7 UMAP
Idents(RMS.wei.pat.fil.noMYOD1) <- RMS.wei.pat.fil.noMYOD1$fusion
DimPlot(RMS.wei.pat.fil.noMYOD1, reduction = "umap_rpca", pt.size = 1, cols = c("#FF8247",
                                                                        
                                                                                "#7A67EE",
                                                                                "#00688B"))


RMS.wei.pat.fil.noMYOD1$fusion <- factor(RMS.wei.pat.fil.noMYOD1$fusion,
                                         levels = c("FN-RMS", "PAX3::FOXO1", "PAX7::FOXO1"))

table(RMS.wei.pat.fil.noMYOD1$fusion)
class(RMS.wei.pat.fil.noMYOD1$fusion)
levels(RMS.wei.pat.fil.noMYOD1$fusion)
#Module scoring FN and FP onto UMAP
RMS.wei.pat.fil.noMYOD1 <- AddModuleScore(RMS.wei.pat.fil.noMYOD1, features = NK.inhib.ligs,
                                          name = "Inhibitory NK score")

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "Inhibitory NK score1", split.by = "newfusion",
            pt.size = 0.5) + theme(legend.position = "right")

NK.inhib.ligs <- list(c("CD274","PDCD1LG2","CLEC2D","LGALS9","ENTPD1", "NT5E","PVR", "PVRL2"))

head(RMS.wei.pat.fil.noMYOD1$`Inhibitory NK score1`)

FeaturePlot(RMS.wei.pat.fil.noMYOD1, features = "Inhibitory NK score1")

RMS.wei.pat.fil.noMYOD1$newfusion <- as.factor(RMS.wei.pat.fil.noMYOD1$newfusion)
