library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(viridis)

#Loading Danielli data
RMS <- readRDS("RMS_atlas_final_20240130.rds")
unique(RMS$fusion)

#removed MYOD1.
RMS.noMYOD1 <- subset(RMS, subset = fusion != "MYOD1")
unique(RMS.noMYOD1$fusion)

RMS.noMYOD1@active.ident <- factor(RMS.noMYOD1@active.ident,
                           levels = c("PAX3::FOXO1", "PAX7::FOXO1", "FN-RMS"))

levels(RMS.noMYOD1$fusion) <- c("PAX3::FOXO1", "PAX7::FOXO1", "FN-RMS")

#Mapping just inhibitory NK cell ligs onto dotplot by fusion status. 
DefaultAssay(RMS.noMYOD1) <- "RNA"
RMS.noMYOD1@active.assay

NKligs.inhib <- c("CD274","PDCD1LG2","CLEC2D","LGALS9","COL18A1","COL4A1","COL4A2","ENTPD1", "NT5E","PVR", "PVRL2")

NKligs.inhib.nocol <- c("CD274","PDCD1LG2","CLEC2D","LGALS9","ENTPD1", "NT5E","PVR", "PVRL2")

NKligs.inhib.dotplot <- DotPlot(RMS.noMYOD1, features = NKligs.inhib, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.inhib.dotplot
ggsave(NKligs.inhib.dotplot, filename= "C:/Danielli data/Danielli data/Danielli NK inhibitory ligands with collagen NO MYOD1 05.05.26.png")

#Nk inhib ligs without the collagen genes. 
NKligs.inhib.nocol.dotplot <- DotPlot(RMS.noMYOD1, features = NKligs.inhib.nocol, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.inhib.nocol.dotplot

#Mapping just activating NK cell ligs onto dotplot by fusion status. 

NKligs.act <- c("PVR", "PVRL2", "CD48", "CD58", "NCR3LG1", 
                "ICAM1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L")


NKligs.act.dotplot <- DotPlot(RMS.noMYOD1, features = NKligs.act, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)
NKligs.act.dotplot

#Rename PAX7::FOXO1 and PAX3::FOXO1 to FP.
RMS.noMYOD1$fusion <- gsub("PAX3::FOXO1", "FP", RMS.noMYOD1$fusion)
unique(RMS.noMYOD1$fusion)

RMS.noMYOD1$fusion <- gsub("PAX7::FOXO1", "FP", RMS.noMYOD1$fusion)
unique(RMS.noMYOD1$fusion)

RMS.noMYOD1$fusion <- gsub("FP", "FP-RMS", RMS.noMYOD1$fusion)
unique(RMS.noMYOD1$fusion)

#NK inhibitory ligands by FP vs FN with collagen.
Idents(RMS.noMYOD1) <- RMS.noMYOD1$fusion

NKligs.inhib.dotplot.FPFN <- DotPlot(RMS.noMYOD1, features = NKligs.inhib, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.inhib.dotplot.FPFN

#NK inhibitory ligands by FP vs FN wihout collagen.
NKligs.inhib.nocol.dotplot.FPFN <- DotPlot(RMS.noMYOD1, features = NKligs.inhib.nocol, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.inhib.nocol.dotplot.FPFN

#NK activating ligands by FP vs FN. 
NKligs.act.dotplot.FPFN <- DotPlot(RMS.noMYOD1, features = NKligs.act, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.act.dotplot.FPFN
