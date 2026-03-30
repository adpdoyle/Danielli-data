library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)

#Loading Danielli data
RMS <- readRDS("RMS_atlas_final_20240130.rds")

#Mapping all NK ligands (apart from HLA) onto Dotplot grouping by fusion status. 
DefaultAssay(RMS) <- "RNA"
RMS@active.assay

NKligs <- c("ENTPD1", "CD48", "CD58", "CD274", "NCR3LG1", "ICAM1", "MICA", "MICB", 
            "PVR", "PVRL2", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L", "NT5E")

Idents(RMS) <- RMS$fusion

NKligs.dotplot <- DotPlot(RMS, features = NKligs) + RotatedAxis()
NKligs.dotplot
ggsave(NKligs.dotplot, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RMS NK ligs dotplot w/ NT5E/CD73 12.3.26.png")
getwd()

#Mapping just inhibitory NK cell ligs onto dotplot by fusion status. 
DefaultAssay(RMS) <- "RNA"
RMS@active.assay


NKligs.inhib <- c("CD274","PDCD1LG2","CLEC2D","LGALS9","COL18A1","COL4A1","COL4A2","ENTPD1", "NT5E","PVR", "PVRL2")

NKligs.inhib.nocol <- c("CD274","PDCD1LG2","CLEC2D","LGALS9","ENTPD1", "NT5E","PVR", "PVRL2")

Idents(RMS) <- RMS$fusion

NKligs.inhib.dotplot <- DotPlot(RMS, features = NKligs.inhib, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

NKligs.inhib.dotplot
ggsave(NKligs.inhib.dotplot, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK inhib ligs with collagen 24.3.26.png")

#Nk inhib ligs without the collagen genes. 
NKligs.inhib.nocol.dotplot <- DotPlot(RMS, features = NKligs.inhib.nocol, scale.min = 0, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features") +
  scale_colour_viridis(option = "viridis", direction = -1)

#Mapping just activating NK cell ligs onto dotplot by fusion status. 
DefaultAssay(RMS) <- "RNA"
RMS@active.assay

NKligs.act <- c("PVR", "PVRL2", "CD48", "CD58", "NCR3LG1", 
                "ICAM1", "MICA", "MICB", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L")

Idents(RMS) <- RMS$fusion

RMS@active.ident <- factor(RMS@active.ident,
                           levels = c("PAX3::FOXO1", "PAX7::FOXO1", "FN-RMS", "MYOD1"))

NKligs.act.dotplot <- DotPlot(RMS, features = NKligs.act, scale.by = "size", scale = FALSE) + RotatedAxis() +
  labs(y = "Fusion status", x= "Features")
NKligs.act.dotplot

#Mapping NK cell ligands by Inhib then activating dotplot by fusion status. 
DefaultAssay(RMS) <- "RNA"
RMS@active.assay

NKligs.inhib.act.ordered <- c("PDCD1LG2","CLEC2D","LGALS9","COL18A1","COL4A1","COL4A2","ENTPD1", "NT5E", "CD274","PVR", "PVRL2", "CD48", "CD58", "NCR3LG1", "ICAM1", "MICA", "MICB",
                      "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L")

Idents(RMS) <- RMS$fusion
RMS@active.ident <- factor(RMS@active.ident,
                           levels = c("PAX3::FOXO1", "PAX7::FOXO1", "FN-RMS", "MYOD1"))

#Changing the colour of the dotplot scale.
NKligs.dotplot.ordered <- DotPlot(RMS, features = NKligs.inhib.act.ordered, scale = FALSE, scale.by = "size") + RotatedAxis() +
  labs(y = "Fusion status", x= "") +
  scale_colour_gradient2(low = "beige", mid = "deeppink", high = "darkorchid4")

#Changing dotplot colour scale to Viridis which is colourblind friendly. 
install.packages("viridis")
library(viridis)

NKligs.dotplot.ordered <- DotPlot(RMS, features = NKligs.inhib.act.ordered, scale = FALSE, scale.by = "size") + RotatedAxis() +
  labs(y = "Fusion status", x= "") +
  scale_colour_viridis(option = "viridis", direction = -1)
NKligs.dotplot.ordered

#connecting github to this git repository. 
library(usethis) 
create_github_token()
library(gitcreds)
gitcreds_set()
use_github()




