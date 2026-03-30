#Loading Danielli data 
library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)
install.packages("SeuratData") #says not available for this version of R?

RMS <- readRDS("RMS_atlas_final_20240130.rds")
head(RMS@meta.data) #It is count data 
RMS
head(RMS$PatientID)

#Looking at metadata for patient ID some of the samples were taken from same patient 
unique(RMS$PatientID)
length(RMS$PatientID)

#How you can see metadata in Seurat
RMS[[]]

VlnPlot(RMS, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3)

#Looking at the labels within the metadata
unique(RMS$fusion) #Has the two different FP-RMS subtypes, FN-RMS and SS-RMS 
#(MYOD1)

#Christina code- visualization of the UMAP of the Seurat object
Idents(RMS) = RMS$fusion

RMS.fusion.umap <- DimPlot(RMS, reduction = "umap_rpca")
plot(RMS.fusion.umap)

RMS.fusion.pca <- DimPlot(RMS,reduction = "pca_rpca")

RMS@assays

#Feb 27th 2026 code- looking at NK ligands in the intergrated variable genes
#Plotting HLA-E on UMAP using FeaturePlot
HLA-E <- FeaturePlot(RMS, features = "HLA-E", reduction = "umap_rpca")
#warning message could not find HLAE in search locations found in RNA assay instead

HLA-C <- FeaturePlot(RMS, features = "HLA-C") #not found at all

#Ridge plot
feature <- c("UBE2C")
RidgePlot(RMS, features =  feature)

#Ridgeplot for HLA-E but has used the values from RNA not integrated??
RidgePlot(RMS, features = "HLA-E")

#Seeing if the NK ligands below can be plotted on UMAP (a lot of them say they are not found but found in RNA instead)
UBE2C <- FeaturePlot(RMS, features = "UBE2C")

ENTPD1 <- FeaturePlot(RMS, features = "ENTPD1")

CD48 <- FeaturePlot(RMS, features = "CD48")

CD58 <- FeaturePlot(RMS, features = "CD58")

CD274 <- FeaturePlot(RMS, features = "CD274")

NCR3LG1 <- FeaturePlot(RMS, features = "NCR3LG1")

HHLA2 <- FeaturePlot(RMS, features = "HHLA2")

HLA-DPB1 <- FeaturePlot(RMS, features = "HLA-DPB1")

#Plotting ICAM1 onto UMAP
ICAM1 <- FeaturePlot(RMS, features = "ICAM1") #no warning message
ICAM1

#Plotting ICAM1 as ridge plot
RidgePlot(RMS, features = "ICAM1")

MICA <- FeaturePlot(RMS, feature = "MICA")

MICB <- FeaturePlot(RMS, features = "MICB")

PVR <- FeaturePlot(RMS, features = "PVR")

PVRL2 <- FeaturePlot(RMS, features = "PVRL2")

ULBP1 <- FeaturePlot(RMS, features = "ULBP1")

ULBP2 <- FeaturePlot(RMS, features = "ULBP2")

ULBP3 <- FeaturePlot(RMS, features = "ULBP3")

RAET1E <- FeaturePlot(RMS, features = "RAET1E")

RAET1G <- FeaturePlot(RMS, features = "RAET1G")

RAET1L <- FeaturePlot(RMS, features = "RAET1L")

#March 2nd 2026
#Looking at NK ligands in RNA assay as they were not part of the variable genes in the data only ICAM1 was.
LayerData(RMS, assay = "integrated", layer = "data")
Layers(RMS)
RMS[["RNA"]]
RMS[["integrated"]]

RMS@reductions$pca_rpca@assay.used

RMS@active.ident

#Setting the default assay to "RNA" rather than "integrated". 
DefaultAssay(RMS) <- "RNA"

#UMAP of "RNA" using umap reduction in seurat object.
RMS@active.assay #check which assay is default now should be "RNA".

#Have to set which group you want the umap to highlight. 
Idents(RMS) <- RMS$fusion

RMS.RNA.fusion.umap <- DimPlot(RMS, reduction = "umap_rpca")
ggsave(RMS.RNA.fusion.umap, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RMS RNA fusion umap.png")

#Plotting NK ligands to the UMAP. 
HLA-E <- FeaturePlot(RMS, features = "HLA-E") #HLA not found.

MICA.RNA <- FeaturePlot(RMS, feature = "MICA")
MICA.RNA
ggsave(MICA.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/MICA.RNA.umap.png")
MICA.cells <- WhichCells(RMS, expression = MICA > 0) #4162 cells with expression over 0

#ENTPD1 also known as CD39.
ENTPD1.RNA <- FeaturePlot(RMS, features = "ENTPD1")
ENTPD1.RNA
ggsave(ENTPD1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/ENTPD1.RNA.umap 2.3.26.png")
ENTPD1.cells <- WhichCells(RMS, expression = ENTPD1 > 0) #2859

#NT5E also known as CD73.
NT5E.RNA <- FeaturePlot(RMS, features = "NT5E")
NT5E.RNA
NT5E.cells <- WhichCells(RMS, expression = NT5E > 0) #2264

CD48.RNA <- FeaturePlot(RMS, features = "CD48")
CD48.RNA
ggsave(CD48.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/CD48.RNA.umap 2.3.26.png")
CD48.cells <- WhichCells(RMS, expression = CD48 > 0) #321


CD58.RNA <- FeaturePlot(RMS, features = "CD58")
CD58.RNA
ggsave(CD58.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/CD58.RNA.umap 2.3.26.png")
CD58.cells <- WhichCells(RMS, expression = CD58 > 0) #10,283


CD247.RNA <- FeaturePlot(RMS, features = "CD247")
CD247.RNA
ggsave(CD247.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/CD247.RNA.umap 2.3.26.png")
CD247.cells <- WhichCells(RMS, expression = CD247 > 0) #2137


NCR3LG1.RNA <- FeaturePlot(RMS, features = "NCR3LG1")
NCR3LG1.RNA
ggsave(NCR3LG1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NCR3LG1.RNA.umap 2.3.26.png")
NCR3LG1.cells <- WhichCells(RMS, expression = NCR3LG1 > 0) #5046

HHLA2.RNA <- FeaturePlot(RMS, features = "HHLA2")
HHLA2.RNA
ggsave(HHLA2.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/HHLA2.RNA.umap 2.3.26.png")
HHLA2.cells <- WhichCells(RMS, expression = HHLA2 > 0) #292

HLA-C.RNA <- FeaturePlot(RMS, features = "HLA-C") #again HLA not found.

HLA-DPB1.RNA <- FeaturePlot(RMS, features = "HLA-DPB1") #again HLA not found.

ICAM1.RNA <- FeaturePlot(RMS, features = "ICAM1")
ICAM1.RNA
ggsave(ICAM1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/ICAM1.RNA.umap 2.3.26.png")
ICAM1.cells <- WhichCells(RMS, expression = ICAM1 > 0) #2,682

MICB.RNA <- FeaturePlot(RMS, features = "MICB")
MICB.RNA
plot(MICB.RNA)
ggsave(MICB.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/MICB.RNA.umap 2.3.26.png")
MICB.cells <- WhichCells(RMS, expression = MICB > 0) #942

PVR.RNA <- FeaturePlot(RMS, features = "PVR")
PVR.RNA
ggsave(PVR.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/PVR.RNA.umap 2.3.26.png")
PVR.cells <- WhichCells(RMS, expression = PVR > 0) #11,662

PVRL2.RNA <- FeaturePlot(RMS, features = "PVRL2")
PVRL2.RNA
ggsave(PVRL2.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/PVRL2.RNA.umap 2.3.26.png")
PVLR2.cells <- WhichCells(RMS, expression = PVLR2 > 0) #23,682

ULBP1.RNA <- FeaturePlot(RMS, features = "ULBP1")
ULBP1.RNA
ggsave(ULBP1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/ULBP1.RNA.umap 2.3.26.png")
ULBP1.cells <- WhichCells(RMS, expression = ULBP1 > 0) #1452

ULBP2.RNA <- FeaturePlot(RMS, features = "ULBP2")
ULBP2.RNA
ggsave(ULBP2.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/ULBP2.RNA.umap 2.3.26.png")
ULBP2.cells <- WhichCells(RMS, expression = ULBP2 > 0) #1162

ULBP3.RNA <- FeaturePlot(RMS, features = "ULBP3")
ULBP3.RNA
ggsave(ULBP3.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/ULBP3.RNA.umap 2.3.26.png")
ULBP3.cells <- WhichCells(RMS, expression = ULBP3 > 0) #1025

RAET1E.RNA <- FeaturePlot(RMS, features = "RAET1E")
RAET1E.RNA
ggsave(RAET1E.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RAET1E.RNA.umap 2.3.26.png")
RAET1E.cells <- WhichCells(RMS, expression = RAET1E > 0) #201

RAET1G.RNA <- FeaturePlot(RMS, features = "RAET1G")
RAET1G.RNA
ggsave(RAET1G.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RAET1G.RNA.umap 2.3.26.png")
RAET1G.cells <- WhichCells(RMS, expression = RAET1G > 0) #125

RAET1L.RNA <- FeaturePlot(RMS, features = "RAET1L")
RAET1L.RNA
ggsave(RAET1L.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RAET1L.RNA.umap 2.3.26.png")
RAET1L.cells <- WhichCells(RMS, expression = RAET1L > 0) #11

#More inhibitory ligands added on 24 March 2026.
#PDCD1LG2/PDL2
PDCD1LG2.RNA <- FeaturePlot(RMS, features = "PDCD1LG2")
PDCD1LG2.RNA
ggsave(PDCD1LG2.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/PDCD1LG2.RNA.umap 24.3.26.png")
PDCD1LG2.CELLS <- WhichCells(RMS, expression = PDCD1LG2 > 0) #number of cells with expression above 0 = 1110

CLEC2D.RNA <- FeaturePlot(RMS, features = "CLEC2D")
CLEC2D.RNA
ggsave(CLEC2D.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/CLEC2D.RNA.umap 24.3.26.png")
CLECD2.CELLS <- WhichCells(RMS, expression = CLECD2 > 0) #9368

LGALS9.RNA <- FeaturePlot(RMS, features = "LGALS9")
LGALS9.RNA
ggsave(LGALS9.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/LGALS9.RNA.umap 24.3.26.png")
LGALS9.CELLS <- WhichCells(RMS, expression = LGALS9 > 0) #955

#Below are the collagen genes which have been found to be upregulated/ worse prognosis.
COL18A1.RNA <- FeaturePlot(RMS, features = "COL18A1")
COL18A1.RNA
ggsave(COL18A1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/COL18A1.RNA.umap 24.3.26.png")
COL18A1.CELLS <- WhichCells(RMS, expression = COL18A1 > 0) #54315

COL4A1.RNA <- FeaturePlot(RMS, features = "COL4A1")
COL4A1.RNA
ggsave(COL4A1.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/COL4A1.RNA.umap 24.3.26.png")
COL4A1.CELLS <- WhichCells(RMS, expression = COL4A1 > 0) #47405

COL4A2.RNA <- FeaturePlot(RMS, features = "COL4A2")
COL4A2.RNA
ggsave(COL4A2.RNA, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/NK ligs Danielli RNA UMAPs 2.3.26/COL4A2.RNA.umap 24.3.26.png")
COL4A2.CELLS <- WhichCells(RMS, expression = COL4A2 > 0) #53872

#Two other collagen genes which are still expressed so not sure whether to look at all of them?
#Would it be that all or most collagen genes are expressed but not all are up or cause worse prognosis?
COL3A1.RNA <- FeaturePlot(RMS, features = "COL3A1")
COL3A1.RNA

COL2A1.RNA <- FeaturePlot(RMS, features = "COL2A1")
COL2A1.RNA


#UMAP of RNA with clusters rather than fusion type. 
Idents(RMS) <- RMS$`Cluster assignment`

RMS.RNA.CLUSTER.umap <- DimPlot(RMS, reduction = "umap_rpca")
RMS.RNA.CLUSTER.umap
ggsave(RMS.RNA.CLUSTER.umap, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RMS RNA CLUSTER umap 2.3.26.png")

#Mapping all NK ligands onto dotplot
DefaultAssay(RMS) <- "RNA"
RMS@active.assay

NKligs <- c("ENTPD1", "CD48", "CD58", "CD274", "NCR3LG1", "ICAM1", "MICA", "MICB", 
            "PVR", "PVRL2", "ULBP1", "ULBP2", "ULBP3", "RAET1E", "RAET1G", "RAET1L", "NT5E")

Idents(RMS) <- RMS$fusion

NKligs.dotplot <- DotPlot(RMS, features = NKligs) + RotatedAxis()
NKligs.dotplot
ggsave(NKligs.dotplot, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RMS NK ligs dotplot w/ NT5E/CD73 12.3.26.png")
getwd()

PVRL2.RNA <- FeaturePlot(RMS, features = "PVRL2")
PVRL2.RNA

#gives you a list of cells which have expression of PVRL2 over 0. 
PVRL2.subset <- WhichCells(RMS, expression = PVRL2 > 0)

NKligs.clusters.dotplot <- DotPlot(RMS, features = NKligs) + RotatedAxis()
NKligs.clusters.dotplot
ggsave(NKligs.clusters.dotplot, filename= "C:/Users/adpd1g23/OneDrive - University of Southampton/Documents/RPC/Danielli data/RPC Danielli/RMS NK ligs clusters dotplot 5.3.26.png")
