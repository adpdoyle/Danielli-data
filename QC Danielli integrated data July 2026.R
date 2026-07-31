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

#Scatter comparing nFeature to percent mito
FeatureScatter(RMS, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "origin", split.by = "origin") 

FeatureScatter(RMS.fil, feature1 = "percent.mt", feature2 = "nFeature_RNA", group.by = "origin", split.by = "origin") + NoLegend()

#Scatter of nFeature and nCount.
FeatureScatter(RMS.fil, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "origin", split.by = "origin") + NoLegend() +
  scale_y_continuous(breaks = seq(0, 8000, by= 1000))
#There is still some cells below the filtering genes threshold of under 1000 for wei et al. and under 400 for Patel et al.
#Do I need to filter this out as well? And if I do, do I need to reintegrate? Check with Dean on monday!.

#Filtered the cells which were less than 1000 for wei et al and less than 400 for patel et al.
RMS.Wei.fil <- subset(RMS.fil, subset = (origin == "Wei et al." & nFeature_RNA >= 1000))

RMS.wei.pat.fil <- subset(RMS.fil, subset = (origin== "Wei et al." & nFeature_RNA >=1000) |
                            (origin == "Patel et al." & nFeature_RNA >= 400) |
                            (origin == "Weng et al.") |
                            (origin == "Danielli et al."))

dim(RMS.wei.pat.fil) #10642 cells now after filtering. 

RMS.wei.pat.fil$newfusion <- mapvalues(x= RMS$fusion, from= c("PAX3::FOXO1", "PAX7::FOXO1"),
                           to= c("FP-RMS", "FP-RMS"))

unique(RMS.wei.pat.fil.noMYOD1$newfusion)

RMS.wei.pat.fil.noMYOD1 <- subset(RMS.wei.pat.fil, subset = newfusion != "MYOD1")

DefaultAssay(RMS.wei.pat.fil.noMYOD1) <- "RNA"

#Aggregate the counts based on the id of the sample and the newfusion status to seurat obj.
Aggcounts <- AggregateExpression(RMS.wei.pat.fil.noMYOD1, assays= "RNA",
                                 group.by = c("id", "newfusion"),
                                 slot = "counts",
                                 return.seurat = TRUE)
Aggcounts$newfusion

Idents(Aggcounts) <- "newfusion"

bulk.DE.FNFP <- FindMarkers(object = Aggcounts,
                            ident.1 = "FN-RMS",
                            ident.2 = "FP-RMS",
                            test.use = "DESeq2")

bulk.DE.FNFP.noNA <- bulk.DE.FNFP[!is.na(bulk.DE.FNFP$p_val_adj), ]

#Make gene rownames into a column so can subset NK ligs out.
library(tibble)
bulk.DE.FNFP.noNA <- rownames_to_column(bulk.DE.FNFP.noNA, var = "Gene_ID")

bulk.DE.FNFP.noNA <- bulk.DE.FNFP.noNA %>% 
  mutate(gene_type = case_when(avg_log2FC > 0 & p_val_adj <= 0.05 ~ "Upregulated",
                               avg_log2FC < 0 & p_val_adj <= 0.05 ~ "Downregulated",
                               TRUE ~ "Not significant"))

#Redid the DGE and has not changed anything!!!. 


FeatureScatter(RMS.wei.pat.fil, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "origin", split.by = "origin") + NoLegend() +
  scale_y_continuous(breaks = seq(0, 8000, by= 1000))

wei.before.fil <- sum(RMS.fil$origin == "Wei et al.")
wei.after.fil <- sum(RMS.Wei.fil$origin == "Wei et al.")
unique(RMS.Wei.fil$origin)
rm(RMS.Wei.fil)

# Count Wei cells that fall outside the thresholds
wei_removed_count <- sum(
  RMS.fil$origin == "Wei et al." &
    (RMS.fil$nFeature_RNA < 1000 |
       RMS.fil$nFeature_RNA > 8000))

wei.before.fil - wei_removed_count
#Number of cells per dataset.
table(RMS$origin)

RMS$origin %>%
  table() %>%
  as.data.frame() %>%
  ggplot(aes(x =., y = Freq)) +
  geom_col(fill = "steelblue") +
  theme_classic() +
  labs(x = "Dataset", y = "Number of Cells") +
  scale_y_continuous(breaks = seq(0, 60000, by= 5000))

# Number of cells per sample 72. 
sort(table(RMS$name))

# Make a table of cell counts per sample 72.
count_table <- as.data.frame(table(RMS$name))
colnames(count_table) <- c("Sample", "Cell Count")

# Write to CSV
write.csv(count_table, "qc_cell_counts.csv", row.names = FALSE)


# Quick check: what's their mito like?
RMS$below_wei_threshold <- RMS$origin == "Wei et al." & RMS$nFeature_RNA < 1000
VlnPlot(RMS, features = "percent.mt", group.by = "below_wei_threshold", pt.size = 0)

RMS.fil$below_wei_threshold <- RMS.fil$origin == "Wei et al." & RMS.fil$nFeature_RNA < 1000
VlnPlot(RMS.fil, features = "percent.mt", group.by = "below_wei_threshold", pt.size = 0)

VlnPlot(RMS.fil, features = "nCount_RNA", group.by = "below_wei_threshold", pt.size=0)
unique(RMS$PatientID)
RMS$name

# See what slots exist in the RNA assay
GetAssayData(RMS, assay = "RNA", layer = "counts")[1:5, 1:5]
GetAssayData(RMS, assay = "RNA", layer = "data")[1:5, 1:5]
GetAssayData(RMS, assay = "RNA", layer =  "scale.data")[1:5, 1:5]
#the counts layer has integers whereas data has decimals which means counts is raw and data is normalized.
#AI says that when you set RNA as default it pulls from the normalized data not the counts. 