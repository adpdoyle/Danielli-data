library(Seurat)
library(tidyverse)
library(plyr)
BiocManager::install("DESeq2")
library(DESeq2)

#Loading Danielli data
RMS <- readRDS("RMS_atlas_final_20240130.rds")
unique(RMS$fusion)

#Renaming the PAX3::FOXO1 and PACX7::FOXO1 to have FP-RMS in newfusion metadata column.
RMS$newfusion <- mapvalues(x= RMS$fusion, from= c("PAX3::FOXO1", "PAX7::FOXO1"),
                           to= c("FP-RMS", "FP-RMS"))
unique(RMS$newfusion)
packageVersion("plyr")

#removed MYOD1.
RMS.noMYOD1 <- subset(RMS, subset = newfusion != "MYOD1")
unique(RMS.noMYOD1$newfusion)
length(RMS$newfusion) #107523
length(RMS.noMYOD1$newfusion) #106023

#Check that the same number of MYOD1 has been removed from above
RMS.MYOD1 <- subset(RMS, subset = newfusion == "MYOD1")
length(RMS.MYOD1$newfusion) #1500 so the above did remove all the MYOD1. 

#Set the default assay to RNA not the integrated which is batch corrected. 
DefaultAssay(RMS.noMYOD1) <- "RNA"

unique(RMS.noMYOD1$PatientID)
unique(RMS.noMYOD1$id)

#Aggregate the counts based on the id of the sample and the newfusion status to seurat obj.
Aggcounts <- AggregateExpression(RMS.noMYOD1, assays= "RNA",
                    group.by = c("id", "newfusion"),
                    slot = "counts",
                    return.seurat = TRUE)

#Aggregate the counts based on the id of the sample and the newfusion status to non-seurat obj.
Aggcounts.mat <- AggregateExpression(RMS.noMYOD1, assays= "RNA",
                                 group.by = c("id", "newfusion"),
                                 slot = "counts",
                                 return.seurat = FALSE)
#22 Apr 2026
Aggcounts$idnewfusion <- paste(Aggcounts$id, Aggcounts$newfusion, sep = "-")
unique(Aggcounts$id)

Idents(Aggcounts) <- "newfusion"
unique(Aggcounts$newfusion)
bulk.DE <- FindMarkers(object = Aggcounts,
                       ident.1 = "FP-RMS",
                       ident.2 = "FN-RMS",
                       test.use = "DESeq2")

head(bulk.DE, n= 15)

bulk.DE.FNFP <- FindMarkers(object = Aggcounts,
                       ident.1 = "FN-RMS",
                       ident.2 = "FP-RMS",
                       test.use = "DESeq2")

#Remove NA from DGE dataframe.
bulk.DE.noNA <- bulk.DE[!is.na(bulk.DE$p_val_adj), ]

bulk.DE.FNFP.noNA <- bulk.DE.FNFP[!is.na(bulk.DE.FNFP$p_val_adj), ]
#Make gene rownames into a column so can subset NK ligs out.
library(tibble)
bulk.DE.noNA <- rownames_to_column(bulk.DE.noNA, var = "Gene_ID")
bulk.DE.FNFP.noNA <- rownames_to_column(bulk.DE.FNFP.noNA, var = "Gene_ID")

#Label up and downregulated and no DE genes. 
bulk.DE.noNA <- bulk.DE.noNA %>% 
  mutate(gene_type = case_when(avg_log2FC > 0 & p_val_adj <= 0.05 ~ "Upregulated",
                               avg_log2FC < 0 & p_val_adj <= 0.05 ~ "Downregulated",
                               TRUE ~ "Not DE"))

bulk.DE.FNFP.noNA <- bulk.DE.FNFP.noNA %>% 
  mutate(gene_type = case_when(avg_log2FC > 0 & p_val_adj <= 0.05 ~ "Upregulated",
                               avg_log2FC < 0 & p_val_adj <= 0.05 ~ "Downregulated",
                               TRUE ~ "Not significant"))

sig.F.up <- bulk.DE.noNA %>% 
  filter(gene_type == "Upregulated")

sig.F.up <- bulk.DE.FNFP.noNA %>% 
  filter(gene_type == "Upregulated")

sig.F.down <- bulk.DE.noNA %>% 
  filter(gene_type == "Downregulated")
sig.F.down <- bulk.DE.FNFP.noNA %>% 
  filter(gene_type == "Downregulated")

NS <- bulk.DE.FNFP.noNA %>% 
  filter(gene_type == "Not significant")

cols.F <- c("Upregulated" = "#D41159", "Downregulated" = "#1A85FF", "Not significant"= "grey")

#Volcano plot for inhibitory NK ligands. 26 June 2026. 
Inhib.ligs <- bulk.DE.FNFP.noNA %>% 
  filter(Gene_ID == "ENTPD1" | Gene_ID == "CD274" | Gene_ID == "PDCD1LG2" | Gene_ID == "CLEC2D" | Gene_ID == "LGALS9" | Gene_ID == "NT5E" | Gene_ID == "PVR" | Gene_ID == "PVRL2")
#NT5E NA in DGE so has been taken out of the bulk.DE.FNFP.noNA.

library(ggrepel)

volc.ligs.inhib <- ggplot(bulk.DE.FNFP.noNA,
                            aes(x= avg_log2FC,
                                y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "FN- vs FP-RMS") +
  geom_point(data = Inhib.ligs,
             colour = "darkblue") +
  geom_label_repel(data = Inhib.ligs,
                   aes(label = Gene_ID),
                   nudge_y = 42,
                   nudge_x = -1,
                   box.padding = 0.5)
volc.ligs.inhib

#Volcano plot for activating ligands. 26 June 2026.
Act.ligs <- bulk.DE.FNFP.noNA %>% 
  filter(Gene_ID == "CD48" | Gene_ID == "CD58" | Gene_ID == "NCR3LG1" | Gene_ID == "ICAM1" | Gene_ID == "MICA" | Gene_ID == "MICB" | Gene_ID == "PVR" | Gene_ID == "PVRL2"|
           Gene_ID == "ULBP1" | Gene_ID == "ULBP2" | Gene_ID == "ULBP3" | Gene_ID == "RAET1E" | Gene_ID == "RAET1G" | Gene_ID == "RAET1L")

volc.ligs.act <- ggplot(bulk.DE.FNFP.noNA,
                          aes(x= avg_log2FC,
                              y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "FN- vs FP-RMS") +
  geom_point(data = Act.ligs,
             colour = "darkblue") +
  geom_label_repel(data = Act.ligs,
                   aes(label = Gene_ID),
                   nudge_y = 62,
                   nudge_x = 0,
                   box.padding = 0.7,
                   max.overlaps = Inf)

volc.ligs.act


ICAM1.ENTPD1 <- bulk.DE.noNA %>% 
  filter(Gene_ID == "ICAM1" | Gene_ID == "ENTPD1")
View(ICAM1.ENTPD1)

library(ggrepel)

volc.ICAM1.ENTPD1 <- ggplot(bulk.DE.noNA,
       aes(x= avg_log2FC,
           y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "FP vs FN") +
  geom_point(data = ICAM1.ENTPD1,
             colour = "darkblue") +
  geom_label_repel(data = ICAM1.ENTPD1,
                   aes(label = Gene_ID),
                   nudge_y = 22,
                   nudge_x = 1)

#volcano plot with top 10 differentially expressed genes. 
top10DE <- head(bulk.DE.noNA, n=10)

volc.top10DE <- ggplot(bulk.DE.noNA,
       aes(x= avg_log2FC,
           y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "FP vs FN") +
  geom_point(data = top10DE,
             colour = "black") +
  geom_label_repel(data = top10DE,
                   aes(label = Gene_ID))

volc.top10DE

#Volcano plot with top 5 up and down reg genes.
top5.up <- bulk.DE.noNA %>% 
  filter(avg_log2FC > 0) %>% 
  arrange(p_val_adj) %>% 
  head(n= 5)

top5.down <- bulk.DE.noNA %>% 
  filter(avg_log2FC < 0) %>% 
  arrange(p_val_adj) %>% 
  head(n= 5)

volc.top5.up.down <- ggplot(bulk.DE.noNA,
                            aes(x= avg_log2FC,
                                y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "FP vs FN") +
  geom_point(data = top5.up,
             colour = "red") +
  geom_point(data = top5.down,
             colour= "darkblue") +
  geom_label_repel(data = top5.up,
                   aes(label = Gene_ID),
                   nudge_x = 1.5,
                   nudge_y = 10) +
  geom_label_repel(data = top5.down,
                   aes(label = Gene_ID),
                   nudge_x = -1,
                   nudge_y = 7)


volc.top5.up.down

#there are 3767 rows with NA, I think this is due to those genes having 0 values which is common in scRNA data. 

#when I looked in the bulk.DE dataframe all the NK ligs apart from ICAM1 and ENTPD1 were non-significant.
#I tried to sense check and ChatGPT said between FP and FN fusion that the FOXO1 and PAX3 and 7 are upregulated but these weren't in bulkDE.
#So don't think that the DE is accurate?.

#Below is making the Aggcounts a dataframe then transposing to get Gene ID as columns and sample ID as rows which is then a large matrix.
#From Aggcount.mat.t can start following that tutorial again. 
Aggcounts.mat.mat <- as.data.frame(Aggcounts.mat)
rm(Aggcounts.mat.mat)

geneID <- row.names(Aggcounts.mat.mat)

Aggcounts.mat.t <- t(Aggcounts.mat.mat)
