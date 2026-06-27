library(Seurat)
library(tidyverse)
library(plyr)
BiocManager::install("DESeq2")
library(DESeq2)

#Loading Danielli data
RMS <- readRDS("RMS_atlas_final_20240130.rds")
unique(RMS$fusion)

#removed MYOD1.
RMS.noMYOD1 <- subset(RMS, subset = fusion != "MYOD1")
unique(RMS.noMYOD1$fusion)

#Remover FN-RMS.
RMS.noMYOD1.noFN <- subset(RMS.noMYOD1, subset = fusion != "FN-RMS")
unique(RMS.noMYOD1.noFN$fusion)



#make newfusion to change PAX3::FOXO1 and PAX7::FOXO1 into PAX3-FOXO1 and PAX7-FOXO1.
RMS.noMYOD1.noFN$newfusion <- mapvalues(x= RMS.noMYOD1.noFN$fusion, from= c("PAX3::FOXO1", "PAX7::FOXO1"),
                           to= c("PAX3-FOXO1", "PAX7-FOXO1"))
unique(RMS.noMYOD1.noFN$newfusion)

#Set the default assay to RNA not the integrated which is batch corrected. 
DefaultAssay(RMS.noMYOD1.noFN) <- "RNA"

unique(RMS.noMYOD1$PatientID)
unique(RMS.noMYOD1$id)

#Aggregate the counts based on the id of the sample and the newfusion status to seurat obj.
Aggcounts <- AggregateExpression(RMS.noMYOD1.noFN, assays= "RNA",
                    group.by = c("id", "newfusion"),
                    slot = "counts",
                    return.seurat = TRUE)
unique(Aggcounts$newfusion)

Idents(Aggcounts) <- "newfusion"
unique(Aggcounts$newfusion)
bulk.DE <- FindMarkers(object = Aggcounts,
                       ident.1 = "PAX3-FOXO1",
                       ident.2 = "PAX7-FOXO1",
                       test.use = "DESeq2")

head(bulk.DE, n= 15)

#Remove NA from DGE dataframe.
bulk.DE.noNA <- bulk.DE[!is.na(bulk.DE$p_val_adj), ]

#Make gene rownames into a column so can subset NK ligs out.
library(tibble)
bulk.DE.noNA <- rownames_to_column(bulk.DE.noNA, var = "Gene_ID")


#Label up and downregulated and no DE genes. 
bulk.DE.noNA <- bulk.DE.noNA %>% 
  mutate(gene_type = case_when(avg_log2FC > 0 & p_val_adj <= 0.05 ~ "Upregulated",
                               avg_log2FC < 0 & p_val_adj <= 0.05 ~ "Downregulated",
                               TRUE ~ "Not DE"))

sig.F.up <- bulk.DE.noNA %>% 
  filter(gene_type == "Upregulated")

sig.F.down <- bulk.DE.noNA %>% 
  filter(gene_type == "Downregulated")

cols.F <- c("Upregulated" = "#D41159", "Downregulated" = "#1A85FF", "Not DE"= "grey")

#NONE of the NK ligands are differentially expressed between PAX3 and PAX7 by searching.
#So the below is still copied from the FN vs FP dge so can visualize this and label the NK ligands to show that they are not DE. 

Inhib.ligs <- bulk.DE.noNA %>% 
  filter(Gene_ID == "ENTPD1" | Gene_ID == "CD274" | Gene_ID == "PDCD1LG2" | Gene_ID == "CLEC2D" | Gene_ID == "LGALS9" | Gene_ID == "NT5E" | Gene_ID == "PVR" | Gene_ID == "PVRL2")
#LGALS9 is NA so doesn't come up on the list. 

library(ggrepel)

volc.Inhib <- ggplot(bulk.DE.noNA,
       aes(x= avg_log2FC,
           y= -log10(p_val_adj))) +
  geom_point(aes(colour = gene_type)) +
  geom_point(data = sig.F.up,
             colour = "#D41159") +
  geom_point(data = sig.F.down,
             colour = "#1A85FF") +
  scale_colour_manual(values = cols.F) +
  labs(colour = "PAX3-FOXO1 vs PAX7-FOXO1") +
  geom_point(data = Inhib.ligs,
             colour = "darkblue") +
  geom_label_repel(data = Inhib.ligs,
                   aes(label = Gene_ID),
                   nudge_y = 22,
                   nudge_x = -3)
volc.Inhib
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
