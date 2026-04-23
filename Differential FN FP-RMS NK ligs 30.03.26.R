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

#Aggregate the counts based on the id of the sample and the newfusion status.
Aggcounts <- AggregateExpression(RMS.noMYOD1, assays= "RNA",
                    group.by = c("id", "newfusion"),
                    slot = "counts",
                    return.seurat = TRUE)

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

#Volcano plot of bulk.DE FP vs FN.
library(ggplot2)
ggplot(bulk.DE, aes(x= avg_log2FC, y= -log10(p_val_adj))) + geom_point()

bulk.DE.noNA$diffexp[bulk.DE.noNA$avg_log2FC > 0 &
                bulk.DE.noNA$p_val_adj < 0.05] <- "upregulated"

bulk.DE.noNA$diffexp[bulk.DE.noNA$avg_log2FC < 0 &
                  bulk.DE.noNA$p_val_adj < 0.05] <- "downregulated"

bulk.DE.noNA$diffexp[bulk.DE.noNA$p_val_adj > 0.05] <- "not DE"

volc1 <- ggplot(bulk.DE.noNA, aes(x= avg_log2FC, y= -log10(p_val_adj), colour = diffexp)) + 
  geom_point(aes(colour = diffexp)) +
  scale_colour_manual(values= mycolours) 
volc1

volc2 <- ggplot(bulk.DE.noNA, aes(x= avg_log2FC, y= -log10(p_val_adj), colour = diffexp)) + 
  geom_point(aes(colour = diffexp)) +
  scale_colour_manual(values= mycolours) +
  geom_point(data = NK.ligs,
             colour = "black") +
  geom_label_repel(data = NK.ligs, aes(label = Gene_ID))

volc2

mycolours <- c("#D41159", "#1A85FF", "grey")

 
names(mycolours) <- c("upregulated", "downregulated", "not DE")
volc1


library(tibble)
library(ggrepel)
bulk.DE <- rownames_to_column(bulk.DE, var = "Gene_ID")
NK.ligs <- bulk.DE %>% 
  filter(Gene_ID == "ENTPD1")
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
