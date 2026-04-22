library(Seurat)
library(tidyverse)

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
                    normalization.method = "LogNormalize")

#get the matrix from the list, it is under RNA as you set assays to RNA. 
Aggcountsmat <- Aggcounts$RNA
