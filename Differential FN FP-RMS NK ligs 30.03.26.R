library(dplyr)
library(Seurat)
library(patchwork)
library(SeuratData)
library(ggplot2)
library(plyr)

#Loading Danielli data
RMS <- readRDS("RMS_atlas_final_20240130.rds")
unique(RMS$fusion)

#Renaming the PAX3::FOXO1 and PACX7::FOXO1 to have FP-RMS in newfusion metadata column.
RMS$newfusion <- mapvalues(x= RMS$fusion, from= c("PAX3::FOXO1", "PAX7::FOXO1"),
                           to= c("FP-RMS", "FP-RMS"))
unique(RMS$newfusion)

RMS.noMYOD1 <- subset(RMS, subset = newfusion != "MYOD1")
unique(RMS.noMYOD1$newfusion)
length(RMS$newfusion) #107523
length(RMS.noMYOD1$newfusion) #106023

#Check that the same number of MYOD1 has been removed from above
RMS.MYOD1 <- subset(RMS, subset = newfusion == "MYOD1")
length(RMS.MYOD1$newfusion) #1500 so the above did remove all the MYOD1. 
