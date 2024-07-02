# Load packages
library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)
library(naniar)

# import data raw counts data
counts_data <- read.csv("Data/GSE183947.csv")


# meta data ~~~~~~~~~ !! Important
res <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
res
class(res)

# meta data
metaData <- pData(phenoData(res[[1]]))

metaData

names(metaData)

# subset of meteData
metaData_sub <- metaData |> 
  select(c(1,10,11,17))

# preprocessing of data
metaData_modified <- metaData_sub |> 
  rename(tissue = characteristics_ch1, metastasis = characteristics_ch1.1) |> 
  mutate(tissue = gsub("tissue: ", "", tissue)) |> 
  mutate(metastasis = gsub("metastasis: ", "", metastasis))
metaData_modified


# reshaping data
counts_data_long <- counts_data |> 
  rename( gene = X) |> # reshape the data
  pivot_longer(-gene, 
               names_to = "samples",
               values_to = "fpkm")
  
# joining Data
counts_final <- counts_data_long |> 
  left_join(metaData_modified, by = c("samples" = "description"))


# Export data
write.csv(counts_final, "Data/GSE183947_counts.csv", row.names = F)


