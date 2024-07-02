# install required bioconductor packages
BiocManager::install(c("airway", "DESeq2", "EnhancedVolcano"), force = TRUE)

# install required packages
install.packages(c("tidyverse", "ggthemes"))


# load packages
library(tidyverse)
library(airway)
library(DESeq2)
library(ggthemes)
library(EnhancedVolcano)


# get the data
data("airway")


# counts data
counts_data  <- assay(airway)

# meta data ~ col_data
col_data <- as.data.frame(colData(airway))


# making sure the row names in "col_data" (metadata) matches to columns names in "Counts_data"
all(colnames(counts_data) %in% rownames(col_data))


# checking are they are in same order??
all(colnames(counts_data) == rownames(col_data))


# prepare the col_data
col_data <- col_data |> 
  dplyr::select(c(2,3)) |> 
  dplyr::rename(dexamethasone = "dex") |> 
  dplyr::rename(cellline = "cell") 

# clean dexamethasone column
col_data$dexamethasone <- gsub("untrt", "Untreated", col_data$dexamethasone)
col_data$dexamethasone <- gsub("trt", "Treated", col_data$dexamethasone)


















