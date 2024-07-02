# install packages
BiocManager::install("recount")
BiocManager::install("apeglm")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("clusterProfiler")
BiocManager::install("msigdbr")

# load packages
library(tidyverse)
library(DESeq2)
library(recount) # used to retrieve rnaseq data directly

# Problem :   Effects of EWS-FLI1 KD in Ewing sarcoma

# RQ: what is the effect of knocking down EWS-FLI1 in Ewing sarcoma? How does it change cellular behavior?

# the dataset we use here is downloaded from recount 
# find a project of interest 

project_info <- abstract_search(query = "Ewing sarcoma")

view(project_info)

# selected study

download_study("SRP015989")

# load the downloaded data
load("SRP015989/rse_gene.Rdata")

# examine the read count
assay(rse_gene)

# counts data
counts_data  <- assay(rse_gene)

# meta data ~ col_data
col_data <- as.data.frame(colData(rse_gene))

# fix the colData to give a column with the appropriate groups
rse_gene$condition <- c(rep("shCTR", 3), rep("shEF1", 4))



# making sure the row names in "col_data" (metadata) matches to columns names in "Counts_data"
all(colnames(counts_data) %in% rownames(col_data))


# checking are they are in same order??
all(colnames(counts_data) == rownames(col_data))

############ Differential gene expression analyis #####################

# step 1: construct a DESeqDataSet data object from rse_data

dds <- DESeqDataSet(
  rse_gene, design = ~condition  # condition
)
dds
plotPCA(dds)

# normalization of dataset
rld <- rlog(dds)
plotPCA(rld)  # PCA plot only drawn from normalized data

# step 2: pre-filtering (keep rows that have at least 10 reads total)

#keep <- rowSums(counts(rlg)) >= 10

# from normalized data set kept filter data again in dds
#dds <- dds[keep, ]

# set reference category !!!!!!!!! Important  !!!!!!!!!!!

#dds$dexamethasone <- relevel(dds$dexamethasone, ref = "Untreated")

# step 3: perform differential gene expression analysis
dds <- DESeq(dds)

# Step 4: Save as results

res <- results(dds)
res |> 
  head()


# Step 5: Exploring results
summary(res)


# Visualization (What? Why? How?)
# MA Plot (M vs A Plot)
# What: An MA plot represents the log-fold change (M) on the y-axis and the average expression (A) on the x-axis for each gene or feature.
# Why: It is used to visualize differential expression between two conditions.The log-fold change (M) gives an idea of the magnitude of change, and the average expression (A) helps identify if the change is dependent on the expression level.
plotMA(dds)



# Interpretations
# 1. Points on the plot represent genes.
# 2. Genes with significant differential expression are often `found at the extremes` of the plot.
# 3. Upregulated genes are at the top
# 4. Downregulated genes at the bottom
# 5. Non-differentially expressed genes are centered around zero on the y-axis.


























