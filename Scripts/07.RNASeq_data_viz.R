# load required packages
library(tidyverse)
library(airway)
library(DESeq2)

library(RColorBrewer)
library(EnhancedVolcano)
library(pheatmap)

# packages for gene filtering
library(genefilter)


# package for poison distribution calculation
library(PoiClaClu)
BiocManager::install("org.Hs.eg.db")

# package for gene annotation
library(org.Hs.eg.db)
library(AnnotationDbi)

# load data
data("airway")

counts_data <- assay(airway)
meta_data <-as.data.frame(colData(airway))

# making sure the sample are equal and in same order

colnames(counts_data)
rownames(meta_data)

all(colnames(counts_data) %in% rownames(meta_data))

all(colnames(counts_data) == rownames(meta_data))


# construction a DESeqDataSetFromMatrix data object
dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = meta_data,
  design = ~dex
)

# plot heatmap of poison distribution between samples
# poison distribution for raw (non-normalizaed) count data
# Use Euclidean distance for data normalized by regularized-logarithm (rlog) or variance stabilization transformation (vst)
counts(dds)
t(counts(dds))

poisd <- PoissonDistance(t(counts(dds)))
sample_poison_matrix <- as.matrix(poisd$dd)
rownames(sample_poison_matrix) <- paste(dds$dex, sep = "-")
colnames(sample_poison_matrix) <- NULL
color <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
  sample_poison_matrix,
  clustering_distance_rows = poisd$dd,
  clustering_distance_cols = poisd$dd,
  color = color
)

# perform differential gene expression analysis
dds <- DESeq(dds)

# obtain result
res <- results(dds)
names(res)
res

# plot average expression versus log2 fold change - points are colored by

plotMA(res)
plotMA(res, ylim = c(-4,4))


# plot histogram of P-values
hist(res$pvalue, breaks = 20, col = "#411443", border = "yellow")

# plot histogram of P-adj values
hist(res$padj, breaks = 20, col = "#A84245", border = "white")

# plot histogram of P-values - improved versoion by filtering out genes with very low expression levels

hist(res$pvalue[res$baseMean > 1], col = "#A84245", border = "yellow")


# add gene annotation to results

anno <- AnnotationDbi::select(org.Hs.eg.db, 
                              rownames(res), 
                              columns=c("ENSEMBL", "ENTREZID", "SYMBOL", "GENENAME"), 
                              keytype="ENSEMBL")

# join annotation table to results join
#2 methods:

anno_res <- cbind(ENSEMBL = rownames(res), res)

# convert res into dataframe
anno_res <- as.data.frame(anno_res)

names(anno_res)


annotation_results <- left_join(anno_res, anno, by = "ENSEMBL" )
view(annotation_results)


# volcano plot
EnhancedVolcano(annotation_results,
                lab = annotation_results$SYMBOL,
                x = "log2FoldChange",
                y = "padj")

## add custom log2FC and adjusted P-value cutoff and size of point sand
EnhancedVolcano(
  annotation_results,
  lab = annotation_results$SYMBOL,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 10e-4,
  FCcutoff = 2,
  pointSize = 1.5,
  labSize = 3.0,
  title = "Untreated vs. Treated"
)


# Adjust  axis limit
EnhancedVolcano(
  annotation_results,
  lab = annotation_results$SYMBOL,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 10e-4,
  FCcutoff = 2,
  pointSize = 1.5,
  labSize = 3.0,
  xlim = c(-5,5),
  ylim = c(0, -log10(10e-10)),
  title = "Untreated vs. Treated"
)


# Modify border and remove grid lines
EnhancedVolcano(
  annotation_results,
  lab = annotation_results$SYMBOL,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 10e-4,
  FCcutoff = 2,
  pointSize = 1.5,
  labSize = 3.0,
  xlim = c(-5,5),
  ylim = c(0, -log10(10e-10)),
  border = "full",
  borderWidth = 1.5,
  borderColour = "yellow",
  gridlines.major = FALSE,
  title = "Untreated vs. Treated"
)

# perform regularized-logarithm transformation (rlog) on the data
rld <- rlog(dds)

# plot  principal component analysis
plotPCA(rld, intgroup = "dex")

# subset genes ~ significant genes

sig_result <- annotation_results[which(annotation_results$padj < 0.01) & abs(annotation_results$log2FoldChange) >= 1 & annotation_results$baseMean >= 20, ]

# DE genes with strongest downregulation
head(sig_result[order(sig_result$log2FoldChange), ])

# DE genes with strongest up-regulation
tail(sig_result[order(sig_result$log2FoldChange), ])

# plot expression of individual genes
# gene with largest positive log2FC
plotCounts(dds, gene = which.max(annotation_results$log2FoldChange), intgroup = "dex")
plotCounts(dds, gene = which.min(annotation_results$log2FoldChange), intgroup = "dex")

# specific gene of interest
plotCounts(dds, gene = "ENSG00000127954", intgroup = "dex")
