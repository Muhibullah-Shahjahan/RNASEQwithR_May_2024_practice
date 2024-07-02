# Required R Packages
install.packages(c('tidyverse', 'ggpubr', 'openxlsx'))

# Required BioC packages
BiocManager::install(c('GEOquery', 'TCGAbiolinks', 'DESeq2', 'airway'))

# Load packages
library(tidyverse)
library(GEOquery)
library(ggpubr)
library(openxlsx)

### Data Manipulation ###
#0. Import
data <- read.csv("Data/GSE183947_fpkm_long.csv")
data



#data exploration
dim(data)
ncol(data)
nrow(data)

head(data)
head(data, 12)

tail(data)

# sampling
sample(data)
sample_n(data, 100)
sample_frac(data, .043)

# missing value
is.na(data)
sum(is.na(data))

# specific package for finding out missing values
install.packages("naniar")
library(naniar)
miss_var_summary(data)
gg_miss_var(data)


#1. Select (used to subsetting columns)
select(data, 1)
select(data, c(1,3,4))

select(data, 1:4)

select(data, gene)
select(data, gene, title)

#2. filter

# filter data for single variable (==, >, <, >=, <=)

names(data)
filter(data, metastasis == "no")


filter(data, fpkm>10)

# using &
filter(data, metastasis == "yes" & fpkm > 10)
# using |
filter(data, metastasis == "yes" | fpkm > 10)



### select and filter
select(data, c(1,3,6)) |> 
  filter(metastasis == "yes" & fpkm == 4)

## multiple filtering criteria ( %in% )
data |> 
  filter(gene %in% c('BRCA1', 'BRCA2', 'TP53')) |> 
  head()



#3. mutate ~ changing variable ~ add new changed or operation column

data |> 
  mutate(fpkm_log = log(fpkm)) |> 
  head()



#4. grouping and summarizing
data |> 
  filter(gene == 'BRCA1' | gene == "BRCA2") |> 
  group_by(tissue) |> 
  summarise(mean(fpkm))


data |> 
  filter(gene == 'BRCA1' | gene == "BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean_fpkm = mean(fpkm)) # can change the column names


# 5. Arrange
data |> 
  filter(gene == 'BRCA1' | gene == "BRCA2") |> 
  group_by(tissue, gene) |> 
  summarise(mean_fpkm = mean(fpkm)) |>  # can change the column names
  arrange(desc(mean_fpkm)) # normally in ascending order











