renv::status()
renv::hydrate()
renv::snapshot(prompt = F)

library(usethis)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(Matrix.utils)

# load sgrna metadata
sgrna_tab <- read_delim("input/CRISPRko_rat_astrocyte_coculture_WTA/sgrna_tab.csv.gz", 
                        delim = "\t") %>%
  mutate(across(everything(), as.factor)) # make gene group a factor

# check number of unique gene groups
sgrna_tab$sgrna_symbol %>% unique %>% length

# load the matrix of counts
mat <- Read10X("input/CRISPRko_rat_astrocyte_coculture_WTA/matrices_filtered")
dimnames(mat)

# create a sce object
sce <- SingleCellExperiment(assays = list(counts = mat), 
                            colData = sgrna_tab)

# create a named vector of gene groups
gene_groups <- purrr::set_names(levels(sce$sgrna_symbol))

length(gene_groups)

mat_agg <- aggregate.Matrix(t(counts(sce)), 
                 groupings = colData(sce)[,c("sgrna_symbol")], 
                 fun = "sum")

mat_agg %>% as_tibble(rownames = "gene_id") %>% write_csv("output/aggregated_counts.csv")

mat_agg %>% as.matrix() %>% View
mat %>% as.matrix() %>% View

# calculate cpm
library(edgeR)

cpm(t(as.matrix(mat_agg))) %>% as_tibble(rownames = "gene_id") %>% write_csv("output/aggregated_cpm.csv")