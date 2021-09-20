library(tidyverse)
library(WebGestaltR)

GLDS101 <- read_csv("GLDS-101_rna_seq_differential_expression.csv")

GLDS101_Sorted <- GLDS101 %>%
  tibble() %>%
  filter(`Adj.p.value_(FLT)v(GC)` < 0.05) %>%
  filter(abs(`Log2fc_(FLT)v(GC)`) >= .305) %>%
  drop_na(SYMBOL) %>%
  filter(!grepl("RIKEN", GENENAME)) %>%
  select(1:7, `Log2fc_(FLT)v(GC)`, `Adj.p.value_(FLT)v(GC)`)

Ensemble <- GLDS101_Sorted$ENSEMBL %>%
  as.data.frame()

WebGestaltR(
  enrichMethod = "ORA",
  organism = "mmusculus",
  enrichDatabase =  'geneontology_Biological_Process',
  enrichDatabaseType = "ensembl_gene_id",
  interestGeneFile = "GTxt.txt",
  interestGeneType = "ensembl_gene_id",
  collapseMethod = "mean",
  referenceSet = "genome",
  isOutput = TRUE,
  outputDirectory = getwd(),
  projectName = NULL,
  gseaPlotFormat = c("png", "svg"),
  setCoverNum = 10,
  hostName = "http://www.webgestalt.org/")
