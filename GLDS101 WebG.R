library(tidyverse)
library(WebGestaltR)


GLDS101 <- read_xlsx("GLDS-101_rna_seq_differential_expression.xlsx")

GLDS101_Sorted <- GLDS101 %>%
  tibble() %>%
  filter(`Adj.p.value_(FLT)v(GC)` < 0.05) %>%
  filter(abs(`Log2fc_(FLT)v(GC)`) >= .305) %>%
  drop_na(SYMBOL) %>%
  filter(!grepl("RIKEN", GENENAME)) %>%
  select(1:5, `Log2fc_(FLT)v(GC)`, `Adj.p.value_(FLT)v(GC)`)

Symbol <- GLDS101_Sorted$SYMBOL %>%
  as.list()

WebGTxt <- read.table("Webg920.txt")

WebGestaltR(
  enrichMethod = "ORA",
  organism = "mmusculus",
  enrichDatabase =  'geneontology_Biological_Process',
  enrichDatabaseType = "genesymbol",
  interestGeneFile = "Webg920.txt",
  interestGeneType = "genesymbol",
  referenceSet = "genome",
  isOutput = TRUE,
  fdrThr = 0.05,
  minNum = 5,
  maxNum = 2000,
  sigMethod = "fdr",
  fdrMethod = "BH",
  outputDirectory = getwd(),
  dagColor = "continuous",
  projectName = NULL,
  gseaPlotFormat = c("png", "svg"),
  setCoverNum = 10,
  hostName = "http://www.webgestalt.org/")
