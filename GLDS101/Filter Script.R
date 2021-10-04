library(tidyverse)

GLDS101 <- read_csv("GLDS-101_rna_seq_differential_expression.csv")

##Created sorted data set named GLDS101_Sorted
GLDS101_Sorted <- GLDS101 %>%
  tibble() %>%
  filter(`Adj.p.value_(FLT)v(GC)` < 0.05) %>%
  filter(abs(`Log2fc_(FLT)v(GC)`) > .3) %>%
  drop_na(SYMBOL) %>%
  filter(!grepl("RIKEN", GENENAME)) %>%
  select(SYMBOL, `Log2fc_(FLT)v(GC)`, `Adj.p.value_(FLT)v(GC)` ) %>%
  rename(AdjpS = `Adj.p.value_(FLT)v(GC)`)
