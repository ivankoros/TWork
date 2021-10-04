##Load packages
library(tidyverse)
library(BiocManager)
library(plotly)
library(ggrepel)

GLDS101 <- read_csv(file.choose())

##Created sorted data set named GLDS101_Sorted
GLDS101_Sorted <- GLDS_101_rna_seq_differential_expression %>%
  tibble() %>%
  filter(`Adj.p.value_(FLT)v(GC)` < 0.05) %>%
  filter(abs(`Log2fc_(FLT)v(GC)`) > .3) %>%
  drop_na(SYMBOL) %>%
  filter(!grepl("RIKEN", GENENAME)) %>%
  select(SYMBOL, `Log2fc_(FLT)v(GC)`, `Adj.p.value_(FLT)v(GC)` ) %>%
  rename(AdjpS = `Adj.p.value_(FLT)v(GC)`)

## Renaming columns
GLDS101 <- GLDS101 %>% rename(Adjp=`Adj.p.value_(FLT)v(GC)`, 
                              Log2f=`Log2fc_(FLT)v(GC)`)

## Creating new column to distinguish between down, up, and notsig
GLDS101 <- GLDS101 %>%
  mutate(
    SaR = case_when(
      Log2f > .3 & Adjp < .05 ~ "Upregulated",
      Log2f < -.3 & Adjp < .05 ~ "Downregulated",
      TRUE ~ "NotSignificant"
    ))

GLDS101 <- GLDS101 %>% mutate(-log10(Adjp))
N
options(ggrepel.max.overlaps = Inf)

## Creating volcano plot named Vol101
Vol101 <- ggplot(GLDS101, aes(x=Log2f, y=-log10(Adjp)), label=SaR) +
  geom_point(aes(color=SaR), size=1) +
  scale_color_manual(values=c('blue', 'grey70', "green3")) +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
  xlim(-5, 5) +
  geom_label_repel(
    data=slice_min(GLDS101, Adjp,n=10), aes(label=SYMBOL))

##Turning ggplot Vol101 into a plotly
ggplotly(Vol101)

## Quick save option for Vol101
ggsave(Vol101, file="Vol1.pdf")

