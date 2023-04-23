if (!require('pacman')) install.packages('pacman')
pacman::p_load(AnnotationDbi, readxlsx, tidyverse, mygene, org.Mm.eg.db, clusterProfiler)

DP_Sham_CD_Sham = read.xlsx("DP Sham vs CD Sham Filtered.xlsx")

geneSet = DP_Sham_CD_Sham$Symbol[DP_Sham_CD_Sham$Log2f]

# I use the org.Mm.eg.db package to convert mouse gene symbols to entrez IDs
# This is necessary for the enrichGO function

GO_results = enrichGO(
  gene = geneSet,
  OrgDb="org.Mm.eg.db",
  key="SYMBOL",
  ont="ALL",
  pAdjustMethod = "BH"
  )

# Plotting GO results
dotplot(GO_results,
        x="GeneRatio",
        showCategory=10,
        split=NULL)

# Upregualted and downregulated DP Sham v CD Sham
DPSham_df = read.xlsx("DP Sham vs CD Sham Filtered.xlsx")

DPSham_Up = data.frame(DPSham_df) %>%
  filter(Log2f>0) %>%
  select('Symbol')

DPSham_Down = data.frame(DPSham_df) %>%
  filter(Log2f<0) %>%
  select('Symbol')



