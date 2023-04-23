# Loading packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(tidyverse, openxlsx, DESeq2, pheatmap, mygene, RColorBrewer, VennDiagram, readxl)



# Venn diagram

setwd("C:/Users/Ivan Korostenskij/Desktop/R/Dried Plum/Data")

# Read in data
CD_IR_v_CD_Sham <- read_excel("CD IR vs CD Sham Filtered.xlsx")
DP_IR_V_DP_Sham <- read_excel("DP IR vs DP Sham Filtered.xlsx")
DP_Sham_V_CD_Sham <- read_excel("DP Sham vs CD Sham Filtered.xlsx")

# Create Venn diagram of DEGs
venn.diagram(
  x = list(
    CD_IR_v_CD_Sham$Symbol,
    DP_IR_V_DP_Sham$Symbol,
    DP_Sham_V_CD_Sham$Symbol
  ),
  category.names = c("CD IR vs CD Sham", "DP IR vs DP Sham", "DP Sham vs CD Sham"),
  filename = "Venn_Diagram.png",
  output = TRUE,
  lwd = 2,
  lty = 'solid',
  fill = c("#1f78b4", "#33a02c", "#e31a1c"),
  cat.cex = 0.7,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  imagetype = "png"
)

# Make a heatmap of DEGs with pheatmap
normcounts <- read_csv("normalized counts.csv") %>%
  column_to_rownames("Identifier")

GLabel <- getGenes(rownames(normcounts),
                   fields = c('symbol'),
                   species = "mouse") %>%
  as.data.frame()

row.names(normcounts) <- NULL

# Adding gene symbols to normalized counts
normcounts <- normcounts %>%
  as.data.frame() %>%
  mutate(Symbol = make.unique(GLabel$symbol)) %>%
  drop_na() %>%
  column_to_rownames('Symbol')

# Comparing normcounts symbols to filtered datasets
cdircdsham <- read.xlsx("CD IR vs CD Sham Filtered.xlsx") %>%
  dplyr::select(Symbol)

dpirdpsham <- read.xlsx("DP IR vs DP Sham Filtered.xlsx") %>%
  dplyr::select(Symbol)

dpshamcdsham <- read.xlsx("DP Sham vs CD Sham Filtered.xlsx") %>%
  dplyr::select(Symbol)

uniqueDEGs <- full_join(cdircdsham,
                        dpirdpsham) %>%
  full_join(., dpshamcdsham)

res <- as.character(uniqueDEGs$Symbol)

# Creating average dataframe to use for heatmap
average_df <-
  data.frame(
    "CD Baseline" = rowMeans(normcounts[, c(1:4, 23)]),
    "DP Baseline" = rowMeans(normcounts[, 12:14]),
    "CD Sham" = rowMeans(normcounts[, 8:11]),
    "CD IR" = rowMeans(normcounts[, 5:7]),
    "DP Sham" = rowMeans(normcounts[, 19:22]),
    "DP IR" = rowMeans(normcounts[, 15:18])
  ) %>%
  data.matrix()

average_df <- average_df[res, ] %>%
  as.data.frame()

# Removing all zeros 
row_sub = apply(average_df, 1, function(row) all(row < 100 ))

average_df = average_df[row_sub,]

average_df %>% 
  as.data.frame() %>%
  rowwise() %>% 
  filter(any(c_across(everything(.)) == 0))

average_df <- rownames_to_column(average_df)

# Save average_df to excel
write.xlsx(average_df,"filterour100counts.xlsx")

# Running pheatmap
cleanColnames <- c("CD Baseline", "DP Baseline", "CD Sham",
                   "CD IR", "DP Sham", "DP IR")

# Make heatmap of DEGs (but with average counts)
pheatmap(
  average_df,
  display_numbers = FALSE,
  scale = "row",
  width = 10,
  height = 10,
  show_rownames = FALSE,
  treeheight_row = 0,
  angle_col = 0,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  labels_col = cleanColnames,
  fontsize_row = 1,
  filename = "phm.pdf"
)
