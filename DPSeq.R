# Checking package requirements
if (!require(BiocManager))
  install.packages('BiocManger')
if (!require(DESeq2))
  BiocManager::install('DESeq2')
library(DESeq2)
if (!require(tidyverse))
  install.packages('tidyverse')
library(tidyverse)
if (!require(openxlsx))
  install.packages('openxlsx')
library(openxlsx)

"
This is a script for running differential expression analysis on a raw counts matrix.
The script will generate a volcano plot and a results matrix.

The script requires the following packages:
DESeq2
tidyverse
openxlsx

The script requires the following files:
Raw Counts Matrix.csv (or .xlsx) - A raw counts matrix with the first column as the identifier column.
Column selection of chosen conditions (see below)
"


# Reading data
df <- read_csv("Raw Counts Matrix.csv") %>%
  column_to_rownames("Identifier")

# Selecting specifc condition columns
comparison <- df %>%
  dplyr::select(15:22)

condition <- c(rep("Cond 1", 4),
               rep("Cond 2", 4))

meta <- data.frame(SampleName = colnames(comparison),
                   Condition = condition)

comparisonDDS <- DESeqDataSetFromMatrix(round(comparison),
                                        meta,
                                        ~ Condition)

# Checking normalization factors
ddsSF <- estimateSizeFactors(comparisonDDS)

ddsNorm <- counts(ddsSF,
                  normalized = TRUE) %>%
  View()

# DESeq analysis and results
dpirdpshamSeq <- DESeq(comparisonDDS)

results <- results(dpirdpshamSeq,
                   contrast = c("Condition", "Cond 1", "Cond 2")) %>%
  as.data.frame() %>%
  rownames_to_column("Identifier")

write.xlsx(results, "My Results Matrix.xlsx",
           overwrite = TRUE)

# Checking volcano plot (function below)
VolFun(
  filename = "My Results Matrix.xlsx",
  Adjpcol = 7,
  Log2fcol = 3,
  Symbolcol = 1,
  outName = "My"
)

# Volcano function
VolFun <-
  function(filename,
           Adjpcol,
           Log2fcol,
           Symbolcol,
           outName) {

    # Checking package requirements
    if (!require(tidyverse))
      install.packages('tidyverse')
    library(tidyverse)
    if (!require(openxlsx))
      install.packages('openxlsx')
    library(openxlsx)
    if (!require(mygene))
      BiocManager::install('mygene')
    library(mygene)
    if (!require(ggrepel))
      install.packages('ggrepel')
    library(ggrepel)
    
    if ((grepl("\\.xlsx$", filename)))
      DF_Raw <- read.xlsx(filename, colNames = TRUE)
    else if (grepl("\\.csv$", filename))
      DF_Raw <- read_csv(filename,
                         show_col_types = FALSE,
                         col_names = TRUE)
    else{
      stop("ERROR: .csv and .xlsx are only acceptable file formats")
    }
    
    DF_Raw <- DF_Raw %>%
      dplyr::rename(Adjp = Adjpcol,
                    Log2f = Log2fcol,
                    SymbolI = Symbolcol) %>%
      mutate(
        SaR = case_when(
          Log2f > .3 & Adjp < .05 ~ "Upregulated",
          Log2f < -.3 & Adjp < .05 ~ "Downregulated",
          TRUE ~ "NotSignificant"
        )
      )
    
    GLabel <-
      getGenes(DF_Raw$SymbolI, fields = c('symbol', 'name')) %>%
      as.data.frame()
    
    DF_Raw <- mutate(DF_Raw,
                     Symbol = GLabel$symbol,
                     Name = GLabel$name)
    
    ggplot(DF_Raw, aes(x = Log2f, y = -log10(Adjp)), label = SaR) +
      geom_point(aes(color = SaR), size = 1) +
      scale_color_manual(values = c('midnightblue', 'grey70', "darkred")) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      xlim(-5, 5) +
      xlab("Log2 fc") +
      geom_text_repel(data = slice_min(DF_Raw, Adjp, n = 10), aes(label = Symbol))
    
    ggsave(paste(outName, "Volcano .jpg"),
           plot = last_plot(),
           device = "jpg")
    
  }
