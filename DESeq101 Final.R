library(tidyverse)
library(tximport)
library(DESeq2)
library(plotly)
library(mygene)
library(compareDF)

#Using tximport to import all 12 raw data genes.results files
raw_file_names <- list.files(path="C:/Users/theas/Desktop/R Desk/Exp1/GLDS101 Raw", full.names = TRUE)

tx101 <- tximport(files =raw_file_names,
                  type = "rsem",
                  ignoreTxVersion = TRUE,
                  txIn = TRUE)

tx101$length[tx101$length == 0] <- 1

#Attaching metadata for DESeq2 to recognize variables

setwd('C:/Users/theas/Desktop/R Desk/Exp1')

Meta <- read_csv("CustomMeta.csv",
                 col_names = TRUE) %>%
  as.data.frame()

#Carrying out DESeq2 normlizations from tximport data
deseq101 <- DESeqDataSetFromTximport(txi = tx101,
                                     colData = Meta,
                                     design = ~ Condition)

VST101 <- varianceStabilizingTransformation(deseq101)

DESeqFinal <- DESeq(deseq101)

#PCA 
ggplotly(plotPCA(VST101, intgroup = "Condition"))

## 3 steps to deseq2
# 1) estimate size factors (normalization)
# 2) estimate dispersions
# 3) apply statistics (Wald test)

EstimSF <- estimateDispersions(EstimSF)

plotDispEsts(EstimSF)

MaybThis <- DESeq(deseq101)

# step 3

EstimSF = nbinomWaldTest(EstimSF)

#results

results_table <- results(DESeqFinal)

data_df <- as.data.frame(results_table)

#Plotting 1 gene that did not pass Cook's test, possible outlier?
plotCounts(EstimSF, gene = "ENSMUSG00000079224",
           intgroup = 'Condition')

#Filter and naming
filter_df1 <- drop_na(data_df)

filter2_df1 <- filter_df1 %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 0.305)

GNames <- getGenes(row.names(filter2_df1), fields=c('symbol', 'name')) %>%
  as.data.frame()

named_df1 <- filter2_df1 %>%
  mutate(Symbol = GNames$symbol,
         Name = GNames$name) %>%
  relocate(c(Symbol, Name), .before = baseMean) %>%
  filter(!grepl("RIKEN", Name))





