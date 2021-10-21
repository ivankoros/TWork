library(DESeq2)
library(tidyverse)
library(openxlsx)

#Reading data
df <- read_csv(file.choose()) %>%
  column_to_rownames("Identifier")

#DP IR vs DP Sham
dfdpirsham <- df %>%
  dplyr::select(15:22)

condition <- c(rep("DP IR", 4),
               rep("DP Sham", 4))

meta <- data.frame(SampleName = colnames(dfdpirsham),
                   Condition = condition)

dpirdpshamDDS <- DESeqDataSetFromMatrix(round(dfdpirsham),
                                        meta,
                                        ~ Condition)

#Checking normalization factors
ddsSF <- estimateSizeFactors(dpirdpshamDDS)

ddsNorm <- counts(ddsSF,
                  normalized = TRUE) %>%
  View()

#DESeq analysis and results

dpirdpshamSeq <- DESeq(dpirdpshamseq)

results <- results(dpirdpshamSeq,
                   contrast = c("Condition", "DP IR", "DP Sham") ) %>%
  as.data.frame() %>%
  rownames_to_column("Identifier")

write.xlsx(results,"my DP IR DP Sham.xlsx",
           overwrite = TRUE)
 