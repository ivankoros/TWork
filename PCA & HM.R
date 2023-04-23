if (!require("pacman")) install.packages("pacman")
pacman::p_load(DESeq2, openxlsx, tidyverse, plotly, ggplot2)

cat("Your current working directory is", getwd())

# Download and choose "Raw Counts Matrix.csv" for below
rawCountsMatrix <- read_csv(file.choose()) %>%
  column_to_rownames("Identifier")

condition <- c(
  rep("Baseline", 4),
  rep("Ex 1", 3),
  rep("Ex 2", 4),
  rep("Baseline 2", 4),
  rep("Ex 3", 4),
  rep("Ex 4", 4)
)

# Creating metadata for DESeq data set creation
metaData <- data.frame(SampleName = colnames(rawCountsMatrix),
                   Condition = condition)

DFforDeseq <- DESeqDataSetFromMatrix(round(rawCountsMatrix),
                                   metaData,
                                   ~ Condition)
# DESeq analysis and results
DESeqDF <- DESeq(DFforDeseq,
                 fitType = "parametric")

results <- results(DESeqDF,
                   tidy = T,
                   independentFiltering = T,
                   alpha = 0.05) %>%
  rename("Identifier" = "row")

write.xlsx(results[order(results$padj),], "Ivan Results.xlsx",
           overwrite = T)

# Plotting PCA
DESeqVST <- varianceStabilizingTransformation(DESeqDF)

pca_results <- plotPCA(DESeqVST, intgroup = "Condition", returnData = TRUE)

pca_results$SampleName <- metaData$SampleName

ggplot(pca_results, aes(PC1, PC2, color = Condition)) +
  geom_point() +
  geom_text(aes(label = SampleName), hjust = 0, vjust = 0, size = 3) +
  scale_color_discrete(name = "Condition") +
  xlab("PC1") +
  ylab("PC2") +
  ggtitle("PCA Plot with Sample Names")






# DESeq to create CD IR v DP Sham

sample_condition =  condition[c(5:7, 19:23)]

metadata_sample_condition <- data.frame(SampleName = colnames(rawCountsMatrix[c(5:7, 19:23)]),
                       Condition = sample_condition)

dds_df <- DESeqDataSetFromMatrix(round(rawCountsMatrix[c(5:7, 19:23)]),
                                          metadata_sample_condition,
                                     ~ Condition)

e# DESeq analysis and results
deseq_out_df <- DESeq(dds_df,
                 fitType = "parametric")

results <- results(deseq_out_df,
                   tidy = T,
                   independentFiltering = T,
                   alpha = 0.05) %>%
  rename("Identifier" = "row")

write.xlsx(results[order(results$padj),], "Results.xlsx",
           overwrite = T)
