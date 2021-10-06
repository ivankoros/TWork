library(DESeq2)
library(tidyverse)
library(openxlsx)
library(mygene)
library(ggrepel)
library(WebGestaltR)
library(edgeR)
library(apeglm)

DP_IR_DP_Sham <- read_tsv("Galaxy103-[Normalized_counts_file_on_data_15,_data_13,_and_others].tabular") %>%
  column_to_rownames('...1')

Condition <- factor(c(rep("DP_Sham", 4), rep("DP_IR", 4)),
                    levels = c("DP_IR", "DP_Sham"))

DP_IR_DP_Sham_Meta <- data.frame(SampleName = colnames(DP_IR_DP_Sham), Condition)

DP_IR_DP_Sham_Load <- DESeqDataSetFromMatrix(round(DP_IR_DP_Sham),
                                             DP_IR_DP_Sham_Meta,
                                             ~Condition)

DP_IR_DP_Sham_DESeq <- DESeq(DP_IR_DP_Sham_Load)
# relevel(DP_IR_DP_Sham_DESeq$Condition, ref = "DP_Sham")

DP_IR_DP_Sham_Results <- results(DP_IR_DP_Sham_Apeglm,
                                 alpha = .05,
                                 contrast = c("Condition", "DP_IR", "DP_Sham"),
                                 pAdjustMethod = "BH",) %>% 
  as.data.frame() %>%
  rownames_to_column("GeneID") %>%
  write_csv("2testdel.csv")

contrast =c("Condition", "DP_IR", "DP_Sham")
DP_IR_DP_Sham_Apeglm <- lfcShrink(DP_IR_DP_Sham_DESeq,
                                   contrast = c("Condition", "DP_IR", "DP_Sham"),
                                   type = "ashr")

DP_IR_DP_Sham_Apeglm <- DP_IR_DP_Sham_Apeglm %>%
  as.data.frame() %>%
  write_csv("2testdel.csv") %>%
  rownames_to_column("")

VolFun(filename = "2testdel.csv",
        Adjpcol = 7,
        Log2fcol = 3,
        Symbolcol = 1,
        outName = "2testdel")





testnotag <- DP_IR_DP_Sham_Results %>%
  as.data.frame() %>%
  filter(padj < .05,
         abs(log2FoldChange) >= 0.305) %>%
  drop_na(padj)

  

testagainst <- read.xlsx("DP IR vs DP Sham DEGs.xlsx") %>%
  filter(`P-value.Benjamini` < .05,
         abs(Log.2.fold.change) >= .305) %>%
  drop_na(`P-value.Benjamini`)







# Sort function
SortFun <- function(filename, Adjpcol, Log2fcol, Symbolcol, outName) {
  
  if((grepl("\\.xlsx$", filename)))
    DF <- read.xlsx(filename, colNames = TRUE)
  else if(grepl("\\.csv$", filename))
    DF <- read_csv(filename,
                   show_col_types = FALSE,
                   col_names = TRUE)
  
  DF <- DF %>%
    dplyr::rename(Adjp=Adjpcol, 
                  Log2f=Log2fcol,
                  SymbolI=Symbolcol) %>%
    mutate(
      SaR = case_when(
        Log2f > .3 & Adjp < .05 ~ "Upregulated",
        Log2f < -.3 & Adjp < .05 ~ "Downregulated",
        TRUE ~ "NotSignificant"
      )) %>%
    filter(Adjp < 0.05) %>%
    filter(abs(Log2f) >= .305) %>%
    drop_na(SymbolI) %>%
    as.data.frame()
  
  GLabel <- getGenes(DF$SymbolI, fields=c('symbol', 'name')) %>%
    as.data.frame()
  
  DF <- mutate(DF,
               Symbol = GLabel$symbol,
               Name = GLabel$name) %>%
    relocate(c(Symbol, Name), .before = SymbolI)
  
  DF <- DF %>%
    filter(!grepl("Rik", Symbol))
  
  write_csv(DF, paste(outName, "Filtered.csv"))
  
}

VolFun <- function(filename, Adjpcol, Log2fcol, Symbolcol, outName) {
  
  if((grepl("\\.xlsx$", filename)))
    DF_Raw <- read.xlsx(filename, colNames = TRUE)
  else if(grepl("\\.csv$", filename))
    DF_Raw <- read_csv(filename,
                       show_col_types = FALSE,
                       col_names = TRUE)
  
  DF_Raw <- DF_Raw %>%
    dplyr::rename(Adjp=Adjpcol, 
                  Log2f=Log2fcol,
                  SymbolI=Symbolcol) %>%
    mutate(
      SaR = case_when(
        Log2f > .3 & Adjp < .05 ~ "Upregulated",
        Log2f < -.3 & Adjp < .05 ~ "Downregulated",
        TRUE ~ "NotSignificant"
      ))
  
  GLabel <- getGenes(DF_Raw$SymbolI, fields=c('symbol', 'name')) %>%
    as.data.frame()
  
  DF_Raw <- mutate(DF_Raw,
                   Symbol = GLabel$symbol,
                   Name = GLabel$name)
  
  ggplot(DF_Raw, aes(x=Log2f, y=-log10(Adjp)), label=SaR) +
    geom_point(aes(color=SaR), size=1) +
    scale_color_manual(values=c('blue', 'grey70', "green3")) +
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    xlim(-5, 5) +
    geom_label_repel(
      data=slice_min(DF_Raw, Adjp,n=10), aes(label= Symbol))
  
  ggsave(paste(outName, "Volcano .pdf"),
         plot = last_plot())
}
