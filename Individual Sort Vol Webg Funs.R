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
if (!require(WebGestaltR))
library(WebGestaltR)

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
    filter(!grepl("Rik", SymbolI))%>%
    as.data.frame()
   
  GLabel <- getGenes(DF$SymbolI, fields=c('symbol', 'name')) %>%
    as.data.frame()
  
  DF <- mutate(DF,
                       Symbol = GLabel$symbol,
                       Name = GLabel$name) %>%
    relocate(c(Symbol, Name), .before = SymbolI)
  
  write_csv(DF, paste(outName, "Filtered.csv"))
  
}

# Volcano function
VolFun <-
  function(filename,
           Adjpcol,
           Log2fcol,
           Symbolcol,
           outName) {
   
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
      scale_color_manual(values = c('blue', 'grey70', "green3")) +
      theme_bw() +
      theme(
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      xlim(-5, 5) +
      geom_label_repel(data = slice_min(DF_Raw, Adjp, n = 10), aes(label = Symbol))
    
    ggsave(paste(outName, "Volcano .pdf"),
           plot = last_plot())
    
  }

# WebG Function
WebG <- function(filename, Adjpcol, Log2fcol, Symbolcol, outName) {
  
    if((grepl("\\.xlsx$", filename)))
      DF_Raw <- read.xlsx(filename, colNames = TRUE)
    else if(grepl("\\.csv$", filename))
      DF_Raw <- read_csv(filename,
                         show_col_types = FALSE,
                         col_names = TRUE)
    
    DF_Sorted <- DF_Raw %>%
      dplyr::rename(Adjp=Adjpcol, 
             Log2f=Log2fcol,
             SymbolI=Symbolcol) %>% 
      filter(Adjp < 0.05) %>%
      filter(abs(Log2f) >= .305) %>%
      drop_na(SymbolI) %>%
      filter(!grepl("Rik", SymbolI)) %>%
      as.data.frame()
    
    GLabel <- getGenes(DF_Sorted$SymbolI, fields=c('symbol', 'name')) %>%
      as.data.frame()
    
    DF_Named <- mutate(DF_Sorted,
                     Symbol = GLabel$symbol,
                     Name = GLabel$name)
    
  Symbol <- DF_Named$Symbol %>%
    as.vector()
  
  WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase =  'geneontology_Biological_Process',
    enrichDatabaseType = "genesymbol",
    interestGene = Symbol,
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
    projectName = paste(outName, "WebG Results"),
    gseaPlotFormat = c("png", "svg"),
    setCoverNum = 10,
    hostName = "http://www.webgestalt.org/")
}
