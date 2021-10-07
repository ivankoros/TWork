library(tidyverse)
library(openxlsx)
library(mygene)
library(ggrepel)
library(WebGestaltR)

# Sort
SVW <- function(filename, Adjpcol, Log2fcol, Symbolcol, outName) {
  
  if((grepl("\\.xlsx$", filename)))
    DF_Raw <- read.xlsx(filename, colNames = TRUE)
  else if(grepl("\\.csv$", filename))
    DF_Raw <- read_csv(filename,
                       show_col_types = FALSE,
                       col_names = TRUE)
  
  dir <- getwd()
  dir.create(paste0(getwd(), "/", outName))
  setwd(paste0(getwd(), "/", outName))
  dirFold <- getwd()
  
  DF_Raw <- DF_Raw %>%
    dplyr::rename(Adjp=Adjpcol, 
                  Log2f=Log2fcol,
                  SymbolI=Symbolcol) %>%
    mutate(
      SaR = case_when(
        Log2f > .3 & Adjp < .05 ~ "Upregulated",
        Log2f < -.3 & Adjp < .05 ~ "Downregulated",
        TRUE ~ "NotSignificant"))
  
  GLabel <- getGenes(DF_Raw$SymbolI, fields=c('symbol', 'name')) %>%
    as.data.frame()
  
  DF_Named <- DF_Raw %>%
    mutate(Symbol = GLabel$symbol,
           Name = GLabel$name) %>%
    relocate(c(Symbol, Name), .before = SymbolI)
  
  DF_Sorted <- DF_Named %>%
    filter(Adjp < 0.05) %>%
    filter(abs(Log2f) >= .305) %>%
    drop_na(SymbolI) %>%
    filter(!grepl("Rik", Symbol))%>%
    as.data.frame()
  
  write.xlsx(DF_Sorted, paste(outName, "Filtered.xlsx"))
  
  # Vol
  ggplot(DF_Named, aes(x=Log2f, y=-log10(Adjp)), label=SaR) +
    geom_point(aes(color=SaR), size=1) +
    scale_color_manual(values=c('blue', 'grey70', "green3")) +
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    xlim(-5, 5) +
    geom_label_repel(
      data=slice_min(DF_Named, Adjp,n=10), aes(label= Symbol))
  
  ggsave(paste(outName, "Volcano .pdf"),
         plot = last_plot())
  
  # WebG
  Symbol <- DF_Sorted$Symbol %>%
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
    projectName = paste(outName, "WebG"),
    gseaPlotFormat = c("png", "svg"),
    setCoverNum = 10,
    hostName = "http://www.webgestalt.org/")
  
  read_tsv(paste0(dirFold,"/Project_", outName, " WebG", "/enrichment_results_", outName, " WebG.txt")) %>%
    write.xlsx(paste0(outName, " Enrichment Results.xlsx"))
  
  setwd(dir)
}
