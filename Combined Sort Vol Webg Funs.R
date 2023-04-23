
" Sort, Volcano, WebGestaltR (SVW)
This function processes a given gene expression data file and outputs
a filtered data set, volcano plot, and gene ontology enrichment analysis.

This currently is set to work with mouse data, but can be easily modified
to work with human data by changing the species in the WebGestaltR function
and changing the gene annotation package in the getGenes function.

Inputs:
filename: The name of the file containing the gene expression data.
Adjpcol: The column number of the adjusted p-value column.
Log2fcol: The column number of the log2 fold change column.
Symbolcol: The column number of the gene symbol column.
outName: The name of the output files.

Outputs:
A filtered data set, volcano plot, and gene ontology enrichment analysis.
"
 
 # Sort
SVW <- function(filename, Adjpcol, Log2fcol, Symbolcol, outName) {

  # Load required packages/installs
  if (!require('pacman')) install.packages('pacman')
  pacman::p_load(tidyverse, openxlsx, mygene, ggrepel, WebGestaltR)
  
  # Read input data from either .xlsx or .csv file
  if((grepl("\\.xlsx$", filename)))
    DF_Raw <- read.xlsx(filename, colNames = TRUE)
  else if(grepl("\\.csv$", filename))
    DF_Raw <- read_csv(filename,
                       show_col_types = FALSE,
                       col_names = TRUE)
  
  # Set up output folder
  dir <- getwd()
  dir.create(paste0(getwd(), "/", outName))
  setwd(paste0(getwd(), "/", outName))
  dirFold <- getwd()
  
  # Process input data
  DF_Raw <- DF_Raw %>%
    dplyr::rename(Adjp=Adjpcol, 
                  Log2f=Log2fcol,
                  SymbolI=Symbolcol) %>%
    mutate(
      SaR = case_when(
        Log2f > .3 & Adjp < .05 ~ "Upregulated",
        Log2f < -.3 & Adjp < .05 ~ "Downregulated",
        TRUE ~ "NotSignificant"))
  
  # Annotate genes with MyGene
  GLabel <- getGenes(DF_Raw$SymbolI,
                     fields=c('symbol', 'name'),
                     species="mouse") %>%
    as.data.frame()
  
  # Add gene information to the input data
  DF_Named <- DF_Raw %>%
    mutate(Symbol = GLabel$symbol,
           Name = GLabel$name) %>%
    relocate(c(Symbol, Name), .before = SymbolI)

  # Filter the input data based on significance and fold change
  DF_Sorted <- DF_Named %>%
    filter(Adjp < 0.05) %>%
    filter(abs(Log2f) >= .305) %>%
    drop_na(Symbol) %>%
    filter(!str_detect(Name, "expressed sequence"),
           !str_detect(Name, "predicted gene"),
           !str_detect(Name, "cDNA"),
           !str_detect(Name, "pseudogene"),
           !str_detect(Name, "(?i)riken",)
           )
  
  # Save the filtered data as an Excel file
  write.xlsx(DF_Sorted, paste(outName, "Filtered.xlsx"))
  
  # Create a volcano plot
  ggplot(DF_Named, aes(x=Log2f, y=-log10(Adjp)), label=SaR) +
    geom_point(aes(color=SaR), size=1) +
    scale_color_manual(values=c('midnightblue', 'grey70', "darkred")) +
    theme_bw() + 
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black")) +
    xlim(-5, 5) +
    xlab("Log2 fc") +
    geom_text_repel(
      data=slice_min(DF_Named, Adjp,n=10), aes(label= Symbol)
    )
  
  # Save the volcano plot as a jpg file
  ggsave(paste(outName, "Volcano .jpg"),
         plot = last_plot(),
         device = "jpg")
  
  # Perform gene ontology enrichment analysis with WebGestaltR
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
  
  # WebGestalt outputs an important .txt file, which I convert to an Excel file
  # and save it for easier viewing  
  read_tsv(paste0(dirFold,"/Project_", outName, " WebG", "/enrichment_results_", outName, " WebG.txt")) %>%
    write.xlsx(paste0(outName, " Enrichment Results.xlsx"))
  
  # Return to the original working directory (not the output folder)
  setwd(dir)
}


# Example
SVW(filename = "Cool file name with Differentially expressed genes.xlsx",
    Adjpcol = 7,
    Log2fcol = 3,
    Symbolcol = 1,
    outName = "CD IR vs CD Sham")

