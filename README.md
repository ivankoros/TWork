# Advanced Genomic Analysis & Visualization with R

This project provides a powerful set of R scripts designed for efficient genomic data analysis and visualization, leveraging state-of-the-art R packages to streamline and enhance your genomics workflow.

## Key Features

- Comprehensive genomic data analysis using DESeq2, WebGestalt, and Limma
- Robust data visualization techniques with ggplot2, plotPCA, and ComplexHeatmap
- Efficient data handling with tidyverse, openxlsx, and WebGestaltR

## Getting Started

These R scripts are designed to simplify the genomic analysis process. To get started, you'll need to have R and the required packages installed.

## Prerequisites

Ensure you have R and the following packages installed:

- BiocManager
- DESeq2
- tidyverse
- openxlsx
- WebGestaltR
- mygene
- ggrepel

To install all the required packages, you can use the `librarian` package. Run the following `R` code:

```R
if (!require(librarian)) {
  install.packages("librarian")
}
librarian::shelf("Bioconductor/DESeq2", "tidyverse", "openxlsx", "jokergoo/ComplexHeatmap", "Bioconductor/WebGestaltR", "BiocManager", "slowkow/ggrepel")
```
## Combined Function (I recommend this)
1. `SVW`: This function streamlines the entire workflow for analyzing genomic data. It takes a `.csv` or `.xlsx` input file, and generates filtered data, a volcano plot, and enrichment results using WebGestaltR.

    - `filename`: The input file in `.csv` or `.xlsx` format.
    - `Adjpcol`: The column name for adjusted p-values in the input file.
    - `Log2fcol`: The column name for log2 fold change values in the input file.
    - `Symbolcol`: The column name for gene symbols in the input file.
    - `outName`: The output name for the generated files.
 
 After installing the required packages, use the SVW function with your specific input parameters:
 ```R
 SVW(filename = "example_data.csv",
    Adjpcol = "adj.P.Val",
    Log2fcol = "logFC",
    Symbolcol = "Gene.symbol",
    outName = "Example_Analysis")
 ```
 This will generate the following output files in a new directory:

 - A filtered data file with significant genes (upregulated, downregulated, and not significant) based on the provided cutoff values.
 - A volcano plot visualizing the significant genes.
 - Enrichment results from WebGestaltR, which provide information about the biological processes associated with the significant genes.
 
## Individual Functions
1. **SortFun**: This function processes and sorts the genomic data.
2. **VolFun**: This function generates a volcano plot from the genomic data.
3. **WebG**: This function carries out WebGestaltR analysis on the genomic data.

## How do I use these?
First, load the necessary libraries and run the required functions in the R script. For example:
```R
source("path/to/R-script")
SortFun(filename, Adjpcol, Log2fcol, Symbolcol, outName)
VolFun(filename, Adjpcol, Log2fcol, Symbolcol, outName)
WebG(filename, Adjpcol, Log2fcol, Symbolcol, outName)
```
Replace the respective arguments with your input data file, the column names/indices for Adjusted P-value, Log2 Fold Change, and Gene Symbol, and the output file name.
## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.

## Contributing

Feel free to submit pull requests or open issues to discuss and propose changes to the code. Please adhere to the standard coding conventions and best practices for R.
