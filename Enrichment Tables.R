if (!require("pacman")) {
  install.packages("pacman")
}
pacman::p_load(tidyverse, openxlsx, WebGestaltR, mygene)

setwd("C:/Users/Ivan Korostenskij/Desktop/R/Dried Plum/Data")

## Making reference set

# Load data
ref_set_id = read_csv("All normalized counts.csv") %>%
  data.frame()

# Removing non-detected (0 count) genes

ref_set_id_non_zero = ref_set_id %>%
  mutate(sumVar = rowSums(across(CD_Baseline_D02:DP__Baseline_D22))) %>%
  select(Identifier, sumVar) %>%
  filter(sumVar != 0)
  
# Gene labeling 
ref_set_get_genes = queryMany(ref_set_id_non_zero$Identifier,
                          fields=c("symbol", "entrezgene"),
                          species="mouse",
                          return.as="DataFrame")

# Removing NA values and outputting results to 'reference_set.txt'
ref_set_entrez = data.frame(ref_set_get_genes$entrezgene) %>%
  drop_na() %>%
  write_delim("reference set.txt",
              col_names = FALSE)


## Making individual upregulated/downregulated files for each condition
# DP Sham v CD Sham data
DPSham_df = read.xlsx("DP Sham vs CD Sham Filtered.xlsx")

DPSham_Up = data.frame(DPSham_df) %>%
  filter(Log2f>0) %>%
  dplyr::select('Symbol') %>%
  write_delim("DP Sham Upregulated Gene Set.txt",
              col_names = FALSE)

DPSham_Down = data.frame(DPSham_df) %>%
  filter(Log2f<0) %>%
  dplyr::select('Symbol') %>%
  write_delim("DP Sham Downregulated Gene Set.txt",
              col_names = FALSE)

# DP IR v DP Sham data
DPIR_df = read.xlsx("DP IR vs DP Sham Filtered.xlsx") %>%
  data.frame()

DPIR_df %>%
  filter(Log2f>0) %>%
  dplyr::select('Symbol') %>%
  write_delim("DP IR Upregulated Gene Set.txt",
              col_names = FALSE)

DPIR_df %>%
  filter(Log2f<0) %>%
  dplyr::select('Symbol') %>%
  write_delim("DP IR Downregulated Gene Set.txt",
              col_names = FALSE)


# CD IR v CD Sham data

CDIR_df = read.xlsx("CD IR vs CD Sham Filtered.xlsx") %>%
  data.frame()

CDIR_df %>%
  filter(Log2f>0) %>%
  dplyr::select('Symbol') %>%
  write_delim("CD IR Upregulated Gene Set.txt",
              col_names = FALSE)

CDIR_df %>%
  filter(Log2f<0) %>%
  dplyr::select('Symbol') %>%
  write_delim("CD IR Downregulated Gene Set.txt",
              col_names = FALSE)

## Top changes in log fold ranked by absolutel value in CD IR vs CD Sham

top_10_CDIR = CDIR_df[order(abs(CDIR_df$Log2f), decreasing = TRUE), ]
write.xlsx(top_10_CDIR, "top_CDIR.xlsx")
view(top_10_CDIR)


## function to make data frame of similarities frame between two data frames

shared_regulation = function(df1, df2, include_mixed=FALSE) {
  
  shared_genes = inner_join(df1, df2, by = "Symbol")
  
  shared_genes = shared_genes %>%
    mutate(
      regulation = case_when(
        Log2f.x > 0 & Log2f.y > 0 ~ "both positive",
        Log2f.x < 0 & Log2f.y < 0 ~ "both negative",
        TRUE ~ "mixed"
      )
    ) %>%
    dplyr::rename(Log2f_dataset1 = Log2f.x, Log2f_dataset2 = Log2f.y) %>% 
    dplyr::select(Symbol, Log2f_dataset1, Log2f_dataset2, regulation)
  
  if (include_mixed) {
    shared_genes = shared_genes %>% filter(regulation != "mixed")
  }
  
  return(shared_genes)
 
}

DP_Sham_df = read.xlsx("DP Sham vs CD Sham Filtered.xlsx")
CD_IR_df = read.xlsx("CD IR vs CD Sham Filtered.xlsx")
DP_IR_df = read.xlsx("DP IR vs DP Sham Filtered.xlsx")

DP_Sham_CD_IR_shared = shared_regulation(DP_Sham_df, CD_IR_df) %>%
  dplyr::arrange(regulation)
DP_IR_CD_IR_shared = shared_regulation(DP_IR_df, CD_IR_df)%>%
  dplyr::arrange(regulation)
DP_IR_DP_Sham_shared = shared_regulation(DP_IR_df, DP_Sham_df)

write.xlsx(DP_IR_CD_IR_shared, "DP IR CD IR shared.xlsx")
write.xlsx(DP_Sham_CD_IR_shared, "DP Sham CD IR Shared.xlsx")

table(DP_IR_DP_Sham_shared$regulation)
table(DP_IR_CD_IR_shared$regulation)
table(DP_Sham_CD_IR_shared$regulation)

library(plotly)

DP_IR_DP_Sham_shared <- DP_IR_DP_Sham_shared %>%
  rename("DP_IR" = Log2f_dataset1, "DP_Sham" = Log2f_dataset2)


ggplot(DP_IR_DP_Sham_shared, aes(x = DP_IR, y = DP_Sham, color = factor(DP_Sham > 0))) + 
  geom_point(size = 3) +
  scale_color_manual(values=c("blue", "red"), labels=c("DP Sham downregulated", "DP Sham upregulated")) +
  labs(title = "Scatter Plot of DP IR vs DP Sham",
       x = "Log2f of DP IR",
       y = "Log2f of DP Sham") +
  theme_bw() +
  guides(color = guide_legend(title = NULL))


### testing
# Read the table into a tibbl
# Count the number of genes in which Log2f_dataset1 > Log2f_dataset2
count <- sum(DP_IR_CD_IR_shared$Log2f_dataset1 > DP_IR_CD_IR_shared$Log2f_dataset2)

# Calculate the percentage
total <- sum(DP_IR_CD_IR_shared$Log2f_dataset1 > 0)

percentage <- count / total * 100

cat(sprintf("%.2f%% of upregulated genes in dataset1 have a higher Log2f fold-change than dataset2.", percentage))

###
# Count the number of genes in which Log2f_dataset1 < Log2f_dataset2
count <- sum(DP_IR_CD_IR_shared$Log2f_dataset1 < DP_IR_CD_IR_shared$Log2f_dataset2)

# Count the total number of genes being compared
total <- sum(DP_IR_CD_IR_shared$Log2f_dataset1 > 0) + sum(DP_IR_CD_IR_shared$Log2f_dataset1 < 0)

# Calculate the percentage
percentage <- count / total * 100

cat(sprintf("%.2f%% of downregulated genes in dataset1 have a lower Log2f fold-change than dataset2.", percentage))



# Rbio API
install.packages("remotes")
remotes::install_github("moosa-r/rbioapi")

view(rba_enrichr_libs())

# DSigDB
# HDSigDB_Mouse_2021
# RNAseq_Automatic_GEO_Signatures_Mouse_Up
# Mouse_Gene_Atlas
# KOMP2_Mouse_Phenotypes_2022

enrich_results = rba_enrichr(
  gene_list = DPSham_Up$Symbol,
  gene_set_library = "KOMP2_Mouse_Phenotypes_2022",
  organism = "human",
  progress_bar = TRUE
)

view(enrich_results)


# DP Sham v CD Sham data###########
DPSham_df = read.xlsx("DP Sham vs CD Sham Filtered.xlsx")

DPSham_Up = data.frame(DPSham_df) %>%
  filter(Log2f>0) %>%
  select('Symbol')

DPSham_up_entrez = queryMany(list(DPSham_df$Symbol),
                             scopes="symbol",
                             fields="entrezgene",
                             species="mouse",
                             return.as="DataFrame") %>%
  data.frame()

DPSham_Down = data.frame(DPSham_df) %>%
  filter(Log2f<0) %>%
  select('Symbol')

################################

cd_ir_looking = CD_IR_df %>%
  dplyr::select(Symbol, Name, Log2f) %>%
  arrange(desc(Log2f))

view(cd_ir_looking)
