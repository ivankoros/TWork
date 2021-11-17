library(tidyverse)
library(gplots)
library(cluster)
library(factoextra)
library(ComplexHeatmap)

install.packages("tidyverse")
install.packages("gplots")
install.packages("cluster")
install.packages("factoextra")

rm(list = ls())
setwd("C:/Users/Ivan Korostenskij/Desktop/R/TWork")

df <- read_csv("porgchamp.csv") %>%
  column_to_rownames(1) %>%
  data.matrix()

log_df <- log10(df + 1)

sc_df <- scale(log_df, scale = TRUE)

#avg
Conc <- data.frame("CD Baseline" = rowMeans(df[,c(1:4, 23)]),
                   "CD IR" = rowMeans(df[,5:7]),
                   "CD Sham" = rowMeans(df[,8:11]),
                   "DP Baseline" = rowMeans(df[,12:14]),
                   "DP IR" = rowMeans(df[,15:18]),
                   "DP Sham" = rowMeans(df[,19:22]))

avg_df <- Conc %>%
  data.matrix()

avg_log_df <- log10(avg_df + 1) %>%
  scale(scale = TRUE)

fviz_nbclust(avg_log_df, kmeans, method = "silhouette")

final <- kmeans(avg_log_df, 2, nstart = 25)
fviz_cluster(final, data = avg_log_df, labelsize = 1)

heatmap(avg_log_df,
        cexCol = 0.9,
        scale = "row",
        col = terrain.colors(256))


#combining 
list <- c(17523,
          102502,
          16590,
          56193,
          72310,
          57349,
          80718,
          225825,
          16416,
          330409,
          218442,
          18476,
          22364,
          77045,
          545486,
          16184,
          16197,
          328780,
          211550,
          68794,
          17754,
          15507,
          12373,
          270210,
          68273,
          16145,
          216831,
          20315,
          83453
) %>%
  as.character()

usedf <- df %>%
  as.data.frame()

select_df <- avg_df[list,] %>%
  as.matrix()

#scale before vcisualization
select_df <- log10(select_df + 1)

data_scaled <- t(scale(t(select_df)))

pallete <- colorRampPalette(c("blue",
                              "white",
                              "red"))

heatmap.2(data_scaled,
          cexCol = 0.9,
          col = pallete(n = 20),
          trace = "none",
          labRow = )



rowlabels <- c("Mpo",
               "Pls1",
               "Kit",
               "Plek",
               "Nkg7",
               "Ppbp",
               "Rab27b",
               "Cd226",
               "Itgb3",
               "Mpl",
               "Cecr2",
               "Serinc5",
               "Pafah1b3",
               "Vpreb3",
               "Bcl7a",
               "Tubb1",
               "Il2ra",
               "Il7r",
               "Prss34",
               "Tifa",
               "Flnc",
               "Map1a",
               "Hspb1",
               "Casq2",
               "Zfp651",
               "Pomgnt1",
               "Igtp",
               "Arhgap44",
               "Cxcl12",
               "Chrdl1"
)

heatmap.2(data_scaled,
          cexCol = 0.9,
          col = pallete(n = 20),
          trace = "none",
          labRow = rowlabels,
          dendrogram = "row")

install.packages("pheatmap")
library(pheatmap)

pheatmap(select_df)
         ,
         display_numbers = TRUE,
         scale = "row")

,
         labels_row = rowlabels)

#filename = "pheatmap.pdf",








