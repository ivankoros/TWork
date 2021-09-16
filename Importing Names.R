library(tidyverse)
library(mygene)
library(openxlsx)

StN_Raw <- read_table("Galaxy92-[DESeq2_result_file_on_data_16,_data_14,_and_others].txt", col_names = FALSE)

names(StN_Raw) <- c('Gene ID',	'Mean Normalized Counts','Log2f',	'SE of Log2f',
               'Wald Statistic',	'PVal',	'BHP')

StN_Sorted <- StN_Raw %>%
  filter(BHP < 0.05) %>%
  filter(abs(Log2f) > .305)

GSymbol <- getGenes(StN_Sorted$`Gene ID`, fields=c('symbol', 'name')) %>%
  as.data.frame()

StN_Lookup <- mutate(StN_Sorted,
              Symbol = GSymbol$symbol,
              Name = GSymbol$name) %>%
  relocate(c(Symbol, Name), .before = `Gene ID`)

write.xlsx(StN_Lookup, "IRvSham mygene R.xlsx", overwrite = TRUE)
