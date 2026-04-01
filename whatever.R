library(tidyr)
library(tidyverse)
library(readr)

data <- readRDS("/opt/omicsdata/datasets/SRP047194.rds")

# coldata_df <- tibble(raw = data$sra.sample_attributes) %>%
#   mutate(id = row_number()) %>%             # optional: keep track of original rows
#   separate_rows(raw, sep = "\\|") %>%       # split key-value pairs into rows
#   separate(raw, into = c("key", "value"), sep = ";;") %>%  # split key and value
#   pivot_wider(names_from = key, values_from = value) %>%   # make wide format
#   select(-id)

data$sra.sample_attributes 

coldata_list <- lapply(data$sra.sample_attributes, function(x) {
  kv <- strsplit(strsplit(x, "\\|")[[1]], ";;")
  keys <- sapply(kv, `[`, 1)
  values <- sapply(kv, `[`, 2)
  df <- as.data.frame(t(values), stringsAsFactors = FALSE)
  colnames(df) <- keys
  df
})
coldata_df <- do.call(rbind, coldata_list)

# [1] "cell type;;iPSC-derived neural progenitors" "donor id;;1123-03" ...

str(as.data.frame(data$sra.sample_attributes))
separate_wider_delim(data = as.data.frame(data$sra.sample_attributes), cols = colnames(data$sra.sample_attributes), delim = "|", names = c('cell type', 'donor ID', 'DIV', 'source name', 'day', 'Control'))                 

# table(rowSums(assay(data, "raw_counts")) >= 15)

# data$sra.sample_attributes
# glimpse(data$sra.sample_attributes)
# as.data.frame(data$sra.sample_attributes)
