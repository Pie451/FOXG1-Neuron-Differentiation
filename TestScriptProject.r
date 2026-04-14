library(tidyr)
library(dplyr)
library(omicsdata)
library(limma)
library(EDASeq)
library(edgeR)
library(plyr)
library(corrplot)
library(ggfortify)
library(SummarizedExperiment)
library(tidyverse)


data <- readRDS("/opt/omicsdata/datasets/SRP047194.rds")																								

class(data)


# extract raw strings into data frame
raw <- data.frame(
  sample = colnames(data),
  raw_attr = as.character(data$sra.sample_attributes)
)

# split key-value pairs for attributes and attribute values
tidy_attr <- raw %>%
  separate_rows(raw_attr, sep = "\\|") %>%
  separate(raw_attr, into = c("key", "value"),
           sep = ";;", fill = "right") %>%
  mutate(key = str_trim(key),
         value = str_trim(value))

# pivot wider so each attribute becomes its own column
# modifying last attribute to specify ASD or Non-ASD and remove redundant data
sample_meta <- tidy_attr %>%
  pivot_wider(names_from = key, values_from = value) %>%
  mutate(condition = if_else(str_detect(source_name, "non-ASD"), "non-ASD", "ASD")) %>%
  dplyr::select(-source_name, -`cell type`) %>%
  dplyr::rename(days_in_vitro = `number of days in vitro in terminal differentiation conditions`) %>%
  filter(days_in_vitro != 0)

# exploratory data analysis

# sanity check on sample_meta data set
glimpse(sample_meta)

# how many samples per condition
sample_meta <- sample_meta %>%
  mutate(days_in_vitro = as.numeric(days_in_vitro))

dplyr::count(sample_meta, days_in_vitro)
dplyr::count(sample_meta, condition)
dplyr::count(sample_meta, condition, days_in_vitro)

# EDA

# sample distribution
library(ggplot2)

ggplot(sample_meta, aes(x = condition, fill = condition)) +
  geom_bar() + 
  facet_wrap(~ days_in_vitro, labeller = label_both) + 
  labs(
    title = "Sample distribution by condition and timepoint",
    x = "Condition",
    y = "Number of sample"
  ) +
  theme_minimal()

# BOXPLOT

raw_count <- assay(data)
# only samples present in sample_meta
raw_count_filtered <- raw_count[, colnames(raw_count) %in% sample_meta$sample]


boxplot(raw_count_filtered[, 1:20],
        main = "Raw counts (first 20 samples)",
        ylab = "Counts",
        las = 2,
        cex.axis = 0.7)

# log-transformed
col <- factor(sample_meta$condition)
levels(col) <- c("red", "blue")
col <- as.character(col)


boxplot(log1p(raw_count_filtered),
        main = "Log-transformed counts",
        ylab = "log1p counts",
        las = 2,
        cex.axis = 0.7,
        col = col)

legend("topright",
       legend = c("ASD", "non-ASD"),
       fill = c("red", "blue"))

# RLE
library(EDASeq)

plotRLE(raw_count_filtered,
        outline = FALSE,
        las = 2,
        ylab = "RLE",
        main = "RLE of raw counts",
        col = col,
        cex.axis = 0.5)

legend("topright",
       legend = c("ASD", "non-ASD"),
       fill = c("red", "blue"))

# RLE - check outlier samples
rle_matrix <- log1p(raw_count_filtered) - rowMedians(log1p(raw_count))
rle_medians <- colMedians(rle_matrix)

# find samples that deviate most from zero
names(rle_medians) <- colnames(raw_count_filtered)
sort(abs(rle_medians), decreasing = TRUE)


# PCA
library(ggfortify)

pca <- prcomp(t(log1p(raw_count_filtered)))
summary(pca)$importance

# by condition
autoplot(pca, data = sample_meta, colour = "condition", 
         main = "PCA colored by condition")

autoplot(pca, data = sample_meta, colour = "days_in_vitro",
         main = "PCA colored by timepoint")


# MA plot
limma::plotMA(log1p(raw_count_filtered),
              main = "MA plot",
              xlab = "A",
              ylab = "M")

# check how many genes have all zero counts
table(rowSums(raw_count_filtered) == 0)

# check how many genes have very low counts
table(rowSums(raw_count_filtered) < 10)

table(rowSums(raw_count_filtered) < 20)

# filter out lowly expressed genes
raw_count_filtered <- raw_count_filtered[rowSums(raw_count_filtered) >= 20, ]

dim(raw_count_filtered)


# Normalization



