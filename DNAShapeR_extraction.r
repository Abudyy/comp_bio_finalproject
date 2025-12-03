# Install libraries
install.packages("DNAshapeR")
install.packages("caret")
install.packages("dplyr")
install.packages("tidyr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")


# Load libraries
library(DNAshapeR)
library(Biostrings)
library(dplyr)
library(tidyr)

library(Biostrings)

fa <- readDNAStringSet("C:/Users/saleh/561_project/data/CTCF_all.fa")
length(fa)

labels <- read.csv("C:/Users/saleh/561_project/data/CTCF_labels.csv")
nrow(labels)


sum(duplicated(names(fa)))     # MUST be 0
sum(duplicated(labels$id))     # MUST be 0




#EXTRACT DNA SHAPES:

library(DNAshapeR)
library(Biostrings)
library(dplyr)
library(tidyr)

fasta_file <- "C:/Users/saleh/561_project/data/CTCF_all.fa"
label_file <- "C:/Users/saleh/561_project/data/CTCF_labels.csv"

# Load sequences
fa <- readDNAStringSet(fasta_file)
ids <- names(fa)

# Shape extraction
shape <- getShape(fasta_file)

# Convert shapes → dataframe
shape_to_df <- function(shape_list) {
  dfs <- list()
  for (name in names(shape_list)) {
    mat <- shape_list[[name]]
    if (!is.null(dim(mat))) {
      df <- as.data.frame(mat)
      colnames(df) <- paste0(name, "_", seq_len(ncol(df)))
      dfs[[name]] <- df
    }
  }
  do.call(cbind, dfs)
}

shape_df <- shape_to_df(shape)

# Attach ID column
shape_df$id <- ids

# Load labels
labels <- read.csv(label_file)

# Merge
full_df <- merge(shape_df, labels, by="id")

# Save final dataset
write.csv(full_df,
          "C:/Users/saleh/561_project/data/CTCF_shapes_new.csv",
          row.names = FALSE)


