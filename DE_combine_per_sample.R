# Combining differential expression analysis results from individual pegRNAs
# to a joint data frame

args <- commandArgs(trailingOnly = TRUE)


setwd(args[1])
intermediate_results_folder <- args[2]
meta_data_file <- args[3]
id <- args[4]

library(scPrimeR)
library(tidyverse)
library(scran)
library(scater)
library(scuttle)
library(stringr)


de_files <- lapply(list.files(file.path(intermediate_results_folder, "temp_DE", id)),function(x) file.path(intermediate_results_folder,"temp_DE", id, x))

peg_names <- str_match(unlist(de_files), "DE_result_[^_]+_(.*)_gene_uncorrected_temp\\.rds")[,2]

# Read individual DE files and standardise colnames (remove ":")
de_list <- lapply(de_files, function(f) {
  df <- readRDS(f)
  colnames(df) <- gsub(":", "", colnames(df))
  df
})

                           
markers_gene <- do.call(rbind,de_list)
                   
markers_pathways <- readRDS(file.path(intermediate_results_folder,
                   paste0("DE_result_",id, "_condition_pathways_uncorrected_temp.rds")))

# Remove intercept columns
markers_gene    <- markers_gene[, !grepl("Intercept", names(markers_gene))]
markers_pathways <- markers_pathways[, !grepl("Intercept", names(markers_pathways))]

# Add pegRNA meta data
meta_data <- read_csv(meta_data_file) %>%
            select(id,gene,mutation) %>% distinct()

meta_data <- meta_data %>%
  mutate(target_gene = gene) %>%
  select(-gene)

markers_gene <- markers_gene %>%
  mutate(DE_gene = gene) %>%
  select(-gene)

markers_gene <- left_join(markers_gene, meta_data, by = c("peg" = "id"))

markers_pathways <- markers_pathways %>%
  mutate(DE_gene = gene) %>% 
  select(-gene)

markers_pathways <- left_join(markers_pathways, meta_data, by = c("peg" = "id"))

# Remove cases of singular convergence
markers_pathways <- markers_pathways %>% filter(!(is_singular))
                   
markers_gene <-  markers_gene %>%  select(-gene_id)         

               
# Control FDR within each peg (Benjamini Bogomolov), as many tests for different pegs are expected
# to have strong positive corrrelations.                   
                   
add_fdr_columns <- function(markers_df, p_pattern = "^p_", group_col = "peg") {
  p_cols <- grep(p_pattern, names(markers_df), value = TRUE)
  
  if (length(p_cols) == 0) {
    warning("No p-value columns found.")
    return(markers_df)
  }
  
  markers_df %>%
    group_by(across(all_of(group_col))) %>%
    mutate(
      across(
        all_of(p_cols),
        ~ p.adjust(.x, method = "BH"),
        .names = "FDR_{gsub('^p_', '', .col)}"
      )
    ) %>%
    ungroup()
}

markers_gene <- add_fdr_columns(markers_gene)
markers_pathways <- add_fdr_columns(markers_pathways)
                   
# Save
write_csv(markers_gene, file.path(intermediate_results_folder,paste0("DE_result_", id, "_genes.csv")))
write_csv(markers_pathways,
          file.path(intermediate_results_folder, paste0("DE_result_", id, "_pathways.csv")))

