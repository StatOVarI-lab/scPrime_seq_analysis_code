library(tidyverse)
library(ggrepel)
library(pheatmap)
library(scPrimeR)

args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])
intermediate_results_folder <- args[2]
meta_data_file <- args[3]
id <- args[4]

fig_dir <- intermediate_results_folder

# Load results
markers_genes <- read_csv(file.path(intermediate_results_folder,
                        paste0("DE_result_", id, "_genes.csv")))
markers_pathway_df <- read_csv(file.path(intermediate_results_folder,
                        paste0("DE_result_", id, "_pathways.csv")))


# Annotate
markers_genes <- markers_genes %>%
  mutate(
    mutation_peg = if_else(
      !is.na(target_gene),
      paste0(target_gene, "_", mutation, ";", peg),
      paste0(target_gene, ";", peg)
    )
  ) %>%
  { if ("logFC_editedTRUE" %in% names(.)) mutate(., logFC = logFC_editedTRUE) else . } %>%
  { if ("FDR_editedTRUE" %in% names(.)) mutate(., FDR = FDR_editedTRUE) else . } %>%
  select(-any_of(c("FDR_editedTRUE", "logFC_editedTRUE")))
   
markers_pathway_df <- markers_pathway_df %>%
  mutate(
    mutation_peg = if_else(
      !is.na(target_gene),
      paste0(target_gene, "_", mutation, ";", peg),
      paste0(target_gene, ";", peg)
    )
  ) %>%
  { if ("logFC_editedTRUE" %in% names(.)) mutate(., logFC = logFC_editedTRUE) else . } %>%
  { if ("FDR_editedTRUE" %in% names(.)) mutate(., FDR = FDR_editedTRUE) else . } %>%
  select(-any_of(c("FDR_editedTRUE", "logFC_editedTRUE")))
fig_dir <- intermediate_results_folder                              


# Volcano plots
dir.create(file.path(fig_dir, paste0("volcano_plots_", id)), showWarnings = FALSE, recursive = TRUE)
                     

pegs_selected <- markers_genes %>%
  filter(FDR <= 0.1) %>%
  count(mutation_peg, name = "n_sig") %>%
  filter(n_sig >= 10) %>%
  pull(mutation_peg)

for (peg_j in pegs_selected) {
  mt <- markers_genes %>% filter(mutation_peg == peg_j)
  make_volcano(
    mt, "logFC", "FDR",
    paste0(peg_j, " ", id),
    file.path(fig_dir, paste0("volcano_plots_", id),
              paste0("volcano_edited_", id, "_", peg_j, ".pdf"))
  )
}


write_csv(markers_genes,
          file.path(intermediate_results_folder, 
            paste0("DE_results_", id, "_genes_improved_colnames.csv")))


# ##############################################################################
# # Volcano plots pathway
# ##############################################################################

pegs_selected_pathways <- markers_pathway_df %>%
  filter(FDR <= 0.1) %>%
  count(mutation_peg, name = "n_sig") %>%
  filter(n_sig >= 1) %>%
  pull(mutation_peg)

for (peg_j in pegs_selected_pathways) {
  mt <- markers_pathway_df %>% filter(mutation_peg == peg_j)
  make_volcano(
    mt, "logFC", "FDR",
    paste0(peg_j, " ", id),
    file.path(fig_dir, paste0("volcano_plots_", id),
              paste0("volcano_edited_pathways", id, "_", peg_j, ".pdf"))
  )
}


write_csv(markers_pathway_df,
          file.path(intermediate_results_folder, 
            paste0("DE_results_", id, "_pathways_improved_colnames.csv")))

