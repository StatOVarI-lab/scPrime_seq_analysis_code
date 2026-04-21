args <- commandArgs(trailingOnly = TRUE)
# argument list is trailing list of path to all SCEs
args <-  c("/lustre/scratch125/casm/staging/team363/ms58/scPrime/intermediate_results_scPrime/sce_A427PE_RMC.rds",
           "/lustre/scratch125/casm/staging/team363/ms58/scPrime/intermediate_results_scPrime/sce_A427PE_Control.rds",
           "/lustre/scratch125/casm/staging/team363/ms58/scPrime/intermediate_results_scPrime/sce_PC9PE_Control.rds")


sce_input_files <- args

library(MAYA)
library(tidyverse)
library(SingleCellExperiment)


sce_list <- lapply(sce_input_files,readRDS)

add_MAYA_scores <- function(sce) {
  logcounts_sce <- logcounts(sce)
  
  activity_summary_hallmark <- MAYA_pathway_analysis(
    expr_mat = logcounts_sce,
    modules_list = "hallmark",
    is_logcpm = TRUE,
    min_cells_pct = 0.05,
    min_genes = 3,
    max_contrib = 0.9,
    compute_umap = FALSE
  )
  activity_summary_kegg <- MAYA_pathway_analysis(
    expr_mat = logcounts_sce,
    modules_list = "kegg",
    is_logcpm = TRUE,
    min_cells_pct = 0.05,
    min_genes = 3,
    max_contrib = 0.9,
    compute_umap = FALSE
  )
  
  mat_hallmark <- activity_summary_hallmark$activity_matrix
  mat_kegg <- activity_summary_kegg$activity_matrix
    
  # Add to colData
  colData(sce) <- cbind(colData(sce), t(mat_hallmark), t(mat_kegg))
  sce
}

remove_MAYA_scores <- function(sce) {# remove any MAYA scores in case they are already present
  cd <- colData(sce)
  # Remove columns that start with "HALLMARK_" or "KEGG_"
  keep_cols <- !grepl("^(HALLMARK_|KEGG_)", colnames(cd))
  colData(sce) <- cd[, keep_cols, drop = FALSE]
  sce
}

sce_list <- lapply(sce_list, remove_MAYA_scores)
sce_list <- lapply(sce_list, add_MAYA_scores)

# Save SCEs
for (i in seq_along(sce_list)) {
  saveRDS(sce_list[[i]], sce_input_files[i])
}