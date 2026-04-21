library(scran)
#' Extract pathway colData column names
get_pathway_colnames <- function(sce) {
  progeny_cols <- grep("^PROGENy_", names(colData(sce)), value = TRUE)
  hallmark_cols <- grep("HALLMARK_", names(colData(sce)), value = TRUE)
  kegg_cols <- grep("KEGG_", names(colData(sce)), value = TRUE)
  c(progeny_cols, hallmark_cols, kegg_cols)
}

#' Filter SCE by peg cell count and classify pegs
#'
#' @param id Sample name (e.g., "A427PE_RMC")
#' @param intermediate_results_folder Path to folder with SCE files
#' @param min_cells Minimum number of cells per peg (default: 20)
#' @return A list with filtered SCE, control_pegs, and pegs_tested
filter_and_classify_pegs <- function(id, intermediate_results_folder, min_cells = 20) {
  sce <- readRDS(file.path(intermediate_results_folder, paste0("sce_", id, ".rds")))
  sce <- sce[, !(grepl(";", sce$peg))]
  sce$lib_class[grepl("Non", sce$peg)] <- "Non-Targeting"
  
  peg_counts <- as.data.frame(colData(sce)) %>%
    dplyr::group_by(peg, lib_class) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_cells >= min_cells)
  
  control_pegs <- peg_counts %>%
    dplyr::filter(lib_class == "synonymous") %>%
    dplyr::pull(peg)
  
  pegs_tested <- peg_counts %>%
    dplyr::filter(lib_class != "synonymous") %>%
    dplyr::pull(peg)
  
  sce <- sce[, sce$peg %in% c(pegs_tested, control_pegs)]
  
  list(
    sce = sce,
    control_pegs = control_pegs,
    pegs_tested = pegs_tested
  )
}

#' Run pathway-level lmer for a given sample id
#'
#' @param id Sample name (e.g., "A427PE_RMC")
#' @param intermediate_results_folder Folder where to save DE table and to locate the SCE
#' @param min_cells Minimum number of cells per peg (default: 20)
#' @return NULL; saves DE results
run_pathway_DE_for_id <- function(id, intermediate_results_folder, min_cells = 20) {
  sce_info <- filter_and_classify_pegs(id, intermediate_results_folder, min_cells)
  sce <- sce_info$sce
  pathway_cols <- get_pathway_colnames(sce)
  
  markers_pathways <- run_pathway_lmer(
    sce = sce,
    groups_tested = sce_info$pegs_tested,
    control_groups = sce_info$control_pegs,
    pathway_cols = pathway_cols,
    group_col = "peg",
    random_col = "puro",
    type_label = paste0("pathway_", id, "_condition")
  )
  
  saveRDS(
    markers_pathways,
    file.path(
      intermediate_results_folder,
      paste0("DE_result_", id, "_condition_pathways_uncorrected_temp.rds")
    )
  )
  NULL
}

