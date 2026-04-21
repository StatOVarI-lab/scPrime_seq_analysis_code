# Robustly assigns puro-barcodes and pegRNAs to cells using a consensus of 
# probabilistic models of of skew-normal distributions, and ensuring agreement
# of puro-barcodes and pegRNAs.

args <- commandArgs(trailingOnly = TRUE)

setwd(args[1])
base_folder_cellranger <- args[2]
intermediate_results_folder <- args[3]

library(scPrimeR)
library(tidyverse)


folders_cellranger_barcodes <- 
paste0(base_folder_cellranger,"/puro_v2/",
      c("A427PE_Control","A427PE_RMC","PC9PE_Control"))

folders_cellranger_pegRNAs <- 
paste0(base_folder_cellranger,"/pegRNAs_v2/",
       c("A427PE_Control","A427PE_RMC","PC9PE_Control"))

folder_puro_peg_assignments <- paste0(intermediate_results_folder,"/puro_peg")

file_name_puro_matrices <- 
paste0(folder_puro_peg_assignments,"/puro_matrices.rds")

file_name_peg_matrices <- 
paste0(folder_puro_peg_assignments,"/peg_matrices.rds")

file_name_sces <- 
paste0(intermediate_results_folder,"/sce_list.rds")

figure_folder <- paste0(folder_puro_peg_assignments,"/figures_peg_puro_assignment")

ids <- c("A427PE_Control","A427PE_RMC","PC9PE_Control")



set.seed(12345)

peg_matrices <- lapply(folders_cellranger_pegRNAs, read_in_filtered_matrix, keep_mRNA=FALSE)
puro_matrices <- lapply(folders_cellranger_barcodes, read_in_filtered_matrix, keep_mRNA=FALSE)

for (k in 1:3){
  colnames(peg_matrices[[k]]) <- paste0(colnames(peg_matrices[[k]]), ids[k])
  colnames(puro_matrices[[k]]) <- paste0(colnames(puro_matrices[[k]]), ids[k])
}
names(puro_matrices) <- ids
names(peg_matrices) <- ids


param_list_peg <- list(
  list(g = 3, thresh = 2),
  list(g = 3, thresh = 3)
)

param_list_puro = param_list_peg


run_full_barcode_assignment_pipeline(
  puro_matrices = puro_matrices,
  peg_matrices = peg_matrices,
  ids = ids,
  param_list_puro,
  param_list_peg,
  output_dir = folder_puro_peg_assignments,
  figure_dir = figure_folder
) 