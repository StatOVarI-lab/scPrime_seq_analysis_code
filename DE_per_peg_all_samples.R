
args <- commandArgs(trailingOnly = TRUE)

intermediate_results_folder <- args[1]
meta_data_file <- args[2]
id <- args[3]

library(scPrimeR)
library(tidyverse)
library(scran)

source("helper_functions.R")

peg_index <- as.integer(Sys.getenv("LSB_JOBINDEX"))

dir.create(file.path(intermediate_results_folder, "temp_DE", id), recursive = TRUE)



  sce_info <- filter_and_classify_pegs(id, intermediate_results_folder, min_cells=20)
  sce <- sce_info$sce
  if (peg_index > length(sce_info$pegs_tested)) {
  return(NULL)
}
  a <- run_gene_DE_nebula_for_peg(
  peg_j = sce_info$pegs_tested[peg_index],               
  sce = sce,                           
  control_pegs = sce_info$control_pegs,        
  id = id,                              
  intermediate_results_folder = 
    file.path(intermediate_results_folder, "temp_DE",id),
  cell_id_col = "puro"                 
)
    
