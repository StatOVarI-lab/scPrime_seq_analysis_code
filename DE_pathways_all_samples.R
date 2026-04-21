# Differential expression results at pathway level


args <- commandArgs(trailingOnly = TRUE)

intermediate_results_folder <- args[1]
meta_data_file <- args[2]
id <- args[3]

library(scPrimeR)
library(tidyverse)
library(scran)

source("helper_functions.R")

run_pathway_DE_for_id(id, intermediate_results_folder, min_cells = 20)   

