# Outlier-based QC for mRNA data, removing cells that are low outliers in 
# terms of the number of detected features or total umi count, or
# high outliers in terms of mitochondrial content. 

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
base_folder_cellranger <- args[2]
intermediate_results_folder <- args[3]

folders_cellranger_barcodes <- 
paste0(base_folder_cellranger,"/puro_v2/",
      c("A427PE_Control","A427PE_RMC","PC9PE_Control"))

folders_cellranger_pegRNAs <- 
paste0(base_folder_cellranger,"/pegRNAs_v2/",
       c("A427PE_Control","A427PE_RMC","PC9PE_Control"))

file_name_mRNA_matrices <- 
paste0(intermediate_results_folder,"/mRNA_matrices.rds")
file_name_mRNA_matrices_QC <- 
paste0(intermediate_results_folder,"/mRNA_matrices_QC.rds")
file_name_sces <- 
paste0(intermediate_results_folder,"/sce_list.rds")

ids <- c("A427PE_Control","A427PE_RMC","PC9PE_Control")

library(scPrimeR)

mRNA_matrices <- list()
for (k in 1:3){
 matrix_temp_barcodes <- read_in_filtered_matrix(folders_cellranger_barcodes[k])
 mRNA_matrices[[k]] <- read_in_filtered_matrix(folders_cellranger_pegRNAs[k])
 colnames(mRNA_matrices[[k]]) <- paste0(colnames(mRNA_matrices[[k]]), ids[k])
}
names(mRNA_matrices) <- ids
saveRDS(mRNA_matrices,file_name_mRNA_matrices)

mRNA_matrices <- readRDS(file_name_mRNA_matrices)

QC_list <- list()
for (j in 1:length(mRNA_matrices)){
  QC_list[[j]] <- QC_mRNA_outlier(mRNA_matrices[[j]],file_name = paste0("QC_",ids[j]))
}
names(QC_list) <- ids
saveRDS(QC_list,file=file_name_mRNA_matrices_QC)

set.seed(42)

QC_list <- readRDS(file_name_mRNA_matrices_QC)

sce_list <- list()
for (j in 1:length(QC_list)){

  sce_list[[j]] = computeSumFactors(QC_list[[j]]$sce)
  sce_list[[j]] <- logNormCounts(sce_list[[j]])
}
names(sce_list) <- ids

saveRDS(sce_list,file=file_name_sces)

cell_counts <- data.frame(
  sample = ids,
  before_QC = sapply(mRNA_matrices, ncol),
  after_QC = sapply(sce_list, ncol)
)


write.csv(cell_counts, file = 
          paste0(intermediate_results_folder,"/cell_counts_summary.csv"), row.names = FALSE)

qc_cutoffs <- data.frame(
  sample = ids,
  lib_size_lower = sapply(QC_list, function(x) x$cutoffs$lib_size["lower"]),
  n_features_lower = sapply(QC_list, function(x) x$cutoffs$n_features["lower"]),
  mito_percent_upper = sapply(QC_list, function(x) x$cutoffs$mito_percent["higher"])
)

write.csv(qc_cutoffs, file = paste0(intermediate_results_folder,
                                    "/qc_cutoffs_summary.csv"), row.names = FALSE)