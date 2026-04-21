# Basic scRNA-seq analysis per sample, including cell cycle and pathway activity scores. 

args <- commandArgs(trailingOnly = TRUE)
#args <-  c("/lustre/scratch125/casm/staging/team363/ms58", "/lustre/scratch125/casm/staging/team363/ms58/scPrime/intermediate_results_scPrime", "/lustre/scratch125/casm/staging/team363/ms58/scPrime/input_files/scRNAseq_zscores.csv", "A427PE_Control") 


setwd(args[1])
intermediate_results_folder <- args[2]
meta_data_file <- args[3]
id <- args[4]


puro_peg_assignment_file <- paste0(intermediate_results_folder, "/puro_peg/merged_puro_peg_assignments_extended_",id,".csv")

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scPrime))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(ggthemes))
suppressPackageStartupMessages(library(progeny))
suppressPackageStartupMessages(library(S4Vectors))

folder_puro_peg_assignments <- paste0(intermediate_results_folder, "/puro_peg")
file_name_sces <- paste0(intermediate_results_folder, "/sce_list.rds")
file_name_sce_processed <- paste0(intermediate_results_folder, "/sce_", id, ".rds")
figure_folder <- paste0("~/scPrime/figures/", id)
dir.create(figure_folder, recursive = TRUE, showWarnings = FALSE)

sce_list <- readRDS(file_name_sces)
sce <- sce_list[[id]]

meta_data <- read_delim(meta_data_file, delim = ",")

if (grepl("RMC", id)) {
  meta_data <- meta_data %>% filter(drug == "RMC6236")
} else if (grepl("A427", id)) {
  meta_data <- meta_data %>% filter(drug == "control", cell_model == "A427")
} else if (grepl("PC9PE", id)) {
  meta_data <- meta_data %>% filter(drug == "control", cell_model == "PC9")
}

puro_peg_assignment <- read_delim(puro_peg_assignment_file, delim = ",")

coldata_df <- as.data.frame(colData(sce))
coldata_df <- coldata_df %>% 
left_join(puro_peg_assignment, by = c("cell_id" = "cell"))

coldata_df <- coldata_df %>% 
left_join(meta_data,  by = c("peg" = "id"))

colData(sce) <- DataFrame(coldata_df)
sce <- sce[,!(is.na(sce$peg))]

sce$single_peg <- !(grepl(";",sce$peg))

coldata_df <- as.data.frame(colData(sce))

df_class <- coldata_df %>%
  select(puro, peg, lib_class, single_peg) %>%
  filter(!(is.na(puro))) %>% 
  filter(single_peg)

nr_pegs <- df_class %>%
  filter(!is.na(peg)) %>%
  group_by(peg, lib_class) %>%
  summarise(nr_cells = n(), .groups = "drop")

p_cells_per_lib_class <- ggplot(nr_pegs, aes(x = nr_cells, color = lib_class)) +
  geom_density(linewidth = 0.5) +
  scale_color_colorblind() +
  theme_classic(base_size = 8) +
  xlab("number of cells per peg") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))

nr_cells_by_lib_class <- nr_pegs %>%
  group_by(lib_class) %>%
  summarise(total_cells = sum(nr_cells), .groups = "drop")

ggsave(p_cells_per_lib_class, file = file.path(figure_folder, "nr_pegs_density.pdf"),
       width = 6, height = 5, units = "cm")


sce <- cell_cycle_scoring(sce)

df_cc <- as_tibble(colData(sce)) %>%
  select(lib_class, phase) %>%
  filter(!is.na(lib_class), !is.na(phase))

cc_count <- df_cc %>%
  group_by(lib_class, phase) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(lib_class) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

p_cc <- ggplot(cc_count, aes(x = lib_class, y = proportion, fill = phase)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("G1" = "darkgrey", "G2M" = "black", "S" = "purple")) +
  theme_bw(base_size = 8) +
  ylab("proportion in phase") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(p_cc, file = file.path(figure_folder, "lib_class_phase_proportion.pdf"),
       width = 6, height = 5, units = "cm")

progeny_scores <- progeny(as.matrix(logcounts(sce)), scale = FALSE)
colnames(progeny_scores) <- paste0("PROGENy_", colnames(progeny_scores))

colData(sce) <- cbind(colData(sce), progeny_scores)

progeny_long <- as_tibble(colData(sce)) %>%
  select(puro, lib_class, starts_with("PROGENy_")) %>%
  pivot_longer(cols = starts_with("PROGENy_"),
               names_to = "pathway", values_to = "score")

progeny_summary <- progeny_long %>%
  group_by(puro, pathway, lib_class) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(norm_score = scale(mean_score)[, 1]) %>%
  ungroup()

p_progeny <- ggplot(progeny_summary, aes(x = norm_score, color = lib_class)) +
  geom_density(linewidth = 0.5) +
  facet_wrap(~ pathway, scales = "free") +
  scale_color_colorblind() +
  theme_classic(base_size = 8) +
  xlab("Normalised PROGENy score (z-score)") +
  ylab("Density") +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(p_progeny, file = file.path(figure_folder, "progeny.pdf"),
       width = 18, height = 15, units = "cm")

saveRDS(sce, file = file_name_sce_processed)