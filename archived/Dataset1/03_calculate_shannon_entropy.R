gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
outdir <- "/media/hieunguyen/HNSD_mini/outdir/LK_data_analysis"

path.to.main.output <- file.path(outdir, "general_outputs")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output") 
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

dataset.names <- list(
  Dataset1 = "1stExp_Kopplin",
  Dataset2 = "211227_Kopplin",
  Dataset3 = "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9",
  Dataset4 = "230316_Kopplin"
)

all.clonedf <- hash()
all.metadata <- hash()
for (dataset.name in names(dataset.names)){
  all.metadata[[dataset.name]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s.metadata_with_cloneInfo.xlsx", dataset.name)))
  all.clonedf[[dataset.name]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s.clonedf.xlsx", dataset.name)))
}

#####---------------------------------------------------------------------------------------------#####
##### helper functions: Calculate Shannon entropy to measure "diversity" of cells/clones over clusters
#####---------------------------------------------------------------------------------------------#####
calculate_shannon_entropy <- function(cell.list, input.metadata){
  ##### by default, we assume that the column "seurat_clusters" contains the active cluster numbers.
  N <- length(unique(input.metadata$seurat_clusters))
  count_clone_in_cluster <- table(subset(input.metadata, input.metadata$barcode %in% cell.list) %>% 
                                    subset(select = c(seurat_clusters))) %>%  as.data.frame()
  count_clone_in_cluster <- count_clone_in_cluster %>% rowwise %>% mutate(Prob = Freq / sum(count_clone_in_cluster$Freq)) %>%
    subset(Prob != 0)
  shannon_entropy <- -sum(count_clone_in_cluster$Prob * log2(count_clone_in_cluster$Prob))/log2(N)
  return(shannon_entropy)
}

# To calculate Shannon entropy for a clone, we extract its list of barcodes (cell.list) and feed the function calculate_shannon_entropy_sampling_cells
calculate_shannon_entropy_clone <- function(input.clone, input.metadata){
  cell.list <- subset(input.metadata, input.metadata$CTaa == input.clone)$barcode
  shannon_entropy <- calculate_shannon_entropy(cell.list, input.metadata)
  return(shannon_entropy)
}

calculate_shannon_entropy_clone(input.clone = "CAMNTGYQNFYF_CTCSVQDTQYF",
                                input.metadata = all.metadata$Dataset1)

#> example: clone = "CAMNTGYQNFYF_CTCSVQDTQYF"
