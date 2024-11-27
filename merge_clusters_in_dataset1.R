gc()
rm(list = ls())

new.pkgs <- c("APackOfTheClones", "svglite")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

path.to.main.src <- "/media/hieunguyen/HNSD01/src/LK_data_analysis"
source(file.path(path.to.main.src, "helper_functions.R"))

# outdir <- "/media/hieunguyen/HNSD_mini/outdir/LK_data_analysis"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/LK_data_analysis"

path.to.main.output <- file.path(outdir, "general_outputs")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output") 
path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.05.output <- file.path(path.to.main.output, "05_output")
path.to.06.output <- file.path(path.to.main.output, "06_output")
path.to.07.output <- file.path(path.to.main.output, "07_output")

dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

dataset.names <- list(
  Dataset1 = "1stExp_Kopplin",
  Dataset2 = "211227_Kopplin",
  Dataset3 = "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9",
  Dataset4 = "230316_Kopplin"
)

dataset.name <- "Dataset1"
s.obj <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode2") %>%
  rowwise() %>%
  mutate(seurat_clusters.org = seurat_clusters) %>%
  mutate(seurat_clusters = ifelse(seurat_clusters.org %in% c(2, 7), 2, seurat_clusters)) %>%
  column_to_rownames("barcode2") 

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, col.name = "seurat_clusters", metadata = meta.data$seurat_clusters)
saveRDS(object = s.obj, file.path(path.to.01.output, sprintf("%s_merge_clusters.rds", dataset.name)))

# DimPlot(object = s.obj, reduction = "RNA_UMAP", group.by = "seurat_clusters", label = TRUE, label.box =TRUE)
