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

all.s.obj <- list()
for (dataset.name in names(dataset.names)){
  all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))
  meta.data <- all.s.obj[[dataset.name]]@meta.data %>% rownames_to_column("barcode")
  clonedf <- data.frame(table(meta.data$CTaa))%>% arrange(desc(Freq))
  colnames(clonedf) <- c("clone", "count")
  if (dataset.name == "Dataset1"){
    reduction.name <- "RNA_UMAP"
  } else {
    reduction.name <- "INTE_UMAP"
  }
  umapdf <- all.s.obj[[dataset.name]]@reductions[[reduction.name]]@cell.embeddings %>% as.data.frame() %>% rownames_to_column("barcode")
  meta.data <- merge(meta.data, umapdf, by.x = "barcode", by.y = "barcode")
  if (dataset.name == "Dataset1"){
    meta.data.GFP <- subset(meta.data, meta.data$seurat_clusters %in% c(3, 5, 6, 8))
    meta.data.CD45 <- subset(meta.data, meta.data$seurat_clusters %in% c(3, 5, 6, 8) == FALSE)
    
    clonedf.GFP  <- data.frame(table(meta.data.GFP$CTaa))%>% arrange(desc(Freq))
    colnames(clonedf.GFP) <- c("clone", "count")
    
    clonedf.CD45  <- data.frame(table(meta.data.CD45$CTaa))%>% arrange(desc(Freq))
    colnames(clonedf.CD45) <- c("clone", "count")
    
    writexl::write_xlsx(meta.data.GFP, file.path(path.to.02.output, sprintf("%s_GFP.metadata_with_cloneInfo.xlsx", dataset.name)))
    writexl::write_xlsx(clonedf.GFP, file.path(path.to.02.output, sprintf("%s_GFP.clonedf.xlsx", dataset.name)))  

    writexl::write_xlsx(meta.data.CD45, file.path(path.to.02.output, sprintf("%s_CD45.metadata_with_cloneInfo.xlsx", dataset.name)))
    writexl::write_xlsx(clonedf.CD45, file.path(path.to.02.output, sprintf("%s_CD45.clonedf.xlsx", dataset.name)))  
    
  } else {
    writexl::write_xlsx(meta.data, file.path(path.to.02.output, sprintf("%s.metadata_with_cloneInfo.xlsx", dataset.name)))
    writexl::write_xlsx(clonedf, file.path(path.to.02.output, sprintf("%s.clonedf.xlsx", dataset.name)))  
  }
}
