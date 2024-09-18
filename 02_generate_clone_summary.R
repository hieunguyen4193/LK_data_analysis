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
  if (dataset.name == "1stExp_Kopplin") {
    if (file.exists(file.path(path.to.01.output, sprintf("1stExp_Kopplin.cluster_res_%s.rds", cluster.resolution))) == FALSE){
      ##### increase the cluster resolution of the first dataset. 
      cluster.resolution <- 1
      
      print(sprintf("CLUSTER RESOLUTIONS: %s", cluster.resolution))
      chosen.assay <- "RNA"
      s.obj <- FindNeighbors(s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
      s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = chosen.seed)
      DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, pt.size = 0.5, label.size = 8) + 
        xlim(-8, 7) + ylim(-4, 5)
      saveRDS(object = s.obj, file = file.path(path.to.01.output, sprintf("1stExp_Kopplin.cluster_res_%s.rds", cluster.resolution)))  
    } else {
      all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("1stExp_Kopplin.cluster_res_%s.rds", cluster.resolution)))
    }
  }
  meta.data <- all.s.obj[[dataset.name]]@meta.data %>% subset(select = -c(barcode)) %>% rownames_to_column("barcode")
  clonedf <- data.frame(table(meta.data$CTaa))%>% arrange(desc(Freq))
  writexl::write_xlsx(meta.data, file.path(path.to.02.output, sprintf("%s.metadata_with_cloneInfo.xlsx", dataset.name)))
  writexl::write_xlsx(clonedf, file.path(path.to.02.output, sprintf("%s.clonedf.xlsx", dataset.name)))
}
