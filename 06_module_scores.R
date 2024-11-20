gc()
rm(list = ls())

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

dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

dataset.names <- list(
  Dataset1 = "1stExp_Kopplin",
  Dataset2 = "211227_Kopplin",
  Dataset3 = "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9",
  Dataset4 = "230316_Kopplin"
)

all.clonedf <- hash()
all.metadata <- hash()
for (dataset.name in names(dataset.names)){
  if (dataset.name == "Dataset1"){
    all.metadata[["Dataset1_CD45"]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s_CD45.metadata_with_cloneInfo.xlsx", dataset.name)))
    all.clonedf[["Dataset1_CD45"]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s_CD45.clonedf.xlsx", dataset.name)))
    
    all.metadata[["Dataset1_GFP"]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s_GFP.metadata_with_cloneInfo.xlsx", dataset.name)))
    all.clonedf[["Dataset1_GFP"]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s_GFP.clonedf.xlsx", dataset.name)))
  } else {
    all.metadata[[dataset.name]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s.metadata_with_cloneInfo.xlsx", dataset.name)))
    all.clonedf[[dataset.name]] <- readxl::read_excel(file.path(path.to.02.output, sprintf("%s.clonedf.xlsx", dataset.name)))    
  }
}

all.entropydf <- hash()
for (dataset.name in names(all.clonedf)){
  all.entropydf[[dataset.name]] <- readxl::read_excel(file.path(path.to.03.output, sprintf("%s.Shannon_entropy.xlsx", dataset.name)))
}

all.s.obj <- list()
for (dataset.name in names(dataset.names)){
  all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))
}

module.gene.list <- list(
  Th1 = c("Il12rb1", "Il12rb2", "Tbx21", "Ifng", "Cxcr3", "Ccl4"),
  Th17 = c("Il17a", "Il17f", "Rorc", "Il23r", "Il1r1", "Ccr6", "Il17re", "Il22", "Ccl20"),
  Th2 = c("Il4", "Il5", "Gata3", "Il17rb"),
  Cytotoxic = c("Gzma", "Gzmb", "Gzmk", "Gzmh", "Gnly"),
  Tr1 = c("Il10", "Gzma", "Irf4", "Batf", "Maf", "Ifng"),
  Treg = c("Foxp3", "Ctla4", "Il10", "Ikzf2", "Tnfrsf18"),
  Tfh = c("Cd200", "Tox2", "Tox", "Cxcr5", "Pdcd1"),
  Tcm = c("Il7r", "Klf2", "Tob1", "S1pr1", "S1pr4", "Ccr7", "Sell", "Tcf7", "Left1"),
  Proliferating = c("Mki67", "Stmn1", "Top2a", "Birc5")
)


for (dataset.name in names(all.s.obj)){
  dir.create(file.path(path.to.06.output, dataset.name), showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- all.s.obj[[dataset.name]]
  for (input.list in names(module.gene.list)){
    DefaultAssay(s.obj) <- "RNA"
    s.obj <- AddModuleScore(object = s.obj, features = list(module.gene.list[[input.list]]), name = sprintf("%s_", input.list), ctrl = 50)
  }
  
  meta.data <- s.obj@meta.data
  fake.module.gene.list <- to_vec(
    for (item in names(module.gene.list)){
      sprintf("%s_1", item)
    }
  )
  library(viridis)
  if ("svglite" %in% installed.packages() == FALSE){
    install.packages("svglite")
  }
  
  feature.plot <- FeaturePlot(object = s.obj, reduction = "INTE_UMAP", ncol = 3, features = fake.module.gene.list, pt.size = 1) &
    scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
  
  violin.plot <- VlnPlot(object = s.obj, features = fake.module.gene.list, pt.size = 0)
  
  heatmapdf <- s.obj@meta.data[, c("seurat_clusters", fake.module.gene.list)] %>%
    group_by(seurat_clusters) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("seurat_clusters")
  
  colors <- rev(RColorBrewer::brewer.pal(n = length(fake.module.gene.list), name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  
  heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
  colnames(heatmapdf.scaled) <- to_vec(
    for (item in colnames(heatmapdf.scaled)){
      str_split(item, "_")[[1]][[1]]
    }
  )
  heatmap.plot <- heatmapdf.scaled %>% rownames_to_column("cluster") %>% 
    pivot_longer(!cluster, names_to = "signature", values_to = "z_score") %>%
    ggplot(aes(x = cluster, y = signature, fill = z_score)) + geom_tile() + 
    scale_fill_distiller(palette = "RdBu") + 
    theme(axis.text = element_text(size = 22))
  
  ggsave(plot = feature.plot, filename = sprintf("feature_plot_module_scores.svg"), 
         path = file.path(path.to.06.output, dataset.name), device = "svg", width = 20, height = 14)
  ggsave(plot = violin.plot, filename = sprintf("violin_plot_module_scores.svg"), 
         path = file.path(path.to.06.output, dataset.name), device = "svg", width = 20, height = 14)
  ggsave(plot = heatmap.plot, filename = sprintf("heatmap_module_scores.svg"), 
         path = file.path(path.to.06.output, dataset.name), device = "svg", width = 10, height = 10)
}

