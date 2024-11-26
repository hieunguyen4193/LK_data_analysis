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
path.to.07.output <- file.path(path.to.main.output, "07_output")

dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

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

# dataset.name <- "Dataset2"
for (dataset.name in names(all.s.obj)){
  print("----------------------------------------------------------")
  print(sprintf("working on dataset %s", dataset.name))
  print("----------------------------------------------------------")
  dir.create(file.path(path.to.07.output, dataset.name, "clone_APOTC"), showWarnings = FALSE, recursive = TRUE)
  
  if (dataset.name == "Dataset1"){
    reduction.name <- "RNA_UMAP"
  } else {
    reduction.name <- "INTE_UMAP"
  }
  dir.create(file.path(path.to.06.output, dataset.name), showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- all.s.obj[[dataset.name]]
  
  ##### install APackOfTheClone package
  install.packages("APackOfTheClones")
  library(APackOfTheClones)
  s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = "INTE_UMAP", clonecall = "CTaa")
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("cell.barcode")
  clonedf <- data.frame(meta.data$CTaa %>% table())
  colnames(clonedf) <- c("clone", "count")
  clonedf <- clonedf %>% arrange(desc(count))
  
  for (cloneid in unique(subset(clonedf, clonedf$count >= 10)$clone)){
    print(sprintf("working on %s", cloneid))
    apotc.clone.plot <- vizAPOTC(s.obj, clonecall = "CTaa", verbose = FALSE, reduction_base = "INTE_UMAP", repulse = TRUE, show_labels = TRUE) %>%
      showCloneHighlight(cloneid, fill_legend = TRUE, )  
    ggsave(plot = apotc.clone.plot,
           filename = sprintf("%s.svg", cloneid),
           path = file.path(path.to.07.output, dataset.name, "clone_APOTC"),
           device = "svg",
           width = 15, 
           height = 10,
           dpi = 300)
  }
}

# umap.plot <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)
