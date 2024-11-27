gc()
rm(list = ls())

new.pkgs <- c("APackOfTheClones", "svglite", 'ggthemes', "ggpubr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)
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
path.to.08.output <- file.path(path.to.main.output, "08_output")

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
  print(sprintf("reading in dataset: %s", dataset.name))
  if (dataset.name == "Dataset1"){
    all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("%s_merge_clusters.rds", dataset.name)))
  } else {
    all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))    
  }
  
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

dataset.name <- "Dataset1"
print("----------------------------------------------------------")
print(sprintf("working on dataset %s", dataset.name))
print("----------------------------------------------------------")
dir.create(file.path(path.to.07.output, dataset.name, "clone_APOTC"), showWarnings = FALSE, recursive = TRUE)

if (dataset.name == "Dataset1"){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

s.obj <- all.s.obj[[dataset.name]]
Idents(s.obj) <- "seurat_clusters"

s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = reduction.name, clonecall = "CTaa")
meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data$CTaa %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))
for (sampleid in unique(s.obj$name)){
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data$name == sampleid & meta.data$CTaa == x))
  }))
}

plot.sampleid <- "CD45exp1_m1"
tmp.clonedf <- clonedf[, c("clone", plot.sampleid)]
colnames(tmp.clonedf) <- c("clone", "count")
tmp.clonedf <- subset(tmp.clonedf, tmp.clonedf$count != 0)

top20.clones <- head(tmp.clonedf, 20) %>% pull(clone)
split.top5.clones <- split(top20.clones, seq(1,4))

colors <- tableau_color_pal(palette = "Tableau 20")(20)
split.colors <- split(colors, seq(1,4))

apotc.clone.plot <- list()
for (i in seq(1,4)){
  apotc.clone.plot[[i]] <- vizAPOTC(s.obj, clonecall = "CTaa", 
                                    verbose = FALSE, 
                                    reduction_base = reduction.name, 
                                    show_labels = TRUE, 
                                    legend_position = "top_right", 
                                    legend_sizes = 2) %>%
    showCloneHighlight(as.character(split.top5.clones[[i]]), fill_legend = TRUE, 
                       color_each = split.colors[[i]]) 
}

ggarrange()

# to-do 
# merge cluster 2 and 7 in dataset 1 in the origin data --> new UMAP, new APOTC. 
# APOTC same color with UMAP. 
# top 20 biggest clones in each sample --> 4 plots. 
# top 20 biggest clones for each cluster --> 4 plots, 5 clones per plot. 

