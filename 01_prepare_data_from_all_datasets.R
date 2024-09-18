gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
if ("Trex" %in% installed.packages()){
  library(Trex)  
}

outdir <- "/media/hieunguyen/HNSD_mini/outdir/LK_data_analysis"

path.to.main.output <- file.path(outdir, "general_outputs")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

dataset.names <- list(
  Dataset1 = "1stExp_Kopplin",
  Dataset2 = "211227_Kopplin",
  Dataset3 = "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9",
  Dataset4 = "230316_Kopplin"
)

path.to.s.obj <- list(
  Dataset1 = file.path(outdir, 
                       "1stExp_Kopplin", 
                       "data_analysis", 
                       "01_output", 
                       "merged_all_first_exp_dataset.rds"),
  Dataset2 = file.path(outdir, 
                       "211227_Kopplin", 
                       "data_analysis",
                       "01_output",
                       "merged_all_second_exp_dataset.annotated.filteredCD45.integrated.rds"),
  Dataset3 = file.path(outdir, 
                       "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9",
                       "1st_round", 
                       "pct_mito_10_1", 
                       "data_analysis", 
                       "01_output",
                       "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9.addVDJ.rds"),
  Dataset4 = file.path(outdir, 
                       "230316_Kopplin",
                       "1st_round",
                       "pct_mito_10_1",
                       "data_analysis",
                       "01_output",
                       "230316_Kopplin.seurat.obj.removed.14_15_16.addedVDJ.integrated.rds")
)

all.s.obj <- list()
for (dataset.name in names(path.to.s.obj)){
  print(sprintf("adding data from the object %s", dataset.name))
  all.s.obj[[dataset.name]] <- readRDS(path.to.s.obj[[dataset.name]])
}

sample.list <- list()
for (dataset.name in names(all.s.obj)){
  sample.list[[dataset.name]] <- unique(all.s.obj[[dataset.name]]$name)
}
  
path.to.storage <- "/media/hieunguyen/HNSD01/storage/LKopplin_data"

all.vdj.files <- Sys.glob(file.path(path.to.storage, "*/*/*", "filtered_contig_annotations.csv"))
sample.names <- unlist(lapply(lapply(all.vdj.files, dirname), basename))
names(all.vdj.files) <- sample.names

##### add VDJ information and modify seurat objects for all datasets
for (dataset.name in names(all.s.obj)){
  cluster.resolution <- 1
  if (dataset.name == "Dataset1"){
    ##### increase the cluster resolution of the first dataset. 
    print(sprintf("CLUSTER RESOLUTIONS: %s", cluster.resolution))
    chosen.assay <- "RNA"
    chosen.seed <- 42
    num.dim.integration <- 25 
    num.PCA <- 25
    num.dim.cluster <- 25
    num.PC.used.in.Clustering <- 25
    
    s.obj <- all.s.obj[[dataset.name]]
    s.obj <- FindNeighbors(s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
    s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = chosen.seed)
    DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, pt.size = 0.5, label.size = 8) + 
      xlim(-8, 7) + ylim(-4, 5)
    saveRDS(object = s.obj, file = file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))  
    all.s.obj[[dataset.name]] <- s.obj
  }
  if (file.exists(file.path(path.to.01.output, sprintf("%s.rds", dataset.name))) == FALSE){
    if (dataset.name == "Dataset1"){
      reduction.name <- "RNA_UMAP"
    } else {
      reduction.name <- "INTE_UMAP"
    }
    all.s.obj[[dataset.name]]$CTaa <- NULL
    all.contig.files <- all.vdj.files[sample.list[[dataset.name]]]
    
    contig_list <- lapply(all.contig.files, vroom, show_col_type = FALSE)
    names(contig_list) <- sample.list[[dataset.name]]
    combined.contigs <- combineTCR(contig_list,
                                   samples = sample.list[[dataset.name]],
                                   ID = sample.list[[dataset.name]],
                                   removeNA=FALSE, 
                                   removeMulti=FALSE, 
                                   cells = "T-AB")
    names(combined.contigs) <- sample.list[[dataset.name]]
    
    all.s.obj[[dataset.name]] <- combineExpression(combined.contigs, all.s.obj[[dataset.name]], cloneCall="aa")
    saveRDS(all.s.obj[[dataset.name]], file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))
  }
}

