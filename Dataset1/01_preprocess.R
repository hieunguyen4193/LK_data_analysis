#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
gc()
rm(list = ls())

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline"
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
PROJECT <- "1stExp_Kopplin"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/LK_data_analysis"

chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
cluster.resolution <- 0.5
num.PC.used.in.Clustering <- 25

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "s8_integration_and_clustering.R"))

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")# keep!
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output") # output for this script only!
dir.create(path.to.01.output, showWarnings = FALSE)

analysis.round.1st.exp <- "2nd"
path.to.first.exp.outdir <- file.path(outdir, PROJECT)
path.to.first.exp <- file.path(path.to.first.exp.outdir, sprintf("%s_round", analysis.round.1st.exp))

#####----------------------------------------------------------------------#####
# collect all output from 1st dataset 
#####----------------------------------------------------------------------#####
all_first_exprs <- Sys.glob(file.path(path.to.first.exp, "*"))
names(all_first_exprs) <-  to_vec(for(exprs in all_first_exprs) str_replace(basename(exprs), sprintf("_%s_round", analysis.round.1st.exp), "")) 

#####----------------------------------------------------------------------#####
# PROCESSING THE FIRST DATASET
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "merged_all_first_exp_dataset.rds")) == FALSE){
  data.list <- list()
  
  for (i in seq_along(all_first_exprs)){
    data.list[[i]] <- readRDS(Sys.glob(file.path(all_first_exprs[[i]], "s8a_output", "*.rds")))
  }
  
  first.exp.s.obj <- merge(x = data.list[[1]],
                           y = unlist(data.list[2:length(data.list)]),
                           merge.data = FALSE,
                           add.cell.ids = names(all_first_exprs),
                           project = PROJECT)
  
  chosen.assay <- "RNA"
  DefaultAssay(first.exp.s.obj) <- chosen.assay
  
  first.exp.s.obj <- NormalizeData(first.exp.s.obj) # ---> use Log Normalized
  first.exp.s.obj <- FindVariableFeatures(first.exp.s.obj, selection.method = "vst")
  first.exp.s.obj <- ScaleData(first.exp.s.obj, features = rownames(first.exp.s.obj))
  
  first.exp.s.obj <- RunPCA(first.exp.s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  first.exp.s.obj <- RunUMAP(first.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                             reduction.name=sprintf("%s_UMAP", chosen.assay), seed.use = chosen.seed)
  first.exp.s.obj <- RunTSNE(first.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                             reduction.name=sprintf("%s_TSNE", chosen.assay), seed.use = chosen.seed)
  
  first.exp.s.obj <- FindNeighbors(first.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  first.exp.s.obj <- FindClusters(first.exp.s.obj, resolution = cluster.resolution, random.seed = chosen.seed)
  
  saveRDS(object = first.exp.s.obj, file = file.path(path.to.01.output, "merged_all_first_exp_dataset.rds")) 
} else {
  first.exp.s.obj <- readRDS(file.path(path.to.01.output, "merged_all_first_exp_dataset.rds"))
}


#####----------------------------------------------------------------------#####
# EOF
#####----------------------------------------------------------------------#####
