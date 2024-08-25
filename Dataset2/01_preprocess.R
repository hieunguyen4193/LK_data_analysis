#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
gc()
rm(list = ls())

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline"

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "helper_functions.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
PROJECT <- "211227_Kopplin"

chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
cluster.resolution <- 0.5
num.PC.used.in.Clustering <- 25
num.PC.used.in.UMAP <- 25

path.to.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(path.to.src, "s8_integration_and_clustering.R"))

outdir <- "/media/hieunguyen/HNSD_MBPro/CRC1382/outdir/LKopplin_OFFICIAL"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")# keep!
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.01.output <- file.path(path.to.main.output, "01_output") # output for this script only!
dir.create(path.to.01.output, showWarnings = FALSE)

analysis.round.1st.exp <- "2nd"
analysis.round.2nd.exp <- "2nd"

path.to.first.exp.outdir <- file.path(outdir ,"1stExp_Kopplin")
path.to.second.exp.outdir <- file.path(outdir, PROJECT)

path.to.first.exp <- file.path(path.to.first.exp.outdir, sprintf("%s_round", analysis.round.1st.exp))
path.to.second.exp <- file.path(path.to.second.exp.outdir, sprintf("%s_round", analysis.round.2nd.exp))

#####----------------------------------------------------------------------#####
# collect all output from 1st dataset 
#####----------------------------------------------------------------------#####
all_first_exprs <- Sys.glob(file.path(path.to.first.exp, "*"))
names(all_first_exprs) <-  to_vec(for(exprs in all_first_exprs) str_replace(basename(exprs), sprintf("_%s_round", analysis.round.1st.exp), "")) 

#####----------------------------------------------------------------------#####
# collect all output from 2nd dataset 
#####----------------------------------------------------------------------#####
all_second_exprs <- Sys.glob(file.path(path.to.second.exp, "*"))
names(all_second_exprs) <-  to_vec(for(exprs in all_second_exprs) str_replace(basename(exprs), sprintf("_%s_round", analysis.round.2nd.exp), "")) 

#####----------------------------------------------------------------------#####
# PROCESSING THE SECOND DATASET 
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "merged_all_second_exp_dataset.rds")) == FALSE){
  data.list <- list()
  for (i in seq_along(all_second_exprs)){
    data.list[[i]] <- readRDS(Sys.glob(file.path(all_second_exprs[[i]], "s8a_output", "*.rds")))
  }
  
  second.exp.s.obj <- merge(x = data.list[[1]], 
                            y = unlist(data.list[2:length(data.list)]),
                            merge.data = FALSE, 
                            add.cell.ids = names(all_second_exprs), 
                            project = PROJECT)
  
  chosen.assay <- "RNA"
  DefaultAssay(second.exp.s.obj) <- chosen.assay
  
  second.exp.s.obj <- NormalizeData(second.exp.s.obj) # ---> use Log Normalized
  second.exp.s.obj <- FindVariableFeatures(second.exp.s.obj, selection.method = "vst")
  second.exp.s.obj <- ScaleData(second.exp.s.obj, features = rownames(second.exp.s.obj))
  
  second.exp.s.obj <- RunPCA(second.exp.s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  second.exp.s.obj <- RunUMAP(second.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                              reduction.name=sprintf("%s_UMAP", chosen.assay), seed.use = chosen.seed)
  second.exp.s.obj <- RunTSNE(second.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                              reduction.name=sprintf("%s_TSNE", chosen.assay), seed.use = chosen.seed)
  saveRDS(object = second.exp.s.obj, file = file.path(path.to.01.output, "merged_all_second_exp_dataset.rds"))
} else {
  second.exp.s.obj <- readRDS(file.path(path.to.01.output, "merged_all_second_exp_dataset.rds"))
}

first.exp.s.obj <- readRDS(file.path(outdir, "1stExp_Kopplin", "data_analysis", "01_output", "merged_all_first_exp_dataset.rds"))

#####----------------------------------------------------------------------#####
# INTEGRATE THE FIRST DATASET AND TRANSFER THEIR LABELS TO THE SECOND
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "merged_all_second_exp_dataset.annotated.rds")) == FALSE){
  chosen.assay <- "RNA"
  integration.anchors <- FindTransferAnchors(reference = first.exp.s.obj, query = second.exp.s.obj,
                                             dims = 1:25, reference.reduction = "RNA_PCA")
  
  ##### Transfer label from anchors to query dataset
  predictions <- TransferData(anchorset = integration.anchors, refdata = first.exp.s.obj$stage,
                              dims = 1:25)
  
  second.exp.s.obj <- AddMetaData(second.exp.s.obj, metadata = predictions$predicted.id, col.name = "prediction")
  
  filtered.second.exp.s.obj <- subset(second.exp.s.obj, prediction != "CD45")
  filtered.second.exp.s.obj <- RunPCA(filtered.second.exp.s.obj, npcs = 25, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  filtered.second.exp.s.obj <- RunUMAP(filtered.second.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay),
                                       dims = 1:25, reduction.name=sprintf("%s_UMAP", chosen.assay))
  filtered.second.exp.s.obj <- RunTSNE(filtered.second.exp.s.obj, reduction = sprintf("%s_PCA", chosen.assay),
                                       dims = 1:25, reduction.name=sprintf("%s_TSNE", chosen.assay))
  filtered.second.exp.s.obj <- s8.integration.and.clustering(s.obj = filtered.second.exp.s.obj, 
                                                             path.to.output = path.to.01.output, 
                                                             save.RDS.s8 = FALSE,
                                                             PROJECT = PROJECT, 
                                                             num.dim.integration = num.dim.integration,
                                                             num.PCA = num.PCA,
                                                             num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                             num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                             cluster.resolution = 0.5,
                                                             my_random_seed = 42,
                                                             umap.method = "uwot",
                                                             genes.to.not.run.PCA = NULL,
                                                             inte_pca_reduction_name = "INTE_PCA", 
                                                             inte_umap_reduction_name = "INTE_UMAP")
  
  saveRDS(object = second.exp.s.obj, file = file.path(path.to.01.output, "merged_all_second_exp_dataset.annotated.rds"))
  saveRDS(object = filtered.second.exp.s.obj, file = file.path(path.to.01.output, "merged_all_second_exp_dataset.annotated.filteredCD45.integrated.rds"))
  
} else {
  second.exp.s.obj <- readRDS(file.path(path.to.01.output, "merged_all_second_exp_dataset.annotated.rds"))
  filtered.second.exp.s.obj <- readRDS(file.path(path.to.01.output, "merged_all_second_exp_dataset.annotated.filteredCD45.integrated.rds"))
}

#####----------------------------------------------------------------------#####
# DIMENSIONAL REDUCTION PLOT
#####----------------------------------------------------------------------#####

# all clusters, integrated
DimPlot(object = filtered.second.exp.s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)

# all clusters, before integration
DimPlot(object = filtered.second.exp.s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)
# Sample GFP_m3 separated!

# group by samples
DimPlot(object = filtered.second.exp.s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name") +
  ggtitle("All 3 samples on UMAP integrated")

DimPlot(object = filtered.second.exp.s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name") +
  ggtitle("All 3 samples on UMAP before integration")

DimPlot(object = first.exp.s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name") +
  ggtitle("All 4 samples on UMAP before integration") # no integration needed, biological difference, not batch effects!

#####----------------------------------------------------------------------#####
# EOF
#####----------------------------------------------------------------------#####
