gc()
rm(list = ls())
my_random_seed <- 42
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "230316_Kopplin"

path.to.storage <- "/home/hieunguyen/CRC1382/storage"

outdir <- "/media/hieunguyen/HNSD_MBPro/CRC1382/outdir/LKopplin_OFFICIAL"

path.to.main.output <- file.path(outdir, PROJECT)

dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

source("/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")
path.to.VDJ.input <- file.path(path.to.storage, PROJECT, "VDJ")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")

summarize_vdj_data(path.to.VDJ.input, 
                   path.to.VDJ.output, 
                   PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, 
                   T_or_B = "T")

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

set.seed(my_random_seed)

source(file.path(path2src, "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

path2input <- file.path(path.to.storage, PROJECT, "GEX")

stage_lst <- list()

stage_lst[["m366"]] <- c( m366 = "m366")
stage_lst[["m367"]] <- c( m366 = "m367")
stage_lst[["m368"]] <- c( m366 = "m368")
stage_lst[["m369"]] <- c( m366 = "m369")

MINCELLS  <- 5
MINGENES  <- 50

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = FALSE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = TRUE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "off",
           s8 = "on",
           s8a = "off",
           s9 = "on")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)


pct_mitoceiling <- 10

input.filter.thresholds <- hash()
##### 1st round
input.filter.thresholds[["1st_round"]] <- list( nFeatureRNAfloor = NULL,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = NULL,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)

##### 2nd round
input.filter.thresholds[["2nd_round"]] <- list( nFeatureRNAfloor = 250,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = 2000,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)

##### 3rd round
input.filter.thresholds[["3rd_round"]] <- list( nFeatureRNAfloor = 2500,
                                                nFeatureRNAceiling = NULL,
                                                nCountRNAfloor = 5000,
                                                nCountRNAceiling = NULL,
                                                pct_mitofloor = NULL,
                                                pct_mitoceiling = pct_mitoceiling,
                                                pct_ribofloor = NULL,
                                                pct_riboceiling = NULL,
                                                ambientRNA_thres = 0.5)


remove_doublet <- TRUE
path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
num.dim.integration <- 25
num.dim.cluster <- 25

# remove_XY_genes <- c("Xist","Jpx","Ftx","Tsx","Cnbp2")
remove_XY_genes <- NULL
filtered.barcodes <- NULL

cluster.resolution <- 1

for (analysis.round in c("1st_round")){
  filter.thresholds <- input.filter.thresholds[[analysis.round]]

  path.to.anno.contigs <- NULL 
  path.to.count.clonaltype <- NULL
  
  path.to.output <- file.path(path.to.main.output, analysis.round, sprintf("pct_mito_%s_%s", pct_mitoceiling, cluster.resolution))
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- run_pipeline_GEX(path2src = path2src,
                            path2input = file.path(path2input),
                            path.to.logfile.dir = file.path(path.to.output, "logs"),
                            stage_lst = stage_lst,
                            path.to.10X.doublet.estimation = path.to.10X.doublet.estimation,
                            MINCELLS = MINCELLS,
                            MINGENES = MINGENES,
                            PROJECT = PROJECT,
                            remove_doublet = remove_doublet,
                            save.RDS = save.RDS,
                            path.to.output = path.to.output,
                            rerun = rerun, 
                            DE.test = "wilcox",
                            num.PCA = num.PCA,
                            num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                            use.sctransform = FALSE,
                            filtered.barcodes = filtered.barcodes,
                            filter.thresholds = filter.thresholds,
                            input.metho = "filterTRAB",
                            path.to.anno.contigs = path.to.anno.contigs,
                            path.to.count.clonaltype = path.to.count.clonaltype,
                            cluster.resolution = cluster.resolution,
                            mode_cell_cycle_scoring = "gene_name",
                            my_random_seed = my_random_seed,
                            umap.method = "uwot",
                            ambientRNA_thres = 0.5,
                            path.to.renv = NULL,
                            features_to_regressOut=NULL,
                            regressOut_mode="alternative",
                            num.dim.integration=num.dim.integration,
                            genes.to.not.run.PCA = NULL,
                            inte_pca_reduction_name = "INTE_PCA",
                            inte_umap_reduction_name = "INTE_UMAP",
                            pca_reduction_name = NULL,
                            umap_reduction_name = NULL,
                            path.to.s3a.source = NULL)
}

writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))  


# EOF