gc()
rm(list = ls())
my_random_seed <- 42

set.seed(my_random_seed)
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "211227_Kopplin"

source("/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")
path.to.storage <- "/home/hieunguyen/CRC1382/storage"
outdir <- "/media/hieunguyen/HNSD_MBPro/CRC1382/outdir/LKopplin_OFFICIAL"

path.to.main.input <- file.path(path.to.storage, PROJECT)
path.to.VDJ.input <- file.path(path.to.main.input, "VDJ")

path.to.main.output <- file.path(outdir, PROJECT)
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")

dir.create(path.to.main.output, showWarnings = FALSE)
dir.create(path.to.VDJ.output, showWarnings = FALSE, recursive = TRUE)

summarize_vdj_data(path.to.VDJ.input, 
                   path.to.VDJ.output, 
                   PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, 
                   T_or_B = "T")

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/LKopplin/20231101_OFFICIAL/Dataset2"
path.to.downstream.rmd <- file.path(path.to.project.src, "preliminary_downstream_analysis")

source(file.path(path2src, "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

for (analysis.round in c("1st", "2nd")){
  
  path2input <- file.path(path.to.main.input, "GEX_single_samples")
  
  path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")
  
  # _____stage lst for single sample_____
  stage_lst <- hash()
  
  stage_lst[["GFP_m1"]] <- c(GFP_m1 = "GFP")
  stage_lst[["GFP_m2"]] <- c(GFP_m2 = "GFP")
  stage_lst[["GFP_m3"]] <- c(GFP_m3 = "GFP")
  
  MINCELLS  <- 0
  MINGENES  <- 0
  
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
             s8 = "off",
             s8a = "on",
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
  
  
  filter.thresholds <- list(nFeatureRNAfloor = NULL,
                            nFeatureRNAceiling = NULL,
                            nCountRNAfloor = NULL, 
                            nCountRNAceiling = NULL,
                            pct_mitofloor = NULL, 
                            pct_mitoceiling = 10,
                            pct_ribofloor = NULL, 
                            pct_riboceiling = NULL,
                            ambientRNA_thres = 0.5)
  
  remove_doublet <- FALSE
  path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")
  
  
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  
  for (sample in names(stage_lst)){
    path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
    dir.create(path.to.output, showWarnings = FALSE)
    path.to.saved.html <- file.path(path.to.output, sprintf("%s_%s_round", sample, analysis.round), "downstream_analysis_report")
    
    path.to.anno.contigs <- file.path(path.to.VDJ.output, sprintf("annotated_contigs_clonaltype_%s.csv", sample))
    path.to.count.clonaltype <- file.path(path.to.VDJ.output, sprintf("count_clonaltype_%s.csv", sample))
    if (analysis.round == "2nd"){
      path.to.filtered.barcodes <- file.path(path.to.main.output, "1st_round")
      filtered.barcodes <- readRDS(file.path(path.to.filtered.barcodes, sprintf("%s_1st_round/remove_barcodes/%s_1st_round_remove_barcodes.rds", sample, sample)))
    } else if (analysis.round == "3rd"){
      path.to.filtered.barcodes.1 <- file.path(path.to.main.output, "1st_round")
      filtered.barcodes.1 <- readRDS(file.path(path.to.filtered.barcodes.1, sprintf("%s_1st_round/remove_barcodes/%s_1st_round_remove_barcodes.rds", sample, sample)))
      
      path.to.filtered.barcodes.2 <- file.path(path.to.main.output, "2nd_round")
      filtered.barcodes.2 <- readRDS(file.path(path.to.filtered.barcodes.2, sprintf("%s_2nd_round/remove_barcodes/%s_2nd_round_remove_barcodes.rds", sample, sample)))
      
      filtered.barcodes <- c(filtered.barcodes.1, filtered.barcodes.2)
    } else if (analysis.round == "4th") {
      path.to.filtered.barcodes.1 <- file.path(path.to.main.output, "1st_round")
      filtered.barcodes.1 <- readRDS(file.path(path.to.filtered.barcodes.1, 
                                               sprintf("%s_1st_round/remove_barcodes/%s_1st_round_remove_barcodes.rds", sample, sample)))
      
      path.to.filtered.barcodes.2 <- file.path(path.to.main.output, "2nd_round")
      filtered.barcodes.2 <- readRDS(file.path(path.to.filtered.barcodes.2, 
                                               sprintf("%s_2nd_round/remove_barcodes/%s_2nd_round_remove_barcodes.rds", sample, sample)))
      
      path.to.filtered.barcodes.3 <- file.path(path.to.main.output, "3rd_round")
      filtered.barcodes.3 <- readRDS(file.path(path.to.filtered.barcodes.3, 
                                               sprintf("%s_3rd_round/remove_barcodes/%s_3rd_round_remove_barcodes.rds", sample, sample)))
      
      filtered.barcodes <- c(filtered.barcodes.1, filtered.barcodes.2, filtered.barcodes.3)
      
    } else {
      filtered.barcodes <- NULL
    }
    s.obj <- run_pipeline_GEX(path2src=path2src,
                              path2input=file.path(path2input, sprintf("GEX_single_%s", sample)),
                              path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample, analysis.round), "logs"),
                              stage_lst=stage_lst[[sample]],
                              path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                              MINCELLS=MINCELLS,
                              MINGENES=MINGENES,
                              PROJECT=PROJECT,
                              remove_doublet=remove_doublet,
                              save.RDS=save.RDS,
                              path.to.output=file.path(path.to.output, sprintf("%s_%s_round", sample, analysis.round)),
                              rerun=rerun,
                              DE.test="wilcox",
                              num.PCA=num.PCA,
                              num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                              num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                              use.sctransform=FALSE,
                              filtered.barcodes=filtered.barcodes,
                              filter.thresholds=filter.thresholds,
                              path.to.anno.contigs=path.to.anno.contigs,
                              path.to.count.clonaltype=path.to.count.clonaltype,
                              input.method = "filterTRAB",
                              my_random_seed = my_random_seed)
    
    dir.create(path.to.saved.html, showWarnings = FALSE)
    
    rmarkdown::render(file.path(path.to.downstream.rmd,
                                sprintf("%s_round", analysis.round),
                                sprintf("%s_round_downstream_analysis_for_%s.Rmd", analysis.round, sample)),
                      params = list(input.datadir = path.to.output),
                      output_dir = path.to.saved.html)
  }
}

#### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))

