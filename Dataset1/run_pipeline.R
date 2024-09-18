gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)

#####----------------------------------------------------------------------#####
# PIPELINE CONFIGURATIONS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/LKopplin_data"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/LK_data_analysis"
PROJECT <- "1stExp_Kopplin"

path.to.main.input <- file.path(path.to.storage, PROJECT)

path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.project.src <- "/media/hieunguyen/HNSD01/src/LK_data_analysis/Dataset1"
path.to.downstream.rmd <- file.path(path.to.project.src, "downstream_analysis.Rmd")

path.to.VDJ.input <- file.path(path.to.main.input, "TCR")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline"

source(file.path(path.to.pipeline.src, "scRNA_VDJ_pipeline", "main_VDJ_pipeline.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "processes_src", "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline", "scRNA_GEX_pipeline.R"))


if ("DescTools" %in% installed.packages() == FALSE){
  install.packages("https://cran.r-project.org/src/contrib/Archive/DescTools/DescTools_0.99.50.tar.gz", repos = NULL, type = "source")
}
#####----------------------------------------------------------------------#####
# VDJ pipeline
#####----------------------------------------------------------------------#####
summarize_vdj_data(path.to.VDJ.input, 
                   path.to.VDJ.output, 
                   PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, 
                   T_or_B = "T")


#####----------------------------------------------------------------------#####
# GEX Pipeline
#####----------------------------------------------------------------------#####
for (analysis.round in c("1st", "2nd")){
  path.to.VDJ.output <- file.path(path.to.main.output, "/VDJ_output")
  
  stage_lst <- hash()
  
  stage_lst[["GFPexp1_m1"]] <- c(GFPexp1_m1 = "GFP")
  stage_lst[["GFPexp1_m2"]] <- c(GFPexp1_m2 = "GFP")
  
  stage_lst[["CD45exp1_m1"]] <- c(CD45exp1_m1 = "CD45")
  stage_lst[["CD45exp1_m2"]] <- c(CD45exp1_m2 = "CD45")
  
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
  
  for (sample.id in names(stage_lst)){
    #####------------------------------------------------------------------#####
    # PATHS CONFIGURATIONS FOR EACH sample.id-RUN
    #####------------------------------------------------------------------#####
    path2input <- file.path(path.to.main.input, "GEX", sprintf("GEX_single_%s", sample.id))
    
    path.to.sample.id.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
    dir.create(path.to.sample.id.output, showWarnings = FALSE)
    
    path.to.saved.html <- file.path(path.to.sample.id.output, sprintf("%s_%s_round", sample.id, analysis.round), "downstream_analysis_report")
    dir.create(path.to.saved.html, showWarnings = FALSE)
    
    path.to.anno.contigs <- file.path(path.to.VDJ.output, sprintf("annotated_contigs_clonaltype_%s.csv", sample.id))
    
    path.to.count.clonaltype <- file.path(path.to.VDJ.output, sprintf("count_clonaltype_%s.csv", sample.id))
    
    if (analysis.round == "2nd"){
      path.to.filtered.barcodes <- file.path(path.to.main.output, "1st_round")
      filtered.barcodes <- readRDS(file.path(path.to.filtered.barcodes, sprintf("%s_1st_round/s9_output/remove_barcodes/%s_1st_round_remove_barcodes.rds", sample.id, sample.id)))
    } else if (analysis.round == "3rd"){
      path.to.filtered.barcodes.1 <- file.path(path.to.main.output, "1st_round")
      filtered.barcodes.1 <- readRDS(file.path(path.to.filtered.barcodes.1, sprintf("%s_1st_round/s9_output/remove_barcodes/%s_1st_round_remove_barcodes.rds", sample.id, sample.id)))
      
      path.to.filtered.barcodes.2 <- file.path(path.to.main.output, "2nd_round")
      filtered.barcodes.2 <- readRDS(file.path(path.to.filtered.barcodes.2, sprintf("%s_2nd_round/s9_output/remove_barcodes/%s_2nd_round_remove_barcodes.rds", sample.id, sample.id)))
      
      filtered.barcodes <- c(filtered.barcodes.1, filtered.barcodes.2)
    } else {
      filtered.barcodes <- NULL
    }
    #####------------------------------------------------------------------#####
    # MAIN GEX PIPELINE RUN
    #####------------------------------------------------------------------#####
    s.obj <- run_pipeline_GEX(path2src=path2src,
                              path2input=path2input,
                              path.to.logfile.dir=file.path(path.to.sample.id.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"),
                              stage_lst=stage_lst[[sample.id]],
                              path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                              MINCELLS=MINCELLS,
                              MINGENES=MINGENES,
                              PROJECT=PROJECT,
                              remove_doublet=remove_doublet,
                              save.RDS=save.RDS,
                              path.to.output=file.path(path.to.sample.id.output, sprintf("%s_%s_round", sample.id, analysis.round)),
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
                              sw =  sw,
                              my_random_seed = my_random_seed)
    
    path.to.save.html <- file.path(outdir, "html_output", PROJECT)
    dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
    rmarkdown::render(path.to.downstream.rmd,
                      params = list(
                        sample.id = sample.id,
                        outdir =  outdir,
                        PROJECT = PROJECT,
                        analysis.round = analysis.round
                      ),
                      output_file = sprintf("downstream_analysis_%s_%s_%s.html", PROJECT, sample.id, analysis.round),
                      output_dir = path.to.save.html)
  }
}

writeLines(capture.output(sessionInfo()), file.path(path.to.main.output, sprintf("%s_sessionInfo.txt", PROJECT)))