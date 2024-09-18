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
path.to.03.output <- file.path(path.to.main.output, "03_output")

dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

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

#####---------------------------------------------------------------------------------------------#####
##### helper functions: Calculate MHI between 2 clones in a dataset
#####---------------------------------------------------------------------------------------------#####

func01_prepare_data <- function(umap, clone1, clone2, sample1, sample2){
  tmpdf1 <- subset(umap, umap$CTaa == clone1 & name == sample1)
  tmpdf2 <- subset(umap, umap$CTaa == clone2 & name == sample2)
  
  if (nrow(tmpdf1)  == 0 | nrow(tmpdf2) == 0){
    stop("No cells found for the given input clone and sample")
  }
  row.names(tmpdf1) <- NULL
  row.names(tmpdf2) <- NULL
  
  count.cell.in.cluster1 <- tmpdf1 %>% subset(select = c(barcode, seurat_clusters)) %>% column_to_rownames("barcode") %>% table() %>% as.data.frame
  colnames(count.cell.in.cluster1) <- c("cluster", "count_1")
  
  count.cell.in.cluster2 <- tmpdf2 %>% subset(select = c(barcode, seurat_clusters)) %>% column_to_rownames("barcode") %>% table() %>% as.data.frame
  colnames(count.cell.in.cluster2) <- c("cluster", "count_2")
  
  return(list(df1 = count.cell.in.cluster1, 
              df2 = count.cell.in.cluster2))
}
func02_calculate_MHI <- function(count.cell.in.cluster1, count.cell.in.cluster2){
  count.mat <- merge(count.cell.in.cluster1, count.cell.in.cluster2, by.x = "cluster", by.y = "cluster")
  
  X <- sum(count.mat$count_1)
  Y <- sum(count.mat$count_2)
  
  
  count.mat <- count.mat %>% rowwise() %>% mutate(xy = count_1*count_2)
  count.mat <- count.mat %>% rowwise() %>% mutate(x2 = count_1^2)
  count.mat <- count.mat %>% rowwise() %>% mutate(y2 = count_2^2)
  
  MHI_numerator <- 2 * sum(count.mat$xy)
  MHI_denominator <- X * Y * ((sum(count.mat$x2)/X^2) + (sum(count.mat$y2)/Y^2))
  
  MHI <- MHI_numerator/MHI_denominator
  return(MHI)
}

func03_pipeline_01_02 <- function(umap, clone1, clone2, sample1, sample2){
  df_list <- func01_prepare_data(umap = umap, clone1 = clone1, clone2 = clone2, sample1 = sample1, sample2 = sample2)
  MHI <- func02_calculate_MHI(count.cell.in.cluster1 = df_list$df1, count.cell.in.cluster2 = df_list$df2)
  return(MHI)
}

##### Example: Calculate the MHI index for "distribution among clusters" clones
##### CAAENTGNYKYVF_CASSLWGGDAETLYF in sample m368 and m369 in Dataset4.
func03_pipeline_01_02(umap = all.metadata[["Dataset4"]],
                      sample1 = "m368",
                      sample2 = "m369",
                      clone1 = "CAAENTGNYKYVF_CASSLWGGDAETLYF",
                      clone2 = "CAAENTGNYKYVF_CASSLWGGDAETLYF"
                      )

