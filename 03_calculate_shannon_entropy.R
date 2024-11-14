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
##### Calculate Shannon Entropy for all clones in each dataset
#####---------------------------------------------------------------------------------------------#####
entropydf <- list()
for (dataset.name in names(all.clonedf)){
  ##### for dataset1, since there are two different compartments of cells in the data, 
  ##### we calculate Shannon entropy for each compartment separately. 
  tmpdf <- all.clonedf[[dataset.name]]
  if (dataset.name == "Dataset1_GFP"){
    restricted_to_clusters <- c(3, 5, 6, 8)
  } else if (dataset.name == "Dataset1_CD45"){
    restricted_to_clusters <- c(0, 1, 2, 4, 9, 10, 11)
  } else {
    restricted_to_clusters <- c()
  }
  tmpdf <- tmpdf %>% rowwise() %>% 
    mutate(Shannon.entropy = calculate_shannon_entropy_for_a_clone(input.clone = clone, 
                                                             input.metadata = all.metadata[[dataset.name]],
                                                             restricted_to_clusters = restricted_to_clusters))      
  entropydf[[dataset.name]] <- tmpdf
  writexl::write_xlsx(tmpdf, file.path(path.to.03.output, sprintf("%s.Shannon_entropy.xlsx", dataset.name)))
}

