gc()
rm(list = ls())

PROJECT <- "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9"
dataset.name <- "Dataset3"

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))

path.to.project.src <- file.path("/home/hieunguyen/CRC1382/src_2023/LKopplin/20231101_OFFICIAL", dataset.name)

outdir <- "/media/hieunguyen/HNSD_MBPro/CRC1382/outdir/LKopplin_OFFICIAL"

path.to.save.html <- file.path(outdir, "html", PROJECT)
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.rmd.filenames <- c(
  "02_GEX_downstream_analysis_after_removing_cluster_9.Rmd",
  "03_MHI_analysis.Rmd",
  "04_Shannon_entropy_and_validation_sampling.Rmd"
)

for (rmd.filename in all.rmd.filenames){
  path.to.Rmd.file <- file.path(path.to.project.src, rmd.filename)
  save.HTML.filename <- str_replace(rmd.filename, ".Rmd", ".html")
  rmarkdown::render(input = path.to.Rmd.file, 
                    output_file = save.HTML.filename,
                    output_dir = path.to.save.html)  
}


