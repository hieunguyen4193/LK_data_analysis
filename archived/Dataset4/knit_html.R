gc()
rm(list = ls())

PROJECT <- "230316_Kopplin"
dataset.name <- "Dataset4"

path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))

path.to.project.src <- file.path("/home/hieunguyen/CRC1382/src_2023/LKopplin/20231101_OFFICIAL", dataset.name)

outdir <- "/media/hieunguyen/HNSD_MBPro/CRC1382/outdir/LKopplin_OFFICIAL"

path.to.save.html <- file.path(outdir, "html", PROJECT)
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.rmd.filenames <- c(
  "01_GEX_data_analysis_dataset4.Rmd",
  "02_VDJ_data_analysis_dataset4_m366_m367.Rmd",
  "02_VDJ_data_analysis_dataset4_m366_m368.Rmd",
  "02_VDJ_data_analysis_dataset4_m366_m369.Rmd",
  "02_VDJ_data_analysis_dataset4_m367_m368.Rmd",
  "02_VDJ_data_analysis_dataset4_m367_m369.Rmd",
  "02_VDJ_data_analysis_dataset4_m368_m369.Rmd",
  "03_MHI_data_analysis_m366_vs_m367.Rmd",
  "03_MHI_data_analysis_m368_vs_m369.Rmd",
  "04_Shannon_entropy_sampling_validation_m366_vs_m367.Rmd",
  "04_Shannon_entropy_sampling_validation_m366_vs_m368.Rmd",
  "04_Shannon_entropy_sampling_validation_m366_vs_m369.Rmd",
  "04_Shannon_entropy_sampling_validation_m367_vs_m368.Rmd",
  "04_Shannon_entropy_sampling_validation_m367_vs_m369.Rmd",
  "04_Shannon_entropy_sampling_validation_m368_vs_m369.Rmd"
)

for (rmd.filename in all.rmd.filenames){
  path.to.Rmd.file <- file.path(path.to.project.src, rmd.filename)
  save.HTML.filename <- str_replace(rmd.filename, ".Rmd", ".html")
  if (file.exists(file.path(path.to.save.html, save.HTML.filename)) == FALSE){
    rmarkdown::render(input = path.to.Rmd.file, 
                      output_file = save.HTML.filename,
                      output_dir = path.to.save.html)      
  } else {
    print(sprintf("File %s is already existed at %s", save.HTML.filename, path.to.save.html))
  }

}


