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
path.to.05.output <- file.path(path.to.main.output, "05_output")
path.to.06.output <- file.path(path.to.main.output, "06_output")
path.to.09.output <- file.path(path.to.main.output, "09_output")

dir.create(path.to.09.output, showWarnings = FALSE, recursive = TRUE)

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

all.entropydf <- hash()
for (dataset.name in names(all.clonedf)){
  all.entropydf[[dataset.name]] <- readxl::read_excel(file.path(path.to.03.output, sprintf("%s.Shannon_entropy.xlsx", dataset.name)))
}

all.s.obj <- list()
for (dataset.name in names(dataset.names)){
  all.s.obj[[dataset.name]] <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))
}


for (input.dataset in names(all.s.obj)){
  s.obj <- all.s.obj[[input.dataset]]
  
  clonedf <- s.obj@meta.data %>% subset(select = c(seurat_clusters, CTaa)) %>%
    subset(is.na(CTaa) == FALSE) 
  
  clone.countdf <- table(clonedf$seurat_clusters, clonedf$CTaa) %>% as.data.frame() %>%
    pivot_wider(names_from = "Var2", values_from = "Freq") %>% 
    column_to_rownames("Var1") %>% 
    t() %>% 
    as.data.frame()
  colnames(clone.countdf) <- to_vec(
    for (i in colnames(clone.countdf)){
      sprintf("cluster%s", i)
    }
  )
  all.clusters <- clone.countdf %>% colnames()
  clone.countdf <- clone.countdf %>% rownames_to_column("CloneID")
  
  mhidf <- data.frame(cluster = all.clusters)
  
  calculate_mhi <- function(tmpdf, input.cluster1, input.cluster2){
    X <- tmpdf[[input.cluster1]] %>% sum()
    Y <- tmpdf[[input.cluster2]] %>% sum()
    
    x <- unlist(lapply(
      tmpdf$CloneID %>% unique(),
      function(x){
        return(subset(tmpdf, tmpdf$CloneID == x)[[input.cluster1]])  
      }
    ))
    
    y <- unlist(lapply(
      tmpdf$CloneID %>% unique(),
      function(x){
        return(subset(tmpdf, tmpdf$CloneID == x)[[input.cluster2]])
      }
    ))
    
    nom <- 2 * sum(x * y)
    det <- (sum(x^2)/X^2) + (sum(y^2)/Y^2)
    mhi <- nom/(X*Y * det)
    return(mhi)
  }
  
  for (input.cluster1 in mhidf$cluster){
    mhidf[[input.cluster1]] <- unlist(lapply(
      mhidf$cluster, function(x){
        calculate_mhi(tmpdf = clone.countdf, input.cluster1 = input.cluster1, input.cluster2 = x)
      }
    ))
  }
  
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)] <- NA
    diag(cormat) <- NA
    return(cormat)
  }
  
  mhidf.upper <- get_upper_tri(mhidf %>% column_to_rownames("cluster")) %>% rownames_to_column("cluster")
  
  mhidf.upper.pivot <- mhidf.upper %>% pivot_longer(!cluster, names_to = "cluster2", values_to = "MHI") 
  mhidf.upper.pivot$cluster <- factor(mhidf.upper.pivot$cluster, levels = all.clusters)
  mhidf.upper.pivot$cluster2 <- factor(mhidf.upper.pivot$cluster2, levels = all.clusters)
  
  mhi.p <- mhidf.upper.pivot %>% 
    subset(is.na(MHI) == FALSE) %>%
    ggplot(aes(x = cluster, y = cluster2, fill = MHI)) +
    geom_tile(color = "black", size=0.8)+
    theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                     size = 12, hjust = 1))+
    theme(axis.text.y = element_text(size = 12))+
    coord_fixed()+
    labs(x = NULL, y = NULL) +
    theme(panel.background = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_blank()) + 
    scale_fill_gradient(high = "red", low = "white")
  ggsave(plot = mhi.p, filename = sprintf("MHI_%s_between_clusters.svg", input.dataset), path = path.to.09.output, width = 14, height = 10, dpi = 300, device = "svg")
}

