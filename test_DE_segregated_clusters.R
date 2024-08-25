gc()
rm(list = ls())
#####----------------------------------------------------------------------#####
##### LIBRARIES 
#####----------------------------------------------------------------------#####
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))

library(ggpubr)

outdir <- "/home/hieunguyen/CRC1382/outdir/LKopplin_OFFICIAL"
path.to.save.output <- file.path(outdir, "check_cluster_seggregation")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

##### PATHS
# Dataset paths
dataset.path <- list(Dataset1 = file.path(outdir, "1stExp_Kopplin"),
                     Dataset2 = file.path(outdir, "211227_Kopplin"),
                     Dataset3 = file.path(outdir, "230215_Kopplin_Pabst_added_NC_000001_merged_zcat_m330_m331_remove_c9"),
                     Dataset4 = file.path(outdir, "230316_Kopplin"))

all.s.obj <- list(
  Dataset1 = readRDS(file.path(dataset.path$Dataset1, "data_analysis/02_output/sobj_with_clusterRes_1.rds")),
  Dataset2 = readRDS(file.path(dataset.path$Dataset2, "data_analysis/01_output/merged_all_second_exp_dataset.annotated.filteredCD45.integrated.rds")),
  Dataset3 = readRDS(file.path(dataset.path$Dataset3, "1st_round/pct_mito_10_1/data_analysis/02_output/230215_Kopplin.seurat.obj.addedShannonEntropy.rds")),
  Dataset4 = readRDS(file.path(dataset.path$Dataset4, "1st_round/pct_mito_10_1/data_analysis/01_output/230316_Kopplin.seurat.obj.removed.14_15_16.addedVDJ.integrated.rds"))
)

special.clusters <- list(
  Dataset2 = c(2),
  Dataset3 = c(3, 5, 8),
  Dataset4 = c(2, 3, 10)
)

for (dataset.name in c("Dataset2", "Dataset3", "Dataset4")){
  if (file.exists(file.path(path.to.save.output, sprintf("%s_cluster_composition.svg", dataset.name))) == FALSE){
    #####--------------------------------------------------------------------#####
    ##### SAMPLE COMPOSITION IN EACH CLUSTER PLOTS
    #####--------------------------------------------------------------------#####
    
    s.obj <- all.s.obj[[dataset.name]]
    
    umap.plot <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "name")
    num.sample.in.dataset <- length(unique(s.obj$name))
    sample.names <- unique(s.obj$name)
    
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>% subset(select = c(barcode, seurat_clusters))
    
    embeddings <- s.obj@reductions$INTE_UMAP@cell.embeddings %>% as.data.frame() %>% rownames_to_column("barcode")
    embeddings <- merge(embeddings, meta.data, by.x = "barcode", by.y = "barcode")
    
    get_x_y_for_cluster <- function(cluster.id){
      tmp.embeddings <- subset(embeddings, embeddings$seurat_clusters == cluster.id)
      embedding.colnames <- to_vec( for (item in colnames(tmp.embeddings)) if (grepl("UMAP", item) == TRUE | grepl("umap", item) == TRUE) item)
      x <- mean(tmp.embeddings[[embedding.colnames[[1]]]])
      y <- mean(tmp.embeddings[[embedding.colnames[[2]]]])
      return(c(x, y))
    }
    
    x <- to_vec(for (item in seq(0, length(unique(s.obj@meta.data$seurat_clusters))- 1)) get_x_y_for_cluster(item)[[1]])
    y <- to_vec(for (item in seq(0, length(unique(s.obj@meta.data$seurat_clusters))- 1)) get_x_y_for_cluster(item)[[2]])
    
    library(RColorBrewer)
    if (num.sample.in.dataset >= 3){
      sample_colors <- brewer.pal(num.sample.in.dataset, "Accent")    
    } else if (num.sample.in.dataset == 2) {
      sample_colors <- c("#f2968f", "#8fcff2")
    }
    
    
    num.clusters <- length(unique(s.obj$seurat_clusters))
    colors <- c(hue_pal()(num.clusters), sample_colors)
    names(colors) <- c(to_vec(for (item in seq(0, num.clusters - 1)) sprintf("%s", item)), sample.names)
    
    meta.data <- s.obj@meta.data %>% subset(select = c(name, seurat_clusters))
    
    count.cells.df <- t(table(meta.data)) %>% as.data.frame() %>% pivot_wider(names_from = "name", values_from = "Freq")
    count.cells.df$total <- rowSums(count.cells.df[, sample.names])
    
    
    for (sample.id in sample.names){
      count.cells.df[[sprintf("pct_%s_in_cluster", sample.id)]] <- count.cells.df[[sample.id]]/count.cells.df$total
    }
    
    count.cells.df$x <- x
    count.cells.df$y <- y
    
    bar.plot.raw.counts <- count.cells.df[, c("seurat_clusters", sample.names)] %>% 
      pivot_longer(!seurat_clusters, names_to = "sample", values_to = "count") %>%
      ggplot(aes(x = seurat_clusters, y = count, fill = sample)) + geom_bar(stat = "identity", position = "dodge") +
      scale_fill_manual(values = sample_colors)  
    
    bar.plot.pct <- count.cells.df[, c("seurat_clusters", to_vec(for (item in colnames(count.cells.df)) if(grepl("pct", item) == TRUE) item))] %>% 
      pivot_longer(!seurat_clusters, names_to = "sample", values_to = "count") %>%
      ggplot(aes(x = seurat_clusters, y = count, fill = sample)) + geom_bar(stat = "identity") +
      scale_fill_manual(values = sample_colors)  
    
    dimplot.p <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE) + guides(fill="none") +
      geom_scatterpie(aes(x = x, y = y), data = count.cells.df[, c("seurat_clusters", sample.names, "x", "y")], cols = sample.names, alpha=0.8) + 
      scale_fill_manual(values = colors[sort(names(colors))]) 
    
    final.plot <- ggarrange(bar.plot.raw.counts,
                            bar.plot.pct,
                            umap.plot, 
                            dimplot.p,
                            ncol = 2, nrow = 2)
    ggsave(plot = final.plot, path = file.path(path.to.save.output), filename = sprintf("%s_cluster_composition.svg", dataset.name), device = "svg", width = 28, height = 20, dpi = 300)
  }
}

#####--------------------------------------------------------------------#####
##### SAMPLE COMPOSITION IN EACH CLUSTER PLOTS
#####--------------------------------------------------------------------#####
for (dataset.name in c("Dataset2", "Dataset3", "Dataset4")){
  s.obj <- all.s.obj[[dataset.name]]
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    mutate(check.cluster = ifelse(seurat_clusters %in% special.clusters[[dataset.name]], "1", "0")) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$check.cluster, col.name = "check.cluster")
  
  if (file.exists(file.path(path.to.save.output, sprintf("%s_diff_markers.filtered.rds", dataset.name))) == FALSE){
    diff.markers <- FindMarkers(object = s.obj, assay = "RNA", test.use = "wilcox", ident.1 = "1", ident.2 = "0", group.by = "check.cluster")
    diff.markers <- diff.markers %>%
      rownames_to_column("Gene") %>% 
      rowwise() %>%
      mutate(abs.logFC = abs(avg_log2FC))
    saveRDS(diff.markers, file.path(path.to.save.output, sprintf("%s_diff_markers.raw.rds", dataset.name)))
    diff.markers <- subset(diff.markers, diff.markers$p_val_adj <= 0.05 & diff.markers$abs.logFC >= 1) 
    saveRDS(diff.markers, file.path(path.to.save.output, sprintf("%s_diff_markers.filtered.rds", dataset.name)))  
  } else {
    diff.markers <- readRDS(file.path(path.to.save.output, sprintf("%s_diff_markers.filtered.rds", dataset.name)))
  }
  diff.markers <- diff.markers %>% arrange(desc(abs.logFC))
  violin.plot <- VlnPlot(object = s.obj, features= head(diff.markers, 16)$Gene, group.by = "check.cluster", pt.size = 0, assay = "RNA")
  writexl::write_xlsx(diff.markers, file.path(path.to.save.output, sprintf("%s_diff_markers.full.xlsx", dataset.name)))
  p <- FeaturePlot(object = s.obj, reduction = "INTE_UMAP", features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "log10GenesPerUMI"), ncol = 2)  
  normal.UMAP <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "check.cluster")
  ggsave(plot = normal.UMAP, filename = sprintf("%s_UMAP_reclustering.svg", dataset.name), path = file.path(path.to.save.output), dpi = 300, device = "svg", width = 14, height = 10)
  ggsave(plot = p, filename = sprintf("%s_QC_metrics.svg", dataset.name), path = file.path(path.to.save.output), dpi = 300, device = "svg", width = 28, height = 30)
  ggsave(plot = violin.plot, filename = sprintf("%s_violin_plot.svg", dataset.name), path = file.path(path.to.save.output), dpi = 300, device = "svg", width = 28, height = 20)
}


#####----------------------------------------------------------------------#####
##### MHI CALCULATION IN SELECTED CLUSTERS ONLY
#####----------------------------------------------------------------------#####
path.to.main.src.dir <- "/home/hieunguyen/CRC1382/src_2023/LKopplin/scRNAseq_VDJ_Kopplin_data"
path.to.rmd.file <- list(Dataset2 = file.path(path.to.main.src.dir, "Dataset2", "05_MHI_in_selected_clusters.Rmd"),
                         Dataset3 = file.path(path.to.main.src.dir, "Dataset3", "05_MHI_analysis_in_selected_clusters.Rmd"),
                         Dataset4_m366_vs_m367 = file.path(path.to.main.src.dir, "Dataset4", "05_MHI_data_analysis_m366_vs_m367.selected_clusters.Rmd"),
                         Dataset4_m368_vs_m369 = file.path(path.to.main.src.dir, "Dataset4", "05_MHI_data_analysis_m368_vs_m369.selected_clusters.Rmd"))

for (dataset.name in names(path.to.rmd.file)){
  path.to.save.HTML.dir <- file.path(path.to.save.output, "MHI_selected_clusters_html")
  dir.create(path.to.save.HTML.dir, showWarnings = FALSE, recursive = TRUE)
  save.HTML.filename <- sprintf("%s_MHI_in_selected_clusters.html", dataset.name)
  rmarkdown::render(input = path.to.rmd.file[[dataset.name]], 
                    output_file = save.HTML.filename,
                    output_dir = path.to.save.HTML.dir)  
}


