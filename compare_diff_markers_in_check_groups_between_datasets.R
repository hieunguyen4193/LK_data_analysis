gc()
rm(list = ls())
#####----------------------------------------------------------------------#####
##### LIBRARIES 
#####----------------------------------------------------------------------#####
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"

source(file.path(path.to.pipeline.src, "processes_src", "helper_functions.R"))
source(file.path(path.to.pipeline.src, "processes_src", "import_libraries.R"))

library(ggpubr)

# just run ONCE per each analysis!!!
# remove.packages("DOSE")
# remove.packages("GOSemSim")
# remove.packages("yulab.utils")
# remove.packages("clusterProfiler")
# remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
# remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
# install.packages("/home/hieunguyen/CRC1382/storage/offline_pkgs/org.Mm.eg.db_3.18.0.tar.gz", type = "sources", repos = NULL)
# install.packages("/home/hieunguyen/CRC1382/storage/offline_pkgs/org.Mm.eg.db_3.18.0.tar.gz", type = "sources", repos = NULL)
# BiocManager::install("org.Mm.eg.db", update = FALSE)
# install.packages("heatmaply")
library(clusterProfiler)
library(org.Mm.eg.db)

# outdir <- "/home/hieunguyen/CRC1382/outdir/LKopplin_OFFICIAL"
outdir <- "/media/hieunguyen/HD0/outdir/CRC1382/LKopplin_OFFICIAL"
path.to.input <- file.path(outdir, "check_cluster_seggregation")
path.to.save.output <- file.path(outdir, "pathway_analysis_cluster_seggregation")
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

raw.diff.markers <- hash()
diff.markers <- hash()
for (dataset.name in names(special.clusters)){
  all.genes <- row.names(all.s.obj[[dataset.name]])
  mito.patterns <- "^mt-|^MT-"
  ribo.patterns <- "^Rpl|^Rps|^RPL|^RPS"
  ribo.genes <- to_vec( for (item in all.genes) if (str_detect(item, ribo.patterns) == TRUE) item)
  mito.genes <- to_vec( for (item in all.genes) if (str_detect(item, mito.patterns) == TRUE) item)
  raw.diff.markers[[dataset.name]] <- readRDS(file.path(path.to.input, sprintf("%s_diff_markers.raw.rds", dataset.name)))
  diff.markers[[dataset.name]] <- readRDS(file.path(path.to.input, sprintf("%s_diff_markers.raw.rds", dataset.name))) %>%
    subset(p_val_adj <= 0.05)
}

if (file.exists(file.path(path.to.save.output, "finished.rds")) == FALSE){
  ora.GOdf <- hash()
  ora.KEGGdf <- hash()
  GSEA.GOdf <- hash()
  GSEA.KEGGdf <- hash()
  
  ora.GO <- hash()
  ora.KEGG <- hash()
  GSEA.GO <- hash()
  GSEA.KEGG <- hash()
  
  #####----------------------------------------------------------------------#####
  ##### ORA
  #####----------------------------------------------------------------------#####
  
  shared.genes <- intersect(intersect(diff.markers$Dataset2$Gene, diff.markers$Dataset3$Gene), diff.markers$Dataset4$Gene)
  mito.shared.genes <- to_vec( for (item in shared.genes) if (str_detect(item, mito.patterns) == TRUE) item)
  ribo.shared.genes <- to_vec( for (item in shared.genes) if (str_detect(item, ribo.patterns) == TRUE) item)
  
  genes.to.do.ORA <- list(
    shared_genes_3_datasets = shared.genes,
    shared_genes_3_datasets_no_MT = setdiff(shared.genes, mito.shared.genes),
    shared_genes_3_datasets_no_Ribo = setdiff(shared.genes, ribo.shared.genes),
    shared_genes_3_datasets_no_MT_no_Ribo =  setdiff(shared.genes, c(ribo.shared.genes, mito.shared.genes)),
    dataset2 = diff.markers$Dataset2$Gene,
    dataset3 = diff.markers$Dataset3$Gene,
    dataset4 = diff.markers$Dataset4$Gene,
    dataset2_no_MT_no_Ribo = setdiff(diff.markers$Dataset2$Gene, c(mito.genes, ribo.genes)),
    dataset3_no_MT_no_Ribo = setdiff(diff.markers$Dataset3$Gene, c(mito.genes, ribo.genes)),
    dataset4_no_MT_no_Ribo = setdiff(diff.markers$Dataset4$Gene, c(mito.genes, ribo.genes))
  )
  
  for (savename in names(genes.to.do.ORA)){
    print(sprintf("Working on %s", savename))
    input.gene.list <- genes.to.do.ORA[[savename]]
    path.to.output <- file.path(path.to.save.output, "ORA", savename)
    dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
    input.gene.list.entrez <- bitr(input.gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
    ora.GO[[savename]] <- enrichGO(gene = input.gene.list,
                                   OrgDb = org.Mm.eg.db,
                                   ont = "ALL",
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.05,
                                   readable = TRUE,
                                   keyType = "SYMBOL",
                                   pAdjustMethod = "BH")
    
    ora.GOdf[[savename]] <- data.frame(ora.GO[[savename]])
    writexl::write_xlsx(ora.GOdf[[savename]], file.path(path.to.output, sprintf("ORA_GO_full_result.%s.xlsx", savename)))
    
    ora.KEGG[[savename]] <-  enrichKEGG(gene = input.gene.list.entrez,
                                        organism     = 'mmu',
                                        pvalueCutoff = 0.05, 
                                        qvalueCutoff = 0.05)
    ora.KEGGdf[[savename]] <- data.frame(ora.KEGG[[savename]])
    writexl::write_xlsx(ora.KEGGdf[[savename]], file.path(path.to.output, sprintf("ORA_KEGG_full_result.%s.xlsx", savename)))
    
    dotplot.ora.go <- dotplot(ora.GO[[savename]], size = "GeneRatio", showCategory=20) + ggtitle("Dotplot for Over-representation analysis")  
    dotplot.ora.kegg <- dotplot(ora.KEGG[[savename]], size = "GeneRatio", showCategory=20) + ggtitle("Dotplot for Gene set enrichment analysis") 
    
    ggsave(plot = dotplot.ora.go, filename = sprintf("ORA_GO_dotplot_geneRatio.svg"), path = path.to.output, device = "svg", width = 20, height = 10, dpi = 300)
    ggsave(plot = dotplot.ora.kegg, filename = sprintf("ORA_KEGG_dotplot_geneRatio.svg"), path = path.to.output, device = "svg", width = 20, height = 10, dpi = 300)
  }
  
  
  #####----------------------------------------------------------------------#####
  ##### GSEA
  #####----------------------------------------------------------------------#####
  
  for (dataset.name in c("Dataset2", "Dataset3", "Dataset4")){
    input.diffdf <- raw.diff.markers[[dataset.name]]
    
    path.to.output <- file.path(path.to.save.output, "GSEA", dataset.name)
    dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
    
    input.diffdf <- input.diffdf %>% arrange(desc(avg_log2FC))
    input.gene.list <- input.diffdf$avg_log2FC
    names(input.gene.list) <- input.diffdf$Gene
    
    GSEA.GO[[dataset.name]] <- gseGO(geneList = input.gene.list,
                                     OrgD = org.Mm.eg.db,
                                     ont = "ALL",
                                     minGSSize = 100,
                                     maxGSSize = 500,
                                     pvalueCutoff = 0.05,
                                     verbose = TRUE,
                                     keyType = "SYMBOL", seed = TRUE)
    GSEA.GOdf[[dataset.name]] <- data.frame(GSEA.GO[[dataset.name]])
    writexl::write_xlsx(GSEA.GOdf[[dataset.name]], file.path(path.to.output, sprintf("GSEA_GO_full_result.%s.xlsx", dataset.name)))
    plot.df <- GSEA.GOdf[[dataset.name]] %>% 
      as.data.frame() %>% 
      subset(select = c(Description, NES))
    
    top10.up.NES <- subset(plot.df, plot.df$NES >= 0) %>% arrange(desc(NES)) %>% head(10)
    top10.down.NES <- subset(plot.df, plot.df$NES < 0) %>% arrange(desc(NES)) %>% tail(10)
    top20.plotdf<- rbind(top10.up.NES, top10.down.NES)
    p <- ggplot(data = top20.plotdf, aes(y = reorder(Description, NES), x = NES)) + geom_bar(stat = "identity")    
    ggsave(plot = p, filename = sprintf("%s_GSEA_GO_NES_plot.svg", dataset.name), path = path.to.output, device = "svg", width = 14, height = 10, dpi = 300)
    
    tmp.full.list <- input.diffdf %>% arrange(desc(avg_log2FC))
    convertdf <- bitr(tmp.full.list$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
    input.gene.list <- tmp.full.list$avg_log2FC
    names(input.gene.list) <- convertdf$ENTREZID
    
    GSEA.KEGG[[dataset.name]] <- gseKEGG(geneList = input.gene.list,
                                         organism = "mmu",
                                         minGSSize = 100,
                                         maxGSSize = 500,
                                         pvalueCutoff = 0.05,
                                         verbose = TRUE, seed = TRUE, keyType = "kegg")
    GSEA.KEGGdf[[dataset.name]] <- data.frame(GSEA.KEGG[[dataset.name]])
    writexl::write_xlsx(GSEA.KEGGdf[[dataset.name]], file.path(path.to.output, sprintf("GSEA_KEGG_full_result.%s.xlsx", dataset.name)))
    plot.df <- GSEA.KEGG[[dataset.name]] %>% 
      as.data.frame() %>% 
      subset(select = c(Description, NES))
    
    top10.up.NES <- subset(plot.df, plot.df$NES >= 0) %>% arrange(desc(NES)) %>% head(10)
    top10.down.NES <- subset(plot.df, plot.df$NES < 0) %>% arrange(desc(NES)) %>% tail(10)
    top20.plotdf<- rbind(top10.up.NES, top10.down.NES)
    p <- ggplot(data = top20.plotdf, aes(y = reorder(Description, NES), x = NES)) + geom_bar(stat = "identity")    
    ggsave(plot = p, filename = sprintf("%s_GSEA_KEGG_NES_plot.svg", dataset.name), path = path.to.output, device = "svg", width = 14, height = 10, dpi = 300)
  }
  
  saveRDS(ora.GO, file.path(path.to.save.output, "ora.GO.rds"))
  saveRDS(ora.KEGG, file.path(path.to.save.output, "ora.KEGG.rds"))
  saveRDS(GSEA.GO, file.path(path.to.save.output, "GSEA.GO.rds"))
  saveRDS(GSEA.KEGG, file.path(path.to.save.output, "GSEA.KEGG.rds"))
  
  saveRDS(ora.GOdf, file.path(path.to.save.output, "ora.GOdf.rds"))
  saveRDS(ora.KEGGdf, file.path(path.to.save.output, "ora.KEGGdf.rds"))
  saveRDS(GSEA.GOdf, file.path(path.to.save.output, "GSEA.GOdf.rds"))
  saveRDS(GSEA.KEGGdf, file.path(path.to.save.output, "GSEA.KEGGdf.rds"))
  
  saveRDS(data.frame(stats = c("finished_all_analysis")), file.path(path.to.save.output, "finished.rds"))
} else {
  ora.GO <- readRDS(file.path(path.to.save.output, "ora.GO.rds"))
  ora.KEGG <- readRDS(file.path(path.to.save.output, "ora.KEGG.rds"))
  GSEA.GO <- readRDS(file.path(path.to.save.output, "GSEA.GO.rds"))
  GSEA.KEGG <- readRDS(file.path(path.to.save.output, "GSEA.KEGG.rds"))
  
  ora.GOdf <- readRDS(file.path(path.to.save.output, "ora.GOdf.rds"))
  ora.KEGGdf <- readRDS(file.path(path.to.save.output, "ora.KEGGdf.rds"))
  GSEA.GOdf <- readRDS(file.path(path.to.save.output, "GSEA.GOdf.rds"))
  GSEA.KEGGdf <- readRDS(file.path(path.to.save.output, "GSEA.KEGGdf.rds"))
}

#####----------------------------------------------------------------------#####
##### COMPARE PATHWAY ANALYSIS RESULTS: GSEA RESULTS
#####----------------------------------------------------------------------#####
GSEA.GO.all3.ids <- intersect(intersect(GSEA.GOdf$Dataset2$ID, GSEA.GOdf$Dataset3$ID), GSEA.GOdf$Dataset4$ID)
GSEA.KEGG.all3.ids <- intersect(intersect(GSEA.KEGGdf$Dataset2$ID, GSEA.KEGGdf$Dataset3$ID), GSEA.KEGGdf$Dataset4$ID)

GSEA.GO.all3.pathways <- subset(GSEA.GOdf$Dataset2, GSEA.GOdf$Dataset2$ID %in% GSEA.GO.all3.ids)[, c("ONTOLOGY", "ID", "Description")] 
colnames(GSEA.GO.all3.pathways) <- c("ONTOLOGY", "GO_ID", "Description")
GSEA.GO.all3.pathways <- GSEA.GO.all3.pathways %>%
  rowwise() %>%
  mutate(NES_Dataset2 = subset(GSEA.GOdf$Dataset2, GSEA.GOdf$Dataset2$ID == GO_ID)$NES) %>%
  mutate(NES_Dataset3 = subset(GSEA.GOdf$Dataset3, GSEA.GOdf$Dataset3$ID == GO_ID)$NES) %>%
  mutate(NES_Dataset4 = subset(GSEA.GOdf$Dataset4, GSEA.GOdf$Dataset4$ID == GO_ID)$NES) %>%
  mutate(abs_avg_NES = abs(NES_Dataset2 + NES_Dataset3 + NES_Dataset4)/3) %>% 
  arrange(desc(abs_avg_NES))
writexl::write_xlsx(GSEA.GO.all3.pathways, file.path(path.to.save.output, "GSEA_GO_all_3_datasets_pathways.xlsx"))

GSEA.KEGG.all3.pathways <- subset(GSEA.KEGGdf$Dataset2, GSEA.KEGGdf$Dataset2$ID %in% GSEA.KEGG.all3.ids)[, c("ID", "Description")] 
colnames(GSEA.KEGG.all3.pathways) <- c("KEGG_ID", "Description")
GSEA.KEGG.all3.pathways <- GSEA.KEGG.all3.pathways %>%
  rowwise() %>%
  mutate(NES_Dataset2 = subset(GSEA.KEGGdf$Dataset2, GSEA.KEGGdf$Dataset2$ID == KEGG_ID)$NES) %>%
  mutate(NES_Dataset3 = subset(GSEA.KEGGdf$Dataset3, GSEA.KEGGdf$Dataset3$ID == KEGG_ID)$NES) %>%
  mutate(NES_Dataset4 = subset(GSEA.KEGGdf$Dataset4, GSEA.KEGGdf$Dataset4$ID == KEGG_ID)$NES) %>%
  mutate(abs_avg_NES = abs(NES_Dataset2 + NES_Dataset3 + NES_Dataset4)/3) %>% 
  arrange(desc(abs_avg_NES))
writexl::write_xlsx(GSEA.KEGG.all3.pathways, file.path(path.to.save.output, "GSEA_KEGG_all_3_datasets_pathways.xlsx"))



