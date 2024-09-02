#####----------------------------------------------------------------------#####
# Assign status for the clonal type of the second dataset
#####----------------------------------------------------------------------#####
assign_status_for_clonaltype_2nd <- function(clone, count.clonaltype, summary.count.clonaltype, thres = 10){
  group1 <- count.clonaltype$GFP_m1$CTaa
  group2 <- count.clonaltype$GFP_m2$CTaa
  group3 <- count.clonaltype$GFP_m3$CTaa
  tmp <- summary.count.clonaltype  %>% column_to_rownames("CTaa")
  tmp <- subset(tmp, row.names(tmp) == clone)
  if (rowSums(tmp) < thres){
    status <- "excluded"
    return(status)
  } else {
    if ((clone %in% group1) & (clone %in% group2) & (clone %in% group3) ){
      # check if the clonal type is in all 3 samples. 
      # condition: sum(counts) >= thres
      status <- "in_all_3_samples"
      return(status)
    } else if ( ( (clone %in% group1) & (clone %in% group2) ) | 
                ( (clone %in% group1) & (clone %in% group3) ) |
                ( (clone %in% group2) & (clone %in% group3) )){
      # else: check if the clonal type is in 2 samples
      # condition: sum(counts) >= thres
      status <- "in_2_samples"
      return(status)
    } else {
      # else: check if the clonal type is in only 1 sample
      # counts >= thres
      status <- "unique_in_1_sample"
      return(status)
    }  
  }
  
}



##### higher Shannon entropy --> more diverse
calculate_shannon_entropy <- function(clone, s.obj){
  N <- length(unique(s.obj$seurat_clusters))
  count_clone_in_cluster <- table(subset(s.obj@meta.data, s.obj@meta.data$CTaa == clone) %>% 
                                    subset(select = c(CTaa, seurat_clusters))) %>%  as.data.frame()
  count_clone_in_cluster <- count_clone_in_cluster %>% rowwise %>% mutate(Prob = Freq / sum(count_clone_in_cluster$Freq)) %>%
    subset(Prob != 0)
  shannon_entropy <- -sum(count_clone_in_cluster$Prob * log2(count_clone_in_cluster$Prob))/log2(N)
  return(shannon_entropy)
}


plot_a_list_of_clonal_types_2nd_dataset <- function(list_of_clonal_types, s.obj, sample = NULL){
  do.call(patchwork::wrap_plots, lapply(list_of_clonal_types, function(x) {
    chosen.ctaa <- x
    m1.cell.names <- row.names(subset(slot(s.obj, "meta.data"),
                                      slot(s.obj, "meta.data")$CTaa == chosen.ctaa &
                                        slot(s.obj, "meta.data")$name == "GFP_m1"))
    m2.cell.names <- row.names(subset(slot(s.obj, "meta.data"),
                                      slot(s.obj, "meta.data")$CTaa == chosen.ctaa &
                                        slot(s.obj, "meta.data")$name == "GFP_m2"))
    m3.cell.names <- row.names(subset(slot(s.obj, "meta.data"),
                                      slot(s.obj, "meta.data")$CTaa == chosen.ctaa &
                                        slot(s.obj, "meta.data")$name == "GFP_m3"))
    cell.names = list(GFP_m1 = m1.cell.names, 
                      GFP_m2 = m2.cell.names,
                      GFP_m3 = m3.cell.names)
    
    p <- DimPlot(object = s.obj,
                 cells.highlight = cell.names,
                 cols.highlight = c("#f77281", "#3773db", "#3dad8c"),
                 cols = "lightgray", order = TRUE,
                 sizes.highlight = 1,
                 label = TRUE,
                 label.box = FALSE, reduction = "INTE_UMAP") +
      theme(plot.title = element_text(size = 8, face = "bold")) +
      ggtitle(sprintf("%s, \nm1: %s, m2: %s, m3: %s", chosen.ctaa, 
                      length(m1.cell.names),
                      length(m2.cell.names),
                      length(m3.cell.names))) 
    return(p)
  })) -> subplot_all_top_20_clonaltypes
  return(subplot_all_top_20_clonaltypes)
}

`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# FUNCTION TO CREATE INTERACTIVE DATA TABLE IN RMARKDOWN REPORTS
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}
