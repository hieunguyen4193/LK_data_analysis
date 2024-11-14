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
#####----------------------------------------------------------------------#####
##### helper functions: Calculate Shannon entropy to measure "diversity" of cells/clones over clusters
#####----------------------------------------------------------------------#####
calculate_shannon_entropy <- function(cell.list, input.metadata, restricted_to_clusters = NA){
  if (length(restricted_to_clusters) != 0){
    input.metadata <- subset(input.metadata, input.metadata$seurat_clusters %in% restricted_to_clusters)
  } 
  ##### by default, we assume that the column "seurat_clusters" contains the active cluster numbers.
  N <- length(unique(input.metadata$seurat_clusters))
  count_clone_in_cluster <- table(subset(input.metadata, input.metadata$barcode %in% cell.list) %>% 
                                    subset(select = c(seurat_clusters))) %>%  as.data.frame()
  count_clone_in_cluster <- count_clone_in_cluster %>% rowwise %>% mutate(Prob = Freq / sum(count_clone_in_cluster$Freq)) %>%
    subset(Prob != 0)
  shannon_entropy <- -sum(count_clone_in_cluster$Prob * log2(count_clone_in_cluster$Prob))/log2(N)
  return(shannon_entropy)
}

# To calculate Shannon entropy for a clone, we extract its list of barcodes (cell.list) and feed the function calculate_shannon_entropy_sampling_cells
calculate_shannon_entropy_for_a_clone <- function(input.clone, input.metadata, restricted_to_clusters = NA){
  cell.list <- subset(input.metadata, input.metadata$CTaa == input.clone)$barcode
  shannon_entropy <- calculate_shannon_entropy(cell.list, input.metadata, restricted_to_clusters)
  return(shannon_entropy)
}
