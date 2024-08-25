func01_prepare_data <- function(umap, clone1, clone2, sample1, sample2){
  tmpdf1 <- subset(umap, umap$CTaa == clone1 & name == sample1)
  tmpdf2 <- subset(umap, umap$CTaa == clone2 & name == sample2)
  
  row.names(tmpdf1) <- NULL
  row.names(tmpdf2) <- NULL
  
  count.cell.in.cluster1 <- tmpdf1 %>% subset(select = c(Barcode, seurat_clusters)) %>% column_to_rownames("Barcode") %>% table() %>% as.data.frame
  colnames(count.cell.in.cluster1) <- c("cluster", "count_1")
  
  count.cell.in.cluster2 <- tmpdf2 %>% subset(select = c(Barcode, seurat_clusters)) %>% column_to_rownames("Barcode") %>% table() %>% as.data.frame
  colnames(count.cell.in.cluster2) <- c("cluster", "count_2")
  
  return(list(df1 = count.cell.in.cluster1, 
              df2 = count.cell.in.cluster2))
}


func01_prepare_data_clonewise <- function(umap, clone1, clone2){
  tmpdf1 <- subset(umap, umap$CTaa == clone1)
  tmpdf2 <- subset(umap, umap$CTaa == clone2)
  
  row.names(tmpdf1) <- NULL
  row.names(tmpdf2) <- NULL
  
  tmpdf1 <- as.data.frame(tmpdf1)
  tmpdf2 <- as.data.frame(tmpdf2)
  
  count.cell.in.cluster1 <- tmpdf1 %>% subset(select = c(Barcode, seurat_clusters)) %>% as.data.frame()
  row.names(count.cell.in.cluster1) <- NULL
  
  count.cell.in.cluster1 <- count.cell.in.cluster1 %>% column_to_rownames("Barcode") %>% table() %>% as.data.frame
  colnames(count.cell.in.cluster1) <- c("cluster", "count_1")
  
  count.cell.in.cluster2 <- tmpdf2 %>% subset(select = c(Barcode, seurat_clusters)) %>% as.data.frame 
  row.names(count.cell.in.cluster2) <- NULL
  
  count.cell.in.cluster2 <- count.cell.in.cluster2 %>% column_to_rownames("Barcode") %>% table() %>% as.data.frame
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


func03_pipeline_01_02 <- function(umap, clone1, clone2, sample1, sample2, mode = "normal"){
  df_list <- func01_prepare_data(umap = umap, clone1 = clone1, clone2 = clone2, sample1 = sample1, sample2 = sample2)
  MHI <- func02_calculate_MHI(count.cell.in.cluster1 = df_list$df1, count.cell.in.cluster2 = df_list$df2)
  return(MHI)
}

func03_pipeline_01_02_clonewise <- function(umap, clone1, clone2){
  df_list <- func01_prepare_data_clonewise(umap = umap, clone1 = clone1, clone2 = clone2)
  MHI <- func02_calculate_MHI(count.cell.in.cluster1 = df_list$df1, count.cell.in.cluster2 = df_list$df2)
  return(MHI)
}
