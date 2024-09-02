plot_clone_on_UMAP <- function(clone){
  m1.cell.names <- row.names(subset(slot(s.obj.2nd, "meta.data"),
                                    slot(s.obj.2nd, "meta.data")$CTaa == clone &
                                      slot(s.obj.2nd, "meta.data")$name == "GFP_m1"))
  m2.cell.names <- row.names(subset(slot(s.obj.2nd, "meta.data"),
                                    slot(s.obj.2nd, "meta.data")$CTaa == clone &
                                      slot(s.obj.2nd, "meta.data")$name == "GFP_m2"))
  m3.cell.names <- row.names(subset(slot(s.obj.2nd, "meta.data"),
                                    slot(s.obj.2nd, "meta.data")$CTaa == clone &
                                      slot(s.obj.2nd, "meta.data")$name == "GFP_m3"))
  cell.names = list(GFP_m1 = m1.cell.names, 
                    GFP_m2 = m2.cell.names,
                    GFP_m3 = m3.cell.names)
  
  p <- DimPlot(object = s.obj.2nd,
               cells.highlight = cell.names,
               cols.highlight = c("#f77281", "#3773db", "#3dad8c"),
               cols = "lightgray", order = TRUE,
               sizes.highlight = 2,
               label = TRUE, label.size = 8,
               label.box = FALSE, reduction = "INTE_UMAP", pt.size = 0.5) + 
    xlim(-7, 10) + ylim(-8, 8) +
    theme(plot.title = element_text(size = 14, face = "bold")) +
    ggtitle(sprintf("%s, \nm1: %s, m2: %s, m3: %s", clone, 
                    length(m1.cell.names),
                    length(m2.cell.names),
                    length(m3.cell.names))) 
  return(p)
}

