showCloneHighlight_modified <- function (apotc_ggplot, clonotype, color_each = TRUE, default_color = "#808080", 
                                         scale_bg = 1, fill_legend = TRUE, input.colors = NULL) 
{
  num_seqs <- length(clonotype)
  sequence_index_map <- hash::hash(clonotype, 1:num_seqs)
  highlighted_ggplot_data <- get_ggplot_data(apotc_ggplot)
  clone_color_vector <- gen_clone_color_vector(color_each, 
                                               clonotype, highlighted_ggplot_data)
  num_matches <- 0
  for (i in seq_len(nrow(highlighted_ggplot_data))) {
    curr_clone <- highlighted_ggplot_data$clonotype[i]
    curr_seq_index <- sequence_index_map[[curr_clone]]
    if (is.null(curr_seq_index)) {
      if (!is.null(default_color)) {
        highlighted_ggplot_data$color[i] <- default_color
      }
      highlighted_ggplot_data$color[i] <- scale_hex_brightness(highlighted_ggplot_data$color[i], 
                                                               scale_bg)
      next
    }
    num_matches <- num_matches + 1
    if (identical(color_each, FALSE)) 
      next
    highlighted_ggplot_data$color[i] <- clone_color_vector[curr_seq_index]
  }
  if (num_matches < num_seqs) {
    warning(ifelse(num_matches == 0, "all", "some"), " ", 
            "input sequences didn't match any clone")
  }
  apotc_ggplot <- set_ggplot_data(apotc_ggplot, highlighted_ggplot_data)
  if (identical(color_each, FALSE) || !fill_legend) 
    return(apotc_ggplot)
  # suppressMessages(apotc_ggplot + ggplot2::scale_fill_identity(guide = "legend", 
  #                                                              name = "clonotype", 
  #                                                              labels = clonotype, 
  #                                                              breaks = clone_color_vector))
  suppressMessages(apotc_ggplot + ggplot2::scale_fill_manual  (guide = "legend",
                                                               name = "clonotype",
                                                               labels = clonotype,
                                                               breaks = clone_color_vector,
                                                               values = input.colors))
}
