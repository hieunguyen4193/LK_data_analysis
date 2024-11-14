gc()
rm(ls())

##### input args
s.obj <- "/path/to/seurat.rds"
all.cells <- colnames(s.obj)

entropydf <- "/path/to/clonedf_with_entropy.csv"
path.to.save.output <- "/path/to/save/output"

entropydf <- subset(entropydf, entropydf$count >= 10)
all.clone.sizes <- sort(unique(entropydf$count))
dir.create(file.path(path.to.save.output, "validation_Shannon_entropy"), showWarnings = FALSE, recursive = TRUE)

all.shannon.entropy <- subset(entropydf, select = c(clone, Shannon.entropy, count))

all.validationdf <- data.frame()
for (num.sampling.cells in all.clone.sizes){
  if (file.exists(file.path(path.to.save.output, "validation_Shannon_entropy", 
                            sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", 
                                    num.sampling.cells, 
                                    selected.dataset))) == FALSE){
    all.sampling.entropies <- c()
    
    for (count in seq(1, 1000)){
      sampling.cells <- sample(all.cells, num.sampling.cells)
      tmp.sampling.shannon.entropy <- calculate_shannon_entropy(cell.list = sampling.cells, 
                                                                input.metadata = all.metadata[[selected.dataset]], 
                                                                restricted_to_clusters = c())
      all.sampling.entropies <- c(all.sampling.entropies, tmp.sampling.shannon.entropy)
    }
    
    samplingdf <- data.frame(data = all.sampling.entropies)
    colnames(samplingdf) <- c("Shannon.entropy")
    samplingdf$clone <- to_vec( for (item in seq(1, 1000)) sprintf("sampling_%s", item))
    samplingdf <- samplingdf[c("clone", "Shannon.entropy")]
    samplingdf$case <- "Sampling"
    samplingdf$count <- num.sampling.cells
    all.shannon.entropy$case <- "real_clone"
    
    validationdf <- rbind(all.shannon.entropy, samplingdf)
    writexl::write_xlsx(validationdf, file.path(path.to.save.output, "validation_Shannon_entropy", sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", num.sampling.cells,selected.dataset)))
  } else {
    validationdf <- readxl::read_xlsx(file.path(path.to.save.output, "validation_Shannon_entropy", sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", num.sampling.cells, selected.dataset)))
  }
  all.validationdf <- rbind(all.validationdf, validationdf)
}

all.validationdf <- all.validationdf %>% 
  rowwise() %>%
  mutate(plot.count = ifelse(case == "real_clone", "real_clone", sprintf("sampling_%s", count)))

all.validationdf$plot.count <- factor(
  all.validationdf$plot.count, 
  levels = c("real_clone", 
             to_vec(for (i in sort(all.clone.sizes)) sprintf("sampling_%s", i)))
)
p.boxplots <- all.validationdf %>% ggplot(aes(x = plot.count, y = Shannon.entropy)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

##### double y axis plots
coeff <- 1000
p.entropy_clonesize <- entropydf %>% ggplot(aes(x = reorder(clone, -count))) + 
  geom_point(aes(y = Shannon.entropy)) +  ylim(0, 1) +
  geom_bar(aes(y = count/coeff), stat = "identity") +  
  theme(axis.text.x = element_blank()) + xlab("Clone") +  
  scale_y_continuous(
    name = "Shannon entropy",
    sec.axis = sec_axis(~.*coeff, name="Number of cells") 
  ) 