gc()
rm(list = ls())


if ("ggpubr" %in% installed.packages() == FALSE){
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-0.tar.gz", type = "source", repos = NULL)
  install.packages("ggpubr")
  install.packages("car")
}
                 
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
library(ggpubr)
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

dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

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

##### change the input dataset here
selected.dataset <- "Dataset4"

s.obj <- all.s.obj[[selected.dataset]]
clonedf <- all.clonedf[[selected.dataset]]
entropydf <- all.entropydf[[selected.dataset]]
entropydf <- subset(entropydf, entropydf$count >= 10)
all.clone.sizes <- sort(unique(entropydf$count))
dir.create(file.path(path.to.05.output, "validation_Shannon_entropy"), showWarnings = FALSE, recursive = TRUE)

all.shannon.entropy <- subset(entropydf, select = c(clone, Shannon.entropy, count))

all.validationdf <- data.frame()
for (num.sampling.cells in all.clone.sizes){
  if (file.exists(file.path(path.to.05.output, "validation_Shannon_entropy", 
                            sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", 
                                    num.sampling.cells, 
                                    selected.dataset))) == FALSE){
    all.cells <- colnames(s.obj)
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
    writexl::write_xlsx(validationdf, file.path(path.to.05.output, "validation_Shannon_entropy", sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", num.sampling.cells,selected.dataset)))
  } else {
    validationdf <- readxl::read_xlsx(file.path(path.to.05.output, "validation_Shannon_entropy", sprintf("validation_Shannon_entropy_sampling_%s.%s.xlsx", num.sampling.cells, selected.dataset)))
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

#####----------------------------------------------------------------------------------------#####
##### double y axis plots
#####----------------------------------------------------------------------------------------#####
coeff <- 1000
p.entropy_clonesize <- entropydf %>% ggplot(aes(x = reorder(clone, -count))) + 
  geom_point(aes(y = Shannon.entropy)) +  ylim(0, 1) +
  geom_bar(aes(y = count/coeff), stat = "identity") +  
  theme(axis.text.x = element_blank()) + xlab("Clone") +  
  scale_y_continuous(
    name = "Shannon entropy",
    sec.axis = sec_axis(~.*coeff, name="Number of cells") 
  ) 

##### Add new plot, 05.12.2024

# Observed vs Simulated dotplot -------------------------------------------

plt <- ggplot(all.validationdf, aes(x = count, y = Shannon.entropy)) +
  geom_jitter(aes(fill = plot.count), shape = 21, size = 5, color = "black", width = 0.2, height = 0) + # Add jitter
  stat_compare_means(
    comparisons = list(c('0', '15')),
    paired = FALSE, 
    method = "t.test", 
    label = "p.signif", 
    position = position_dodge(width = 0.8), 
    hide.ns = FALSE,
    step.increase = 0,
    tip.length = 0,
    size = 10, # Adjust size of the significance labels
    bracket.size = 0.5
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, family = 'Helvetica'),
    strip.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 0, 10),
    legend.title = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1)) +
  scale_x_discrete(
    breaks = unique(all.validationdf$count),
    labels = function(x) ifelse(x == '0', "n=433", x) # Remove label at 0
  ) +
  scale_fill_manual(values = c('Simulated' = 'grey', 'Observed' = 'orange')) +
  labs(
    y = "Shannon entropy",
    x = "Selected simulated clone sizes"
  )


# Clone size entropy plot -------------------------------------------------

coeff <- 400
plt <- ggplot(real_entropy, aes(x = rank, y = size_of_clones / coeff, group = 1)) +
  geom_line(color = "black", alpha = 0.7, size = 1) +
  geom_area(fill = "lightgrey") +
  geom_point(aes(y = shannon_entropy), color = "darkblue", alpha = 1, size = 2) +
  labs(
    y = "Shannon entropy",
    x = "Clones ranked by size"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 20, family = "Helvetica"),
    axis.ticks.x = element_line(size = 0.5),
    axis.ticks.y.left = element_line(color = "darkblue", size = 0.5), 
    axis.ticks.y.right = element_line(color = "black", size = 0.5), 
    axis.text.y.left = element_text(color = "darkblue"),
    axis.title.y.left = element_text(color = "darkblue"), 
    axis.text.y.right = element_text(color = "black"),    
    axis.title.y.right = element_text(color = "black"),   
    panel.grid = element_blank(),
    axis.line.y.left = element_line(color = "black", linewidth = 0.5), 
    axis.line.y.right = element_line(color = "black", linewidth = 0.5), 
    axis.line.x = element_line(color = "black", linewidth = 0.5)
  ) +
  scale_x_discrete(
    breaks = as.character(seq(0, length(real_entropy$rank), by = 50)),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.1)),
    name = "Shannon entropy",
    sec.axis = sec_axis(~.*coeff, name = "Number of cells")
  )
plt

