require(limma)
require(parallel)
require(ggplot2)
require(ggpubr)
require(future.apply)
require(car)
require(ridge)
require(e1071)
require(preprocessCore)
require(foreach)  
require(doParallel)
require(data.table)
require(ggalluvial)
require(pheatmap)
require(gridExtra)
require(reshape2)
require(circlize)

#### main analysis ####
dir.create("figure_reorgnized/Fig6")
outpdf <- "figure_reorgnized/Fig6/immune_content.pdf"
pFilter <- 0.05

exp <- expr_df[["gse32472"]] %>% as.matrix()
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
data <- data[rowMeans(data) > 0, ]

out <- rbind(ID = colnames(data), data)
out <- cbind(ID = rownames(out), out)

fwrite(out, file = "figure_reorgnized/Fig6/uniq_symbol.txt", sep = "\t", 
       quote = FALSE, col.names = FALSE, row.names = TRUE)

# Run main CIBERSORT code to get immune cell infiltration results
source("cibersort/CIBERSORT_multithread.R", encoding = "utf-8")
Sys.time()
results <- system.time(CIBERSORT("cibersort/ref.txt", 
                                 "figure_reorgnized/Fig6/uniq_symbol.txt", perm = 1000, QN = TRUE))
Sys.time()

# Read and filter CIBERSORT results
immune <- read.table("figure_reorgnized/Fig6/CIBERSORT_Results.txt", sep = "\t", 
                     header = TRUE, row.names = 1, check.names = FALSE)
immune_sig <- immune[immune[, "P-value"] < pFilter, ]
immune_sig <- as.matrix(immune_sig[, 1:(ncol(immune_sig) - 3)])

# Visualization of immune cell content
plot_data <- t(immune_sig)
col <- rainbow(nrow(plot_data), s = 0.75, v = 0.7)

# ##### stack barplot R #####
# pdf(outpdf, height = 10, width = 20)
# par(las = 1, mar = c(8, 5, 4, 16), mgp = c(3, 0.1, 0), cex.axis = 1.5)
# 
# a1 <- barplot(plot_data, col = col, yaxt = "n", ylab = "Relative Percent", xaxt = "n", cex.lab = 1.8)
# a2 <- axis(2, tick = FALSE, labels = FALSE)
# axis(2, a2, paste0(a2 * 100, "%"))
# axis(1, a1, labels = FALSE)
# 
# par(srt = 60, xpd = TRUE)
# text(a1, -0.03, colnames(plot_data), adj = 1, cex = 0.18)
# par(srt = 0)
# 
# ytick2 <- cumsum(plot_data[, ncol(plot_data)])
# ytick1 <- c(0, ytick2[-length(ytick2)])
# legend(par('usr')[2] * 0.98, par('usr')[4], legend = rownames(plot_data), 
#        col = col, pch = 15, bty = "n", cex = 1.3)
# dev.off()

# I prefer ggplot2

##### Fig6A stack barplot #####
d3_category <- pal_d3("category20b")(20)

d3_pals <- c(d3_category, "grey", "steelblue")

plot_data %>% 
  as.data.frame() %>% 
  rownames_to_column("immune") %>% 
  pivot_longer(-immune, names_to = "gse_id", values_to = "percentage") %>% 
  left_join(meta_df[["gse32472"]], by = "gse_id") %>% 
  # group_by(immune) %>%
  # mutate(median_score = median(percentage, na.rm = TRUE)) %>%
  # ungroup() %>%
  # mutate(immune = fct_reorder(immune, median_score, .desc = TRUE)) %>%
  mutate(immune = factor(immune, levels = rownames(plot_data))) %>% 
  ggplot(aes(gse_id, percentage, fill = immune)) +
  geom_col(width = 1, color = "black", size = 0.1) +
  facet_grid(~BPD_status, space = "free", scales = "free_x", switch = "x") +
  labs(
    x = NULL,
    y = "Relative Percentage (%)",
    title = NULL, 
    fill = "Immune Cell Types"
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12), 
    axis.title = element_text(size = 16), 
    strip.background = element_blank(),
    strip.placement = "outside", 
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.75, "line"),
    legend.box = "vertical",
    legend.position = "right", 
    strip.text = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1))+
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) + 
  scale_fill_manual(values = col)

ggsave("figure_reorgnized/Fig6/Fig6A_barplot_ggplot_BPD_status.pdf", height = 10, width = 20)

#### hub gene - immune infiltration correlation ####
# Convert immune score names to readable format with consistent spacing
format_immune_names <- function(x) {
  x <- gsub("\\.", " ", x)  # Replace dots with spaces
  x <- gsub("\\s+", " ", x)  # Replace multiple spaces with single space
  x <- trimws(x)  # Remove leading and trailing whitespace
  
  # Capitalize key words
  x <- sub("cells", "Cells", x, ignore.case = TRUE)
  x <- sub("memory", "Memory", x, ignore.case = TRUE)
  x <- sub("naive", "Naive", x, ignore.case = TRUE)
  x <- sub("activated", "Activated", x, ignore.case = TRUE)
  x <- sub("resting", "Resting", x, ignore.case = TRUE)
  x <- sub("helper", "Helper", x, ignore.case = TRUE)
  x <- sub("regulatory", "Regulatory", x, ignore.case = TRUE)
  x <- sub("follicular", "Follicular", x, ignore.case = TRUE)
  x <- sub("dendritic", "Dendritic", x, ignore.case = TRUE)
  x <- sub("mast", "Mast", x, ignore.case = TRUE)
  x <- sub("macrophages", "Macrophages", x, ignore.case = TRUE)
  x <- sub("monocytes", "Monocytes", x, ignore.case = TRUE)
  x <- sub("neutrophils", "Neutrophils", x, ignore.case = TRUE)
  x <- sub("eosinophils", "Eosinophils", x, ignore.case = TRUE)
  
  return(x)
}

gene_interest <- overlapping_genes[["cell_death"]] %>% 
  unlist() %>% as.vector()

genes_expr <- expr_df[["gse32472"]] %>% 
  t() %>% as.data.frame() %>% 
  select(any_of(gene_interest)) %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]]) %>% 
  filter(BPD_status == "BPD") %>% 
  select(-gse_id, -BPD_status, -BPD_type)

immune_scores <- immune_sig %>% 
  as.data.frame() %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]]) %>% 
  filter(BPD_status == "BPD") %>% 
  select(-gse_id, -BPD_status, -BPD_type)

# Check and fix column names to avoid issues with NA column names
colnames(genes_expr) <- make.names(colnames(genes_expr), unique = TRUE)
colnames(immune_scores) <- make.names(colnames(immune_scores), unique = TRUE)

# Calculate correlation matrix and p-values
correlation_matrix <- cor(as.matrix(immune_scores), as.matrix(genes_expr), use = "pairwise.complete.obs")
p_values_matrix <- matrix(NA, nrow = ncol(immune_scores), ncol = ncol(genes_expr),
                          dimnames = list(colnames(immune_scores), colnames(genes_expr)))

for (i in 1:ncol(immune_scores)) {
  for (j in 1:ncol(genes_expr)) {
    test <- cor.test(immune_scores[[i]], genes_expr[[j]], use = "pairwise.complete.obs")
    p_values_matrix[i, j] <- test$p.value
  }
}

# Replace any NA/NaN/Inf values
correlation_matrix[is.na(correlation_matrix) | is.infinite(correlation_matrix)] <- 0
p_values_matrix[is.na(p_values_matrix) | is.infinite(p_values_matrix)] <- 1

# Melt and combine correlation and p-values
correlation_melted <- melt(correlation_matrix)
p_values_melted <- melt(p_values_matrix)

combined_data <- correlation_melted %>% 
  rename(correlation = value) %>% 
  left_join(p_values_melted, by = c("Var1", "Var2")) %>% 
  rename(p_value = value) %>%
  mutate(label = sprintf("%.2f\n(%.2g)", correlation, p_value))

# Reshape for heatmap annotation
annotation_matrix <- acast(combined_data, Var1 ~ Var2, value.var = "label")

# Format immune cell names
formatted_immune_names <- format_immune_names(colnames(immune_scores))

##### FigS1 Full heatmap of correlation #####
dir.create("figure_reorgnized/Supplementary_Figure")

plot <- pheatmap(correlation_matrix,
                 color = colorRampPalette(c("blue", "white", "red"))(100), 
                 cluster_rows = TRUE, 
                 cluster_cols = TRUE, 
                 treeheight_row = 0, 
                 treeheight_col = 0, 
                 main = "", 
                 fontsize_row = 12, 
                 fontsize_col = 12, 
                 display_numbers = annotation_matrix,
                 number_color = "black", 
                 labels_row = formatted_immune_names,
                 labels_col = colnames(genes_expr), 
                 angle_col = 90, 
                 legend = TRUE,
                 legend_breaks = seq(-1, 1, by = 0.2), 
                 legend_labels = sprintf("%.1f", seq(-1, 1, by = 0.2)), 
                 clustering_method = "complete")
pdf("figure_reorgnized/Supplementary_Figure/FigS1_gene_immune_full_heatmap.pdf", height = 15, width = 20)
print(plot)
dev.off()

##### Fig6D hub gene vs immune traits correlation #####
selected_genes <- c("THBS1", "ALOX5", "ACSL1")
correlation_subset <- correlation_matrix[, selected_genes, drop = FALSE]
p_values_subset <- p_values_matrix[, selected_genes, drop = FALSE]

row_names <- format_immune_names(rownames(correlation_subset))
rownames(correlation_subset) <- row_names
rownames(p_values_subset) <- row_names

sig_matrix <- matrix(
  sapply(p_values_subset, function(p) {
    if (p < 0.001) return("\n***")
    else if (p < 0.01) return("\n**")
    else if (p < 0.05) return("\n*")
    else return("")
  }), 
  nrow = nrow(p_values_subset),
  dimnames = dimnames(p_values_subset)
)

text_annotation <- matrix(
  sprintf("%.2f%s", correlation_subset, sig_matrix),
  nrow = nrow(correlation_subset),
  dimnames = dimnames(correlation_subset)
)

pdf("figure_reorgnized/Fig6/Fig6D_gene_immune_complexheatmap_no_dendro.pdf", height = 7.5, width = 5)
Heatmap(correlation_subset,
        name = "Correlation",
        col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        row_names_side = "left",
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(text_annotation[i, j], x, y, gp = gpar(fontsize = 10))
        },
        row_names_gp = gpar(fontsize = 12),
        column_names_gp = gpar(fontsize = 12),
        column_names_rot = 45,
        heatmap_legend_param = list(
          title = "Correlation",
          at = seq(-1, 1, by = 0.2),
          labels = sprintf("%.1f", seq(-1, 1, by = 0.2))
        ))
dev.off()

#### lm ####
plot_data <- data.frame(
  THBS1 = as.numeric(expr_df[["gse32472"]]["THBS1",]),
  ALOX5 = as.numeric(expr_df[["gse32472"]]["ALOX5",]),
  ACSL1 = as.numeric(expr_df[["gse32472"]]["ACSL1",]),
  Neutrophils = as.numeric(immune_sig[, "Neutrophils"]),
  T.cells.CD8 = as.numeric(immune_sig[, "T cells CD8"])
)

# function
create_correlation_plot <- function(cell_type, gene, data, color) {
  orig_cell_name <- if(cell_type == "T.cells.CD8") "T.cells.CD8" else "Neutrophils"
  r_value <- correlation_matrix[orig_cell_name, gene]
  p_value <- p_values_matrix[orig_cell_name, gene]
  
  p <- ggplot(data, aes(x = .data[[gene]], y = .data[[cell_type]])) +
    geom_point(shape = 17, color = color, size = 2, alpha = 0.6) +
    geom_smooth(method = "lm", 
                color = "grey30", 
                fill = "grey70",
                formula = y ~ x) +
    labs(y = format_immune_names(cell_type),
         x = gene) +
    annotate("text", 
             x = min(data[[gene]]), 
             y = max(data[[cell_type]]),
             label = sprintf("R=%.2f, p=%.2g", r_value, p_value),
             hjust = 0, vjust = 1) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  p_with_marginal <- ggExtra::ggMarginal(p, 
                                         type = "histogram", 
                                         fill = color, 
                                         alpha = 0.6)
  
  return(p_with_marginal)
}

genes <- c("THBS1", "ALOX5", "ACSL1")
cell_types <- c("Neutrophils", "T.cells.CD8") 
colors <- c("darkred", "darkgrey")

dir.create("figure_reorgnized/Fig6/Fig6E2J/")

##### Fig6E-J #####
for(cell in cell_types) {
  for(gene in genes) {
    p <- create_correlation_plot(
      cell_type = cell,
      gene = gene,
      data = plot_data,
      color = if(cell == "Neutrophils") colors[1] else colors[2]
    )
    save_name <- gsub("\\.", "_", cell)
    ggsave(
      filename = file.path("figure_reorgnized/Fig6/Fig6E2J/", 
                           sprintf("correlation_%s_%s.pdf", 
                                   save_name,
                                   gene)),
      plot = p,
      width = 6,
      height = 6,
      device = "pdf"
    )
  }
}

#### Fig6B immune infiltration comparison ####
immune_cell_data <- t(immune_sig) %>%
  as.data.frame() %>%
  rownames_to_column("immune_cell") %>%
  pivot_longer(-immune_cell, names_to = "sample_id", values_to = "proportion") %>%
  left_join(meta_df[["gse32472"]], by = c("sample_id" = "gse_id"))

max_value <- max(immune_cell_data$proportion)
significance_y <- max_value

significance_tests <- immune_cell_data %>%
  group_by(immune_cell) %>%
  summarise(
    p_value = wilcox.test(proportion ~ BPD_status)$p.value
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    text_y = significance_y 
  )

# Create final plot with significance
ggplot(immune_cell_data, aes(x = immune_cell, y = proportion)) +
  geom_boxplot(aes(fill = BPD_status), outlier.shape = 16, outlier.size = 0.6) +
  geom_text(data = significance_tests, 
            aes(x = immune_cell, y = text_y, label = significance),
            inherit.aes = FALSE,
            size = 4) +
  scale_fill_manual(values = c("green3", "darkred")) +
  scale_y_continuous(limits = c(0, significance_y + 0.01)) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.position = "top", 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("figure_reorgnized/Fig6/Fig6B_immune_proportion_boxplot.pdf", height = 6, width = 12)

#### immune trait - PCD pattern correlation ####
immune_data <- immune_sig %>%
  as.data.frame() %>%
  rownames_to_column("gse_id")

create_gene_set <- function(genes, set_name) {
  GeneSet(genes,
          geneIdType = SymbolIdentifier(),
          setName = set_name)
}

gene_sets <- list()
for(name in names(death_genes)) {
  gene_sets[[name]] <- create_gene_set(death_genes[[name]], name)
}

gene_set_collection <- GeneSetCollection(gene_sets)
gmtFile <- "figure_reorgnized/Fig6/death_genes.gmt"
toGmt(gene_set_collection, gmtFile)

cell_death_gmt <- getGmt(gmtFile, 
                         collectionType = BroadCollection(category = "c3"), 
                         geneIdType = SymbolIdentifier())

mat <- as.matrix(expr_df[["gse32472"]]) %>% 
  avereps() %>% 
  normalizeBetweenArrays()

gsva_params <- gsvaParam(
  expr = mat,
  geneSets = cell_death_gmt,
  maxDiff = FALSE
)

gsva_results <- gsva(
  param = gsva_params,
  verbose = TRUE
)

output <- rbind(id = colnames(gsva_results), gsva_results)
write.table(output, file="figure_reorgnized/Fig6/output_death_genes.txt", sep="\t", quote=FALSE, col.names=FALSE)

death_scores <- fread("figure_reorgnized/Fig6/output_death_genes.txt") %>% 
  column_to_rownames("id") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]], by = "gse_id")
str(death_scores)

death_data <- death_scores %>%
  filter(BPD_type != "mild_BPD") %>%
  select(-BPD_type, -BPD_status)

merged_data <- immune_data %>%
  inner_join(death_data, by = "gse_id") %>%
  select(-gse_id)

pcd_columns <- colnames(death_data)[-1]
immune_columns <- colnames(immune_sig)

correlation_matrix <- cor(merged_data, method = "spearman", use = "pairwise.complete.obs")
cor_subset <- correlation_matrix[immune_columns, pcd_columns]
valid_rows <- rowSums(is.na(cor_subset)) < ncol(cor_subset)
valid_cols <- colSums(is.na(cor_subset)) < nrow(cor_subset)

correlation_subset <- cor_subset[valid_rows, valid_cols]

p_values <- matrix(NA, nrow = sum(valid_rows), ncol = sum(valid_cols))
rownames(p_values) <- rownames(correlation_subset)
colnames(p_values) <- colnames(correlation_subset)

for(i in rownames(correlation_subset)) {
  for(j in colnames(correlation_subset)) {
    test_result <- cor.test(merged_data[[i]], merged_data[[j]], 
                            method = "spearman", exact = FALSE)
    p_values[i,j] <- test_result$p.value
  }
}

significance_matrix <- matrix(
  sapply(p_values, function(p) {
    if(is.na(p)) return("")
    if(p < 0.001) return("***")
    if(p < 0.01) return("**")
    if(p < 0.05) return("*")
    return("")
  }),
  nrow = nrow(p_values)
)

cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.2f", correlation_subset[i, j]), 
            x = x, 
            y = y + height*0.2,  
            gp = gpar(fontsize = 8))
  if(significance_matrix[i, j] != "") {
    grid.text(significance_matrix[i, j], 
              x = x, 
              y = y - height*0.2,  
              gp = gpar(fontsize = 8))
  }
}

column_names <- colnames(correlation_subset)
column_names <- gsub("_", " ", column_names)

##### Fig6C correlation heatmap #####
ht <- Heatmap(correlation_subset,
              name = "Correlation",
              col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
              cell_fun = cell_fun,
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              column_title = "Correlation between Cell Death and Immune Infiltration",
              column_title_gp = gpar(fontsize = 12),
              column_names_rot = 30,
              column_names_side = "top",
              column_names_gp = gpar(fontsize = 10),
              row_names_gp = gpar(fontsize = 10),
              column_dend_side = "bottom", 
              show_row_names = TRUE,
              show_column_names = TRUE,
              column_labels = column_names,
              border = TRUE,
              heatmap_legend_param = list(
                title = "Correlation",
                title_gp = gpar(fontsize = 10),
                labels_gp = gpar(fontsize = 8),
                at = seq(-1, 1, 0.5),
                labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
              ))

pdf("figure_reorgnized/Fig6/Fig6C_cell_death_immune_correlation.pdf", width = 12, height = 8)
draw(ht)
dev.off()

