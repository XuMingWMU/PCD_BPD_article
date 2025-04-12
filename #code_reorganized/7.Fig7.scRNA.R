require(parallel)
require(future)
require(future.apply)
require(patchwork)
require(cowplot)
require(Seurat, lib.loc = "/usr/local/lib/R/site-library")
require(ggalluvial)
require(ggpubr)
require(ggsci)
require(glue)
require(paletteer)
require(ggrepel)
require(readxl)
require(writexl)
require(openxlsx)
require(data.table)
require(readr)
require(coin)
require(psych)
require(spam, lib.loc = "/usr/local/lib/R/site-library")
require(tidyverse)
require(ComplexHeatmap)

options(future.globals.maxSize = 60000 * 1024^2)
options(Seurat.object.assay.version = "v5")
options(future.rng.onMisuse = "ignore")
plan(multisession, workers = 10)
set.seed(5201314)

#### data input ####

rough <- fread("scRNA_matrix/GSE151974_raw_umi_matrix_postfilter.csv.gz") %>%
  column_to_rownames("V1")
metaRef <- fread("scRNA_matrix/GSE151974_cell_metadata_postfilter.csv.gz") %>%
  column_to_rownames("V1") %>%
  rownames_to_column("sampleID")
roughObj <- CreateSeuratObject(
  counts = rough, project = "BPD",
  min.cells = 0
)
fullMeta <- roughObj@meta.data %>%
  rownames_to_column("sampleID") %>%
  left_join(metaRef %>% select(!c(orig.ident:nFeature_RNA)), by = "sampleID")
roughObj@meta.data <- roughObj@meta.data %>%
  rownames_to_column("sampleID") %>%
  select(sampleID) %>%
  left_join(fullMeta %>% select(sampleID, nCount_RNA:percent.mito, Oxygen, Sample)) %>%
  column_to_rownames("sampleID")

sam.name <- "figure_reorgnized/Fig7/test"
if (!dir.exists(sam.name)) {
  dir.create(sam.name, recursive = T)
}

roughObj[["percent.mt"]] <- PercentageFeatureSet(roughObj,
                                                 pattern = "^MT-|^mt-"
)

VlnPlot(roughObj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, raster = FALSE
)

gene.freq <- do.call("cbind", tapply(roughObj@meta.data$nFeature_RNA,
                                     roughObj@meta.data$Oxygen,
                                     quantile,
                                     probs = seq(0, 1, 0.05)
))
rna.freq <- do.call("cbind", tapply(roughObj@meta.data$nCount_RNA,
                                    roughObj@meta.data$Oxygen,
                                    quantile,
                                    probs = seq(0, 1, 0.05)
))
mt.freq <- do.call("cbind", tapply(roughObj@meta.data$percent.mt,
                                   roughObj@meta.data$Oxygen,
                                   quantile,
                                   probs = seq(0, 1, 0.05)
))
freq.combine <- as.data.frame(cbind(gene.freq, rna.freq, mt.freq))
colnames(freq.combine) <- c(
  paste(colnames(gene.freq), "Gene", sep = "_"),
  paste(colnames(rna.freq), "RNA", sep = "_"),
  paste(colnames(mt.freq), "MT", sep = "_")
)
write.table(freq.combine, file = paste0(sam.name, "/QC-gene_frequency.txt"), quote = F, sep = "\t")
View(freq.combine)

plot1 <- FeatureScatter(roughObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(roughObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2), legend = "none")

#### preprocess ####

cat("Before filter :", nrow(roughObj@meta.data), "cells\n")
table(roughObj@meta.data$Oxygen)

# 10%
roughObj1 <- subset(roughObj,
                    subset =
                      nFeature_RNA > 500 &
                      nCount_RNA > 1000 &
                      nCount_RNA < 20000 &
                      percent.mt < 10
)
cat("After filter :", nrow(roughObj1@meta.data), "cells\n")
table(roughObj1@meta.data$Oxygen)
# 5%
roughObj2 <- subset(roughObj,
                    subset =
                      nFeature_RNA > 500 &
                      nCount_RNA > 1000 &
                      nCount_RNA < 20000 &
                      percent.mt < 5
)
cat("After filter :", nrow(roughObj2@meta.data), "cells\n")
table(roughObj2@meta.data$Oxygen)

saveRDS(roughObj1, file = "rds_backup/BPD_Mt_thre_10.rds")
saveRDS(roughObj2, file = "rds_backup/BPD_Mt_thre_5.rds")
saveRDS(roughObj, file = "rds_backup/BPD_orig.rds")

bpd_obj <- roughObj1

bpd_obj <- NormalizeData(bpd_obj,
                         normalization.method = "LogNormalize",
                         scale.factor = 10000
)
bpd_obj <- FindVariableFeatures(bpd_obj,
                                selection.method = "vst",
                                nfeatures = 1000
)
bpd_obj <- ScaleData(
  object = bpd_obj,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("Oxygen", "percent.mito")
)
bpd_obj <- RunPCA(
  object = bpd_obj,
  features = VariableFeatures(bpd_obj),
  verbose = F, npcs = 50
)

plan(multisession, workers = 20)

print(Sys.time())
bpd_obj <- JackStraw(bpd_obj, num.replicate = 100, dims = 40)
bpd_obj <- ScoreJackStraw(bpd_obj, dims = 1:40)
print(Sys.time())

pdf(paste0("./", sam.name, "/PCA-JackStrawPlot_40.pdf"), width = 6, height = 5)
JackStrawPlot(object = bpd_obj, dims = 1:40)
dev.off()

pdf(paste0("./", sam.name, "/PCA-ElbowPlot.pdf"), width = 6, height = 5)
ElbowPlot(bpd_obj, ndims = 40)
dev.off()

dim.use <- 1:20

bpd_obj <- FindNeighbors(bpd_obj, dims = dim.use)

res <- 0.25
bpd_obj <- FindClusters(bpd_obj, resolution = res)
bpd_obj <- RunUMAP(bpd_obj,
                   dims = dim.use,
                   do.fast = TRUE
)

all.markers <- FindAllMarkers(bpd_obj,
                              only.pos = TRUE,
                              min.pct = 0.3, logfc.threshold = 0.25
)

write.table(all.markers,
            file = glue("{sam.name}/diff_gene_seurat_clusters_{max(dim.use)}PC_res{res}.txt"),
            sep = "\t", quote = F, row.names = F
)

marker.sig <- all.markers %>%
  mutate(Ratio = round(pct.1 / pct.2, 3)) %>%
  filter(p_val_adj <= 0.05)

n <- 20
top <- marker.sig %>%
  group_by(cluster) %>%
  top_n(n, avg_log2FC)
write.xlsx(top, glue("{sam.name}/top_gene_mes_seurat_clusters_res.xlsx"))

DimPlot(
  object = bpd_obj,
  group.by = "seurat_clusters",
  pt.size = 0.2, reduction = "umap", 
  label = T
)

#### annotation ####
cell_type_ref <- data.frame(
  seurat_clusters = factor(0:20),
  cell_type = c(
    "Col13a1+ fib", "gCap", "AT2", "Myofibroblast", "Alv Macrophage", "CD4+ T cell", 
    "Monocyte", "Neutrophil", "Col14a1+ fib", "Pericytes/SMC", "aCap", "B cells", 
    "Int Macrophage", "AT1", "Lipofibroblast", "Art endothelial cells", 
    "Airway epithelial cells", "Fibroblast progenitor", "Fibroblast", 
    "Vein endothelial cells", "Mast basophil"
  )
)

meta_bkup <- bpd_obj@meta.data

meta_new <- bpd_obj@meta.data %>%
  left_join(cell_type_ref, by = "seurat_clusters")
table(meta_new$cell_type, meta_new$seurat_clusters)

rownames(meta_new) <- rownames(meta_bkup)
bpd_obj@meta.data <- meta_new

#### Fig7A dimplot ####
cluster_labels <- setNames(
  glue("{cell_type_ref$seurat_clusters}: {cell_type_ref$cell_type}"),
  as.character(cell_type_ref$seurat_clusters)
)

cluster_colors <- scales::hue_pal()(length(cluster_labels))

DimPlot(
  object = bpd_obj,
  group.by = "seurat_clusters",
  split.by = "Oxygen",
  pt.size = 0.2, 
  reduction = "umap", 
  label = TRUE
)+ 
  scale_color_manual(
    values = cluster_colors,
    labels = unname(cluster_labels)
  ) + 
  labs(title = "UMAP with Cluster Labels and Cell Type Annotations", color = "Cell Type") + 
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 10), 
    axis.title = element_text(size = 10), 
    axis.text = element_text(size = 8), 
    plot.title = element_text(hjust = 0.5), 
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(
    color = guide_legend(
      ncol = 1, 
      byrow = TRUE,
      override.aes = list(size = 4), 
      key.size = 1.5 
    )
  )

ggsave("figure_reorgnized/Fig7/Fig7A_grouped_UMAP_annotated.pdf", width = 12, height = 7)

#### Fig7B heatmap of PCD pattern ####
death_genes_formatted <- list()
for(name in names(death_genes)) {
  death_genes_formatted[[name]] <- sapply(death_genes[[name]], function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
  death_genes_formatted[[name]] <- unique(death_genes_formatted[[name]]) 
}

gene_sets <- list()
for(name in names(death_genes_formatted)) {
  gene_sets[[name]] <- GeneSet(death_genes_formatted[[name]],
                               geneIdType = SymbolIdentifier(),
                               setName = name)
}

gene_set_collection <- GeneSetCollection(gene_sets)
gmtFile <- "figure_reorgnized/Fig7/test/death_genes_mice.gmt"
dir.create(dirname(gmtFile), recursive = TRUE, showWarnings = FALSE)
toGmt(gene_set_collection, gmtFile)

unique_samples <- unique(bpd_obj@meta.data$Sample)

expr_matrices <- list()

for(sample in unique_samples) {
  cells <- WhichCells(bpd_obj, expression = Sample == sample)
  sample_matrix <- GetAssayData(bpd_obj, slot = "data")[, cells]
  avg_expr <- Matrix::rowMeans(sample_matrix)
  expr_matrices[[sample]] <- avg_expr
}

expr_matrix <- do.call(cbind, expr_matrices)

cell_death_gmt <- getGmt(gmtFile, 
                         collectionType = BroadCollection(category = "c3"),
                         geneIdType = SymbolIdentifier())

scaled_matrix <- t(scale(t(expr_matrix)))
gsva_params <- gsvaParam(
  expr = scaled_matrix,
  geneSets = cell_death_gmt,
  kcdf = "Gaussian",
  minSize = 1,
  maxSize = 500,
  maxDiff = FALSE
)

gsva_results <- gsva(param = gsva_params, verbose = TRUE)

breaks <- seq(-2, 2, length.out = 100)
colors <- colorRampPalette(c("navy", "white", "red"))(length(breaks) - 1)

desired_order <- c("P3_Normoxia", "P3_Hyperoxia", 
                   "P7_Normoxia", "P7_Hyperoxia",
                   "P14_Normoxia", "P14_Hyperoxia")
gsva_results <- gsva_results[, desired_order]
col_info <- data.frame(
  Time = factor(sub("_.*", "", desired_order), levels = c("P3", "P7", "P14")),
  Condition = factor(sub(".*_", "", desired_order), levels = c("Normoxia", "Hyperoxia"))
)
rownames(col_info) <- desired_order
all_colors <- pal_simpsons()(8)  
condition_colors <- c("Normoxia" = all_colors[1], "Hyperoxia" = all_colors[2])
time_colors <- c("P3" = all_colors[3], "P7" = all_colors[4], "P14" = all_colors[5])
column_ha <- HeatmapAnnotation(
  Condition = col_info$Condition,
  Time = col_info$Time,
  col = list(
    Condition = condition_colors,
    Time = time_colors
  )
)

pdf("figure_reorgnized/Fig7/Fig7B_PCD_heatmap_in_diff_time_point.pdf", height = 10, width = 9)
Heatmap(gsva_results,
        name = "Score",
        col = colorRampPalette(c("navy", "white", "red"))(100),
        show_column_names = TRUE,
        cluster_columns = FALSE,
        column_names_rot = 45,
        top_annotation = column_ha,
        column_split = factor(col_info$Time, levels = c("P3", "P7", "P14")),
        column_gap = unit(2, "mm"),
        border = TRUE)
dev.off()

#### Fig7C hub gene dotplot ####
DotPlot(bpd_obj, features = convert_gene_case(gene_interest), 
        group.by = "cell_type",
        # split.by = "Oxygen",
        cols = "RdYlBu")&
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        axis.title = element_blank())
ggsave("figure_reorgnized/Fig7/Fig7C_full_candidate_gene_dotplot.pdf", width = 12, height = 7)

#### PCD boxplot ####
cell_scores <- matrix(nrow = ncol(bpd_obj), ncol = length(gene_sets))
colnames(cell_scores) <- names(gene_sets)

for(i in seq_along(gene_sets)) {
  genes <- gene_sets[[i]]@geneIds
  genes <- genes[genes %in% rownames(bpd_obj)]
  if(length(genes) > 0) {
    cell_scores[,i] <- colMeans(GetAssayData(bpd_obj, slot = "data")[genes,])
  }
}

scores_df <- as.data.frame(cell_scores)
scores_df$CellType <- cells_scores$CellType
scores_df$Sample <- cells_scores$Sample

scores_long <- tidyr::pivot_longer(
  scores_df,
  cols = -c(CellType, Sample),
  names_to = "PCD_Type",
  values_to = "Score"
)

scores_long$TimePoint <- factor(sub("_.*", "", scores_long$Sample), 
                                levels = c("P3", "P7", "P14"))
scores_long$Oxygen <- factor(sub(".*_", "", scores_long$Sample), 
                             levels = c("Normoxia", "Hyperoxia"))

format_pcd_name <- function(x) {
  x <- gsub("_", " ", x)
  words <- str_split(x, " ")[[1]]
  if(length(words) > 3) {
    paste(paste(words[1:2], collapse = " "), 
          paste(words[-(1:2)], collapse = " "), 
          sep = "\n")
  } else {
    x
  }
}

scores_long <- scores_long %>%
  mutate(PCD_Type_formatted = sapply(PCD_Type, format_pcd_name))

pcd_max_y <- scores_long %>%
  group_by(PCD_Type_formatted) %>%
  summarise(max_y = max(Score, na.rm = TRUE) * 1, 
            .groups = "drop")
stat.test <- scores_long %>%
  group_by(TimePoint, PCD_Type_formatted, CellType) %>%
  summarise(
    data_n = n(),
    n1 = sum(Oxygen == "Normoxia"),
    n2 = sum(Oxygen == "Hyperoxia"),
    p = tryCatch({
      if(n1 >= 3 && n2 >= 3) {
        test_result <- wilcox.test(Score ~ Oxygen, exact = FALSE)
        test_result$p.value
      } else {
        NA
      }
    }, error = function(e) NA),
    .groups = "drop"
  ) %>%
  mutate(
    p.signif = case_when(
      is.na(p) | data_n < 6 ~ "",
      p > 0.05 ~ "ns",
      p <= 0.05 & p > 0.01 ~ "*",
      p <= 0.01 & p > 0.001 ~ "**",
      p <= 0.001 & p > 0.0001 ~ "***",
      p <= 0.0001 ~ "****"
    ),
    xmin = as.numeric(factor(CellType)) - 0.2,
    xmax = as.numeric(factor(CellType)) + 0.2
  ) %>%
  left_join(pcd_max_y, by = "PCD_Type_formatted") %>%
  mutate(y.position = max_y)

##### FigS2 PCD statistics boxplot #####
ggplot(scores_long, aes(x = CellType, y = Score, fill = Oxygen)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.size = 0.5) +
  facet_grid(PCD_Type_formatted ~ TimePoint, scales = "free_y") +
  geom_bracket(
    data = stat.test,
    aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.signif),
    inherit.aes = FALSE,
    size = 0.5,
    tip.length = 0.005
  ) +
  scale_fill_manual(values = c("Normoxia" = "#619CFF", "Hyperoxia" = "#F8766D")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "top",
    strip.text.y = element_text(size = 8),
    panel.spacing = unit(0.5, "lines"), 
    axis.title = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    strip.text.x = element_text(size = 14), 
    strip.text.y.right = element_text(size = 14)
  ) +
  labs(
    x = NULL,
    y = "Enrichment score",
    fill = "Condition"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

ggsave("figure_reorgnized/Supplementary_Figure/FigS2_PCD_stat_boxplot_days_specified.pdf", height = 45, width = 25)

##### Fig7D hub gene related PCD boxplot #####
selected_pcds <- c("apoptosis", "ferroptosis", "anoikis")
scores_long_filtered <- scores_long %>% 
  filter(PCD_Type %in% selected_pcds)
stat_test_filtered <- stat.test %>%
  filter(PCD_Type_formatted %in% selected_pcds)

ggplot(scores_long_filtered, aes(x = CellType, y = Score, fill = Oxygen)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7, outlier.size = 0.5) +
  facet_grid(PCD_Type_formatted ~ TimePoint, scales = "free_y") +
  geom_bracket(
    data = stat_test_filtered,
    aes(xmin = xmin, xmax = xmax, y.position = y.position, label = p.signif),
    inherit.aes = FALSE,
    size = 0.5,
    tip.length = 0.005
  ) +
  scale_fill_manual(values = c("Normoxia" = "#619CFF", "Hyperoxia" = "#F8766D")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    legend.position = "top",
    strip.text.y = element_text(size = 8),
    panel.spacing = unit(0.5, "lines"), 
    axis.title = element_text(size = 14), 
    legend.title = element_text(size = 14), 
    legend.text = element_text(size = 12), 
    strip.text.x = element_text(size = 14), 
    strip.text.y.right = element_text(size = 14)
  ) +
  labs(
    x = NULL,
    y = "Enrichment score",
    fill = "Condition"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))

ggsave("figure_reorgnized/Fig7/Fig7D_PCD_hub_gene_stat_boxplot_days_specified.pdf", height = 16, width = 20)

