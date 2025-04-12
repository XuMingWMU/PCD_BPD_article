require(WGCNA)
require(STRINGdb)
require(igraph)
require(ggraph)
require(tidygraph)

allowWGCNAThreads(nThreads = max(1, floor(parallel::detectCores() / 4)))

#### Step 1: preprocess ####
expr_data <- expr_df[["gse32472"]] %>% 
  as.matrix() %>% avereps() %>% t() %>% as.data.frame()

#### Step 2: filter ####
expr_data <- expr_data %>%
  mutate(mean_expression = rowMeans(expr_data)) %>% 
  select(where(~ mean(.x) > 0.1)) %>%  
  select(-mean_expression) %>% 
  avereps() %>% 
  as.data.frame()

#### Step 3: match metadata ####
rownames(meta_df[["gse32472"]]) <- meta_df[["gse32472"]][["gse_id"]]
traitRows <- match(rownames(expr_data), rownames(meta_df[["gse32472"]]))
datTraits <- meta_df[["gse32472"]][traitRows, , drop = FALSE] %>%
  rename(Disease_status = "BPD_status") %>% 
  select(-gse_id) %>%
  mutate(across(everything(), as.numeric))

#### Step 4: top genes ####
top_gene_count = 5000
top_genes <- order(apply(expr_data, 2, mad), decreasing = TRUE)[1:top_gene_count]
expr_data <- expr_data[, top_genes]

#### Step 5: cluster ####
##### Fig3A sample cluster #####
sampleTree <- expr_data %>% dist() %>% hclust(method = "average")
dir.create("figure_reorgnized/Fig3")
pdf(file = "figure_reorgnized/Fig3/Fig3A_sample_clustering.pdf", width = 10, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "", cex.main = 2)
dev.off()

#### Step 6: pick soft threshold ###
##### Fig3B soft threshold #####
cat("Picking soft threshold", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
power_range = 1:20
sft <- pickSoftThreshold(expr_data, powerVector = power_range, verbose = 5)

pdf(file = "figure_reorgnized/Fig3/Fig3B_soft_threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "Scale-Free Topology Model Fit, signed R^2", 
     type = "n", main = "Scale independence")
abline(h = 0.8, col = "red")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = power_range, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", 
     type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = power_range, col = "red")
dev.off()

softPower <- sft$powerEstimate
adjacency <- adjacency(expr_data, power = softPower)
cat("Done picking", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

#### Step 7: module detection ####
minModuleSize = 25
mergeCutHeight = 0.25
net <- blockwiseModules(
  expr_data,
  power = softPower,
  maxBlockSize = 6000,
  TOMType = "unsigned", minModuleSize = minModuleSize,
  reassignThreshold = 0, mergeCutHeight = mergeCutHeight,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = FALSE, verbose = 3
)

##### Fig3D module gene color #####
mergedColors <- labels2colors(net$colors)
pdf("figure_reorgnized/Fig3/Fig3C_module_gene_colors.pdf", width = 12, height = 7)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#### Step 8: module trait relationship ####
MEs <- moduleEigengenes(expr_data, colors = labels2colors(net$colors))$eigengenes
datTraits_numeric <- datTraits %>% dplyr::mutate(across(everything(), as.numeric))
moduleTraitCor <- cor(MEs, datTraits_numeric, use = "p", method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(expr_data))

##### Fig3C module clustering #####
ME_tree <- hclust(dist(t(MEs)), method = "average")
pdf(file = "figure_reorgnized/Fig3/Fig3C_module_eigengene_clustering.pdf", 
    width = 8, height = 5)
par(mar = c(5, 5, 3, 1))
plot(ME_tree,
     main = "Clustering of module eigengenes",
     xlab = "",
     sub = "",
     ylab = "Height",
     hang = 0.03,
     cex = 0.9)
abline(h = mergeCutHeight, col = "red", lty = 2)
dev.off()

##### Fig3E module-trait relationship heatmap #####
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 3), ")")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("figure_reorgnized/Fig3/Fig3E_module_trait_relationships.pdf", width = 4, height = 6)
par(mar = c(4, 6, 3, 2))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits_numeric),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.7, 
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
dev.off()

##### Fig3F-I module-trait correlation #####
mm_gs_dir <- "figure_reorgnized/Fig3/FigF3I_module_membership_vs_significance"
dir.create(mm_gs_dir, showWarnings = FALSE)
all_modules <- unique(mergedColors)
all_modules <- all_modules[all_modules != "grey"]
traits <- c("Disease_status", "BPD_type")
for(trait in traits) {
  trait_dir <- glue("{mm_gs_dir}/{trait}")
  dir.create(trait_dir, showWarnings = FALSE)
  
  GS <- as.numeric(cor(expr_data, datTraits_numeric[[trait]], use = "p"))
  GS_pvalue <- as.numeric(corPvalueStudent(GS, nrow(expr_data)))
  
  module_cors <- data.frame(
    Module = all_modules,
    Correlation = numeric(length(all_modules)),
    P_value = numeric(length(all_modules))
  )
  
  for(i in seq_along(all_modules)) {
    module <- all_modules[i]
    module_genes <- which(mergedColors == module)
    
    module_column <- match(paste0("ME", module), names(MEs))
    MM <- abs(cor(expr_data[, module_genes], MEs[, module_column], use = "p", method = "spearman"))
    MM_pvalue <- corPvalueStudent(MM, nrow(expr_data))
    
    module_GS <- abs(GS[module_genes])
    
    cor_value <- cor(MM, module_GS, use = "p", method = "spearman")
    cor_pvalue <- corPvalueStudent(cor_value, length(MM))
    
    module_cors$Correlation[i] <- cor_value
    module_cors$P_value[i] <- cor_pvalue
    
    pdf(glue("{trait_dir}/MM_GS_{module}_module.pdf"), width = 8, height = 8)
    par(mar = c(5, 5, 4, 2) + 0.1)
    
    plot(MM, module_GS,
         col = module,
         pch = 20,
         main = paste0("Module membership vs. gene significance\n",
                       "cor=", signif(cor_value, 2), ", p=", format(cor_pvalue, scientific = TRUE, digits = 2)),
         xlab = paste(module, "module"),
         ylab = paste0("Gene Significance for ", trait),
         xlim = c(0, 1), 
         ylim = c(0, max(module_GS) * 1.1),
         cex = 0.8)
    
    abline(lm(module_GS ~ MM), col = module, lty = 2)
    abline(h = 0.2, col = "gray", lty = 3)
    abline(v = 0.8, col = "gray", lty = 3)
    
    legend("topleft", 
           legend = c("Genes", "Regression line", "Threshold"),
           col = c(module, module, "gray"),
           pch = c(20, NA, NA),
           lty = c(NA, 2, 3),
           bty = "n")
    
    dev.off()
  }
  
  pdf(glue("{trait_dir}/all_modules_summary.pdf"), width = 12, height = 6)
  par(mar = c(8, 5, 4, 2))
  barplot_pos <- barplot(module_cors$Correlation,
                         names.arg = module_cors$Module,
                         col = module_cors$Module,
                         main = paste0("Module-trait correlations for ", trait),
                         ylab = "MM-GS correlation",
                         las = 2)
  
  abline(h = 0, lty = 2)
  significant <- module_cors$P_value < 0.05
  text(barplot_pos,
       module_cors$Correlation + sign(module_cors$Correlation) * 0.05,
       ifelse(significant, "*", ""),
       cex = 1.5)
  
  dev.off()
  
  summary_stats <- data.frame(
    Module = module_cors$Module,
    Correlation = module_cors$Correlation,
    P_value = module_cors$P_value,
    Gene_count = sapply(module_cors$Module, function(module) {
      sum(mergedColors == module)
    }),
    Mean_GS = sapply(module_cors$Module, function(module) {
      module_genes <- which(mergedColors == module)
      mean(abs(GS[module_genes]))
    })
  )
  
  summary_stats <- summary_stats[order(-abs(summary_stats$Correlation)), ]
  write.csv(summary_stats,
            file = glue("{mm_gs_dir}/module_summary_{trait}.csv"),
            row.names = FALSE)
  
  write.csv(module_cors, 
            file = glue("{mm_gs_dir}/module_correlations_{trait}.csv"), 
            row.names = FALSE)
}

#### Step 9: extract modules ####
geneInfo <- tibble(probes = colnames(expr_data), moduleColor = mergedColors) %>%
  arrange(moduleColor)
write.table(geneInfo, file = "figure_reorgnized/Fig3/Fig3C_all_genes.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#### Step 10: intersection with DEGs ####
gene_wgcna <- list()
wgcna_all_gene <- read.table("figure_reorgnized/Fig3/Fig3C_all_genes.txt", header = T)
gene_venn[["up"]][["wgcna"]] <- wgcna_all_gene %>% 
  filter(moduleColor %in% c("blue", "brown", "green", "midnightblue", 
                            "pink", "purple", "red", "tan", "yellow")) %>% 
  pull(probes)
gene_venn[["up"]][["death"]] <- overlapping_genes[["cell_death"]] %>% 
  unlist() %>% unique()
gene_venn[["down"]][["wgcna"]] <- gene_venn[["up"]][["wgcna"]]
gene_venn[["down"]][["death"]] <- gene_venn[["up"]][["death"]]

gene_up_list_append <- gene_up_list
gene_up_list_append[["WGCNA"]] <- gene_venn[["up"]][["wgcna"]]
gene_up_list_append[["PCD"]] <- death_genes %>% unlist() %>% as.vector() %>% unique()
gene_up_list_append[["gse32472_mild"]] <- NULL
gene_down_list_append <- gene_down_list
gene_down_list_append[["WGCNA"]] <- gene_venn[["down"]][["wgcna"]]
gene_down_list_append[["PCD"]] <- death_genes %>% unlist() %>% as.vector() %>% unique()
gene_down_list_append[["gse32472_mild"]] <- NULL

pdf("figure_reorgnized/Fig3/Fig3J_up_reg_with_wgcna_upset_plot.pdf", height = 4, width = 10)
create_upset_plot(gene_up_list_append, angle = 60)
dev.off()

pdf("figure_reorgnized/Fig3/Fig3K_down_reg_with_wgcna_upset_plot.pdf", height = 4, width = 10)
create_upset_plot(gene_down_list_append, angle = 60)
dev.off()

Reduce(intersect, gene_up_list_append)
lapply(death_genes, function(genes) intersect(Reduce(intersect, gene_up_list_append), genes))

