#### pkgs ####
require(data.table)
require(glue)
require(readxl)
require(limma)
require(DESeq2)
require(clusterProfiler, lib.loc = "/usr/local/lib/R/site-library")
require(org.Hs.eg.db)
require(DOSE)
require(ggvenn)
require(tidyverse)
require(ggsci)
require(ggrepel)
require(ComplexHeatmap)
require(circlize)

#### overview ####
for (name in names(meta_df)) {
  print(name)
  print(str_df[[name]])
  print(table(meta_df[[name]]$BPD_status))
}

## choose gse32472 for initial exploration

#### differential analysis ####

# chip
run_limma_analysis <- function(
    expression_matrix, metadata, prefix, padj = "BH", voom = FALSE,
    grouping_column = "BPD_status", reference_level = "non_BPD") {
  valid_subjects <- colnames(expression_matrix)
  col_data <- metadata %>%
    filter(gse_id %in% valid_subjects) %>%
    arrange(match(gse_id, valid_subjects))
  
  col_data[[grouping_column]] <- as.factor(col_data[[grouping_column]])
  if (!is.null(reference_level)) {
    col_data[[grouping_column]] <- relevel(col_data[[grouping_column]], ref = reference_level)
  }
  
  design <- model.matrix(~ col_data[[grouping_column]])
  
  if (voom) {
    expression_matrix[expression_matrix < 0] <- 0
    expression_matrix <- round(expression_matrix)
    
    v <- voom(expression_matrix, design, plot = FALSE)
    fit <- lmFit(v, design)
  } else {
    fit <- lmFit(expression_matrix, design)
  }
  
  fit <- eBayes(fit)
  results <- topTable(fit, coef = 2, number = Inf, adjust.method = padj)
  
  output_file <- glue("results/DE/{prefix}_limma_{ifelse(voom, 'voom_', '')}differential_expression_results.csv")
  write.csv(results, file = output_file)
  
  return(results)
}

# raw counts
run_deseq2_analysis <- function(
    expression_matrix, metadata, prefix,
    grouping_column = "BPD_status", reference_level = "non_BPD", padj = "BH") {
  valid_subjects <- colnames(expression_matrix)
  col_data <- metadata %>%
    filter(gse_id %in% valid_subjects) %>%
    arrange(match(gse_id, valid_subjects))
  
  col_data[[grouping_column]] <- as.factor(col_data[[grouping_column]])
  if (!is.null(reference_level)) {
    col_data[[grouping_column]] <- relevel(col_data[[grouping_column]], ref = reference_level)
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(expression_matrix),
    colData = col_data,
    design = as.formula(glue("~ {grouping_column}"))
  )
  
  dds <- DESeq(dds)
  
  results <- results(dds, alpha = 0.05, pAdjustMethod = padj)
  
  output_file <- glue("results/DE/{prefix}_deseq2_differential_expression_results.csv")
  write.csv(as.data.frame(results), file = output_file)
  
  return(as.data.frame(results))
}

# run !!!!!!!
diff_res <- list()

for (name in names(str_df)) {
  print(name)
  meta_temp <- meta_df[[name]] %>%
    select(gse_id, BPD_status)
  
  if (str_detect(str_df[[name]], "raw counts")) {
    expr_temp <- expr_df[[name]] %>%
      mutate(across(everything(), as.integer))
    diff_res[[name]] <- run_deseq2_analysis(
      expression_matrix = expr_temp, metadata = meta_temp,
      prefix = glue("DE_{name}")
    )
  } else {
    diff_res[[name]] <- run_limma_analysis(
      expression_matrix = expr_df[[name]], metadata = meta_temp,
      prefix = glue("DE_{name}")
    )
  }
}

diff_res_sig <- list()
for (name in names(str_df)) {
  prefix <- glue("DE_{name}")
  if (str_detect(str_df[[name]], "raw counts")) {
    diff_res_sig[[name]] <- diff_res[[name]] %>%
      filter(padj < 0.05)
    output_file <- glue("results/DE/{prefix}_sig_deseq2_differential_expression_results.csv")
  } else {
    diff_res_sig[[name]] <- diff_res[[name]] %>%
      filter(adj.P.Val < 0.05)
    output_file <- glue("results/DE/{prefix}_sig_limma_differential_expression_results.csv")
  }
  diff_res_sig[[name]] %>%
    write.csv(output_file)
}

# cell death
cell_death_ref1 <- read_xlsx("ref/pcd_rel_genes.xlsx")
cell_death_ref2 <- read_xlsx("ref/10495_2024_1960_MOESM1_ESM/Table S1.xlsx")

names(cell_death_ref1) <- names(cell_death_ref1) %>%
  tolower() %>%
  str_replace_all(" ", "_") %>%
  str_replace_all("-", "_")
names(cell_death_ref2) <- names(cell_death_ref2) %>%
  tolower() %>%
  str_replace_all(" ", "_") %>%
  str_replace_all("-", "_")

names(cell_death_ref1)
names(cell_death_ref2)

cell_death_ref <- bind_cols(cell_death_ref1, cell_death_ref2[, setdiff(names(cell_death_ref2), names(cell_death_ref1))])
write_xlsx(cell_death_ref, "ref/cell_death_finalized.xlsx")

death_genes <- list()
for (death in names(cell_death_ref)) {
  death_genes[[death]] <- cell_death_ref[[death]] %>% na.omit()
}

log_trans <- function(x) {
  x <- log2(x + 1)
}
diff_res[["gse220135_log2"]] <- run_limma_analysis(
  expression_matrix = expr_df[["gse220135"]] %>%
    mutate(across(everything(), log_trans)),
  metadata = meta_df[["gse220135"]], prefix = "DE_gse220135_log2"
)
diff_res_sig[["gse220135_log2"]] <- diff_res[["gse220135_log2"]] %>%
  filter(adj.P.Val < 0.05)

#### overlapping genes ####

# cleaning
expr_df_sel <- list()
meta_df_sel <- list()
str_df_sel <- list()

for(status in c("mild", "moderate", "severe")){
  meta_df_sel[[glue("gse32472_{status}")]] <- meta_df[["gse32472"]] %>%
    filter(BPD_type %in% c(glue("{status}_BPD"), "non_BPD"))
  expr_df_sel[[glue("gse32472_{status}")]] <- expr_df[["gse32472"]] %>%
    select(any_of(meta_df_sel[[glue("gse32472_{status}")]]$gse_id))
  str_df_sel[[glue("gse32472_{status}")]] <- "chip"
}

for(day in c("Day0", "Day14", "Day28")){
  meta_df_sel[[glue("gse220135_{day}")]] <- meta_df[["gse220135"]] %>%
    filter(time == day)
  expr_df_sel[[glue("gse220135_{day}")]] <- expr_df[["gse220135"]] %>%
    select(any_of(meta_df_sel[[glue("gse220135_{day}")]]$gse_id)) %>% 
    mutate(across(everything(), as.integer))
  str_df_sel[[glue("gse220135_{day}")]] <- "raw counts"
}

for(day2 in c("Day28", "Day14")){
  for(day1 in c("Day14", "Day0")){
    if(day2 != day1){
      meta_df_sel[[glue("gse220135_{day2}_vs_{day1}")]] <- meta_df[["gse220135"]] %>%
        filter(time %in% c(day2, day1), 
               BPD_status == "BPD")
      expr_df_sel[[glue("gse220135_{day2}_vs_{day1}")]] <- expr_df[["gse220135"]] %>%
        select(any_of(meta_df_sel[[glue("gse220135_{day2}_vs_{day1}")]]$gse_id)) %>% 
        mutate(across(everything(), as.integer))
      str_df_sel[[glue("gse220135_{day2}_vs_{day1}")]] <- "raw counts"
    }
  }
}

for(name in names(meta_df_sel)){
  if(str_df_sel[[name]] == "raw counts"){
    meta_df_sel[[glue("{name}_log2")]] <- meta_df_sel[[name]]
    expr_df_sel[[glue("{name}_log2")]] <- expr_df_sel[[name]] %>% 
      mutate(across(everything(), log_trans))
    str_df_sel[[glue("{name}_log2")]] <- "log2 raw counts"
  }
}

for(name in names(meta_df_sel)){
  print(name)
  print(table(meta_df_sel[[name]]$BPD_status))
  print(table(meta_df_sel[[name]]$time))
}

# analysis
diff_res_sel <- list()
diff_res_sel_sig <- list()

for (name in names(meta_df_sel)) {
  print(name)
  geo_id <- str_extract(name, "gse\\d+")
  degree <- str_split(name, "_", simplify = TRUE)[2]
  
  if (str_df_sel[[name]] == "chip") {
    meta_temp <- meta_df[[geo_id]] %>%
      select(gse_id, BPD_type)
    grouping_col <- "BPD_type"
    ref_lvl <- "non_BPD"
  } else if (str_detect(name, "vs")) {
    meta_temp <- meta_df[[geo_id]] %>%
      select(gse_id, time) 
    day_vals <- str_extract_all(name, "Day\\d+")[[1]]
    ref_lvl <- day_vals[2]
    grouping_col <- "time"
  } else {
    meta_temp <- meta_df[[geo_id]] %>%
      select(gse_id, BPD_status)
    grouping_col <- "BPD_status"
    ref_lvl <- "non_BPD"
  }
  
  if (str_df_sel[[name]] == "raw counts") {
    diff_res_sel[[name]] <- run_deseq2_analysis(
      expression_matrix = expr_df_sel[[name]],
      metadata = meta_temp,
      grouping_column = grouping_col,
      reference_level = ref_lvl,
      prefix = glue("DE_{name}")
    )
    diff_res_sel_sig[[name]] <- diff_res_sel[[name]] %>%
      filter(padj < 0.05)
    output_file <- glue("results/DE/specific_DE_{name}_sig_deseq2_differential_expression_results.csv")
  } else {
    diff_res_sel[[name]] <- run_limma_analysis(
      expression_matrix = expr_df_sel[[name]],
      metadata = meta_temp,
      grouping_column = grouping_col,
      reference_level = ref_lvl,
      prefix = glue("DE_{name}")
    )
    diff_res_sel_sig[[name]] <- diff_res_sel[[name]] %>%
      filter(adj.P.Val < 0.05)
    output_file <- glue("results/DE/specific_DE_{name}_sig_limma_differential_expression_results.csv")
  }
  
  diff_res_sel_sig[[name]] %>% write.csv(output_file)
}

gene_list <- list(up = list(), 
                  down = list())

gene_list_sel <- list(up = list(), 
                      down = list())

for (name in names(diff_res_sel_sig)) {
  res <- diff_res_sel_sig[[name]]
  
  if ("logFC" %in% names(res)) {
    up_genes <- res %>%
      filter(logFC > 0) %>%
      rownames()
    down_genes <- res %>%
      filter(logFC < 0) %>%
      rownames()
    
    up_genes_sel <- res %>%
      filter(logFC > 0.5) %>%
      rownames()
    down_genes_sel <- res %>%
      filter(logFC < -0.5) %>%
      rownames()
    
  } else if ("log2FoldChange" %in% names(res)) {
    up_genes <- res %>%
      filter(log2FoldChange > 0) %>%
      rownames()
    down_genes <- res %>%
      filter(log2FoldChange < 0) %>%
      rownames()
    
    up_genes_sel <- res %>%
      filter(log2FoldChange > 1) %>%
      rownames()
    down_genes_sel <- res %>%
      filter(log2FoldChange < -1) %>%
      rownames()
  }
  
  gene_list[["up"]][[name]] <- up_genes
  gene_list[["down"]][[name]] <- down_genes
  
  gene_list_sel[["up"]][[name]] <- up_genes_sel
  gene_list_sel[["down"]][[name]] <- down_genes_sel
}

print(names(diff_res_sel_sig))
print(gene_list)

## 32472
gene_venn <- list()

gene_venn[["up"]][["gse32472"]] <- list(
  mild = gene_list_sel[["up"]][["gse32472_mild"]], 
  moderate = gene_list_sel[["up"]][["gse32472_moderate"]], 
  severe = gene_list_sel[["up"]][["gse32472_severe"]]
)

gene_venn[["down"]][["gse32472"]] <- list(
  mild = gene_list_sel[["down"]][["gse32472_mild"]], 
  moderate = gene_list_sel[["down"]][["gse32472_moderate"]], 
  severe = gene_list_sel[["down"]][["gse32472_severe"]]
)

ggvenn(
  gene_venn[["up"]][["gse32472"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 5,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse32472_limma_up_reg_lfc_thre_0.5.pdf", height = 4, width = 4)

ggvenn(
  gene_venn[["down"]][["gse32472"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 5,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse32472_limma_down_reg_lfc_thre_0.5.pdf", height = 4, width = 4)


common_genes <- list()

common_genes[["up"]][["gse32472"]] <- Reduce(intersect, gene_venn[["up"]][["gse32472"]])
common_genes[["down"]][["gse32472"]] <- Reduce(intersect, gene_venn[["down"]][["gse32472"]])

## 220135

for (reg in c("up", "down")) {
  
  gene_venn[[reg]][["gse220135"]] <- list(
    day = list(
      raw = list(),
      log2 = list()
    ),
    status = list(
      raw = list(),
      log2 = list()
    )
  )
  
  for (name in names(gene_list_sel[[reg]])) {
    if (str_detect(name, "220135")) {
      
      if (str_detect(name, "vs")) {
        if (str_detect(name, "log2")) {
          gene_venn[[reg]][["gse220135"]][["day"]][["log2"]][[name]] <- gene_list_sel[[reg]][[name]]
        } else {
          gene_venn[[reg]][["gse220135"]][["day"]][["raw"]][[name]] <- gene_list_sel[[reg]][[name]]
        }
      } else {
        if (str_detect(name, "log2")) {
          gene_venn[[reg]][["gse220135"]][["status"]][["log2"]][[name]] <- gene_list_sel[[reg]][[name]]
        } else {
          gene_venn[[reg]][["gse220135"]][["status"]][["raw"]][[name]] <- gene_list_sel[[reg]][[name]]
        }
      }
    }
  }
}

gene_venn[["up"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day14"]] <- NULL
gene_venn[["down"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day14"]] <- NULL

ggvenn(
  gene_venn[["up"]][["gse220135"]][["status"]][["raw"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse220135_deseq2_up_reg_same_day_lfc_thre_1.5.pdf", height = 4, width = 4)

ggvenn(
  gene_venn[["up"]][["gse220135"]][["day"]][["raw"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse220135_deseq2_up_reg_different_day_lfc_thre_1.5.pdf", height = 4, width = 6)

ggvenn(
  gene_venn[["down"]][["gse220135"]][["status"]][["raw"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse220135_deseq2_down_reg_same_day_lfc_thre_1.pdf", height = 4, width = 4)

ggvenn(
  gene_venn[["down"]][["gse220135"]][["day"]][["raw"]],
  fill_color = c("red", "yellow", "blue"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
)
ggsave("figures/venn/DE_gse220135_deseq2_down_reg_different_day_lfc_thre_1.pdf", height = 4, width = 6)

common_genes[["up"]][["gse220135"]] <- list(
  status = list(
    raw = Reduce(intersect, gene_venn[["up"]][["gse220135"]][["status"]][["raw"]]),
    log2 = Reduce(intersect, gene_venn[["up"]][["gse220135"]][["status"]][["log2"]])
  ),
  day = list(
    raw = Reduce(intersect, gene_venn[["up"]][["gse220135"]][["day"]][["raw"]]),
    log2 = Reduce(intersect, gene_venn[["up"]][["gse220135"]][["day"]][["log2"]])
  )
)

common_genes[["down"]][["gse220135"]] <- list(
  status = list(
    raw = Reduce(intersect, gene_venn[["down"]][["gse220135"]][["status"]][["raw"]]),
    log2 = Reduce(intersect, gene_venn[["down"]][["gse220135"]][["status"]][["log2"]])
  ),
  day = list(
    raw = Reduce(intersect, gene_venn[["down"]][["gse220135"]][["day"]][["raw"]]),
    log2 = Reduce(intersect, gene_venn[["down"]][["gse220135"]][["day"]][["log2"]])
  )
)

# combining 2 GEO datasets
create_upset_plot <- function(test_list, angle = 0) {
  # Get unique gene names
  gene_names <- unique(unlist(test_list))
  
  # Create binary matrix
  gene_set_data <- data.frame(gene = gene_names)
  
  for (set_name in names(test_list)) {
    gene_set_data[[set_name]] <- gene_set_data$gene %in% test_list[[set_name]]
  }
  
  # Convert to matrix format
  m <- make_comb_mat(
    gene_set_data %>% 
      select(-gene) %>% 
      as.matrix()
  )
  
  # Create the upset plot
  p <- UpSet(
    m,
    set_order = names(test_list),
    comb_order = order(comb_size(m), decreasing = TRUE),
    top_annotation = upset_top_annotation(
      m, 
      numbers_rot = angle,
      annotation_name_side = "left",
      show_annotation_name = TRUE,
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 10),
      height = unit(5, "cm")  # 增加柱状图高度
    ),
    right_annotation = upset_right_annotation(
      m, 
      numbers_rot = 0,
      add_numbers = TRUE,
      numbers_gp = gpar(fontsize = 8),
      width = unit(3, "cm")
    ),
    row_gap = unit(0.5, "mm"),  # 减小行间距
    column_gap = unit(0.5, "mm"),  # 减小列间距
    name = "Genes"
  )
  
  # Print the plot
  print(p)
  
  # Return common genes
  common_genes <- Reduce(intersect, test_list)
  return(common_genes)
}

gene_up_list <- list(
  gse32472_mild = gene_venn[["up"]][["gse32472"]][["mild"]], 
  gse32472_moderate = gene_venn[["up"]][["gse32472"]][["moderate"]], 
  gse32472_severe = gene_venn[["up"]][["gse32472"]][["severe"]], 
  gse220135_day28_vs_day0 = unique(
    c(
      gene_venn[["up"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day28_vs_Day0_log2"]],
      gene_venn[["up"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day0"]]
    )
  ), 
  gse220135_day14_vs_day0 = unique(
    c(
      gene_venn[["up"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day14_vs_Day0_log2"]],
      gene_venn[["up"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day14_vs_Day0"]]
    )
  )
)

gene_down_list <- list(
  gse32472_mild = gene_venn[["down"]][["gse32472"]][["mild"]], 
  gse32472_moderate = gene_venn[["down"]][["gse32472"]][["moderate"]], 
  gse32472_severe = gene_venn[["down"]][["gse32472"]][["severe"]], 
  gse220135_day28_vs_day0 = unique(
    c(
      gene_venn[["down"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day28_vs_Day0_log2"]],
      gene_venn[["down"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day0"]]
    )
  ), 
  gse220135_day14_vs_day0 = unique(
    c(
      gene_venn[["down"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day14_vs_Day0_log2"]],
      gene_venn[["down"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day14_vs_Day0"]]
    )
  )
)

common_gene_combined <- list()

common_gene_combined[["up"]] <- create_upset_plot(gene_up_list)
ggsave("figures/venn/#up_reg_upset_plot.pdf", height = 8, width = 12)

common_gene_combined[["down"]] <- create_upset_plot(gene_down_list)
ggsave("figures/venn/#down_reg_upset_plot.pdf", height = 8, width = 12)

common_gene_combined[["up"]] <- Reduce(
  intersect, 
  list(
    # gse32472_mild = gene_venn[["down"]][["gse32472"]][["mild"]], 
    gse32472_moderate = gene_venn[["up"]][["gse32472"]][["moderate"]], 
    gse32472_severe = gene_venn[["up"]][["gse32472"]][["severe"]], 
    gse220135_day28_vs_day0 = unique(
      c(
        gene_venn[["up"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day28_vs_Day0_log2"]],
        gene_venn[["up"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day0"]]
      )
    ), 
    gse220135_day14_vs_day0 = unique(
      c(
        gene_venn[["up"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day14_vs_Day0_log2"]],
        gene_venn[["up"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day14_vs_Day0"]]
      )
    )
  )
)

common_gene_combined[["down"]] <- Reduce(
  intersect, 
  list(
    # gse32472_mild = gene_venn[["down"]][["gse32472"]][["mild"]], 
    gse32472_moderate = gene_venn[["down"]][["gse32472"]][["moderate"]], 
    gse32472_severe = gene_venn[["down"]][["gse32472"]][["severe"]], 
    gse220135_day28_vs_day0 = unique(
      c(
        gene_venn[["down"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day28_vs_Day0_log2"]],
        gene_venn[["down"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day28_vs_Day0"]]
      )
    ), 
    gse220135_day14_vs_day0 = unique(
      c(
        gene_venn[["down"]][["gse220135"]][["day"]][["log2"]][["gse220135_Day14_vs_Day0_log2"]],
        gene_venn[["down"]][["gse220135"]][["day"]][["raw"]][["gse220135_Day14_vs_Day0"]]
      )
    )
  )
)

common_gene_sel <- unique(
  c(common_gene_combined[["up"]], common_gene_combined[["down"]])
)
overlapping_genes[["cell_death"]] <- lapply(death_genes, function(genes) intersect(common_gene_sel, genes))

#### Fig 2A volcanoplot ####

upregulated_by_logFC <- diff_res[["gse32472"]] %>%
  arrange(desc(logFC)) %>%
  head(6) %>% 
  rownames()

downregulated_by_logFC <- diff_res[["gse32472"]] %>%
  arrange(logFC) %>%
  head(6) %>% 
  rownames()

upregulated_by_pvalue <- diff_res[["gse32472"]] %>%
  arrange(P.Value) %>%
  head(6) %>% 
  rownames()

downregulated_by_pvalue <- diff_res[["gse32472"]] %>%
  arrange(P.Value) %>%
  tail(6) %>% 
  rownames()

genes_to_label <- c(upregulated_by_logFC, downregulated_by_logFC, upregulated_by_pvalue, downregulated_by_pvalue) %>%
  unique()

valid_genes <- diff_res[["gse32472"]] %>% 
  rownames_to_column("gene_name") %>% 
  mutate(
    change = case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Up", 
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down", 
      .default = "None"), 
    valid = ifelse(change == "None", "n", "y")
    ) %>% 
  filter(valid == "y")

genes_to_label <- genes_to_label[genes_to_label %in% valid_genes$gene_name]

diff_res[["gse32472"]] %>% 
  rownames_to_column("gene_name") %>% 
  mutate(
    pformat = -log10(adj.P.Val), 
    change = case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Up", 
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down", 
      .default = "None"), 
    change = factor(change, levels = c("Up", "Down", "None")), 
    gene_labels = ifelse(gene_name %in% genes_to_label, 
                         gene_name, "")) %>% 
  ggplot(aes(logFC, pformat))+
  geom_point(aes(color = change))+
  theme_minimal()+
  geom_text_repel(aes(label = gene_labels, color = change), vjust = 0, hjust = 0, size = 3,
                  box.padding = unit(0.35, "lines"), 
                  point.padding = unit(0.6, "lines"), 
                  show.legend = FALSE, 
                  max.overlaps = Inf)+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color="#999999")+
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color="#999999")+
  theme_bw()+
  scale_color_manual(values = c("#FC4E07", "#00AFBB", "#999999"))+
  theme(panel.grid = element_blank(),
        legend.position = c(0.01,0.97),
        legend.justification = c(0,1), 
        legend.key.size = unit(0.2, "cm")
  )+
  scale_x_continuous(limits = c(-1.5, 2), 
                     breaks = seq(-1.5, 2, by = 0.5))

ggsave("figure_reorgnized/Fig2/Fig2A_volcanoplot_32472.pdf", height = 8, width = 10)

#### Fig 2B DEGs in 32472 ####

top_genes_32472 <- diff_res_sig[["gse32472"]] %>%
  mutate(gene = rownames(diff_res_sig[["gse32472"]])) %>%
  arrange(desc(abs(logFC))) %>%
  group_by(sign(logFC)) %>%
  slice_head(n = 100) %>%
  pull(gene)

mat_32472 <- expr_df[["gse32472"]][top_genes_32472, ] %>% as.matrix()
mat_32472 <- t(scale(t(mat_32472)))

anno_list <- list(
  BPD = c("non_BPD" = "#00BFC4", "BPD" = "#F8766D"),
  BPD_type = c(
    "non_BPD" = "#00BFC4", 
    "mild_BPD" = "#7CAE00",
    "moderate_BPD" = "#C77CFF",
    "severe_BPD" = "#F8766D"
  )
)

# 创建注释对象
ha <- HeatmapAnnotation(
  df = list(
    BPD = meta_df[["gse32472"]]$BPD_status,
    BPD_type = meta_df[["gse32472"]]$BPD_type
  ),
  col = anno_list
)

# 创建热图
gse32472_heatmap <- Heatmap(mat_32472,
                            name = "z-score",
                            top_annotation = ha,
                            column_split = meta_df[["gse32472"]]$BPD_status,
                            show_column_names = FALSE,
                            show_row_names = F,
                            cluster_columns = F,
                            cluster_rows = T,
                            row_names_gp = gpar(fontsize = 6),
                            column_title = "GSE32472 Top DEGs",
                            col = colorRamp2(c(-6, 0, 6), c("steelblue", "white", "red"))
)

pdf("figure_reorgnized/Fig2/Fig2B_DEGs_32472.pdf", width = 6, height = 4)
draw(gse32472_heatmap)
dev.off()

#### Fig2C PCD pattern heatmap in 32472 ####

# Get DEGs and filter for PCD genes
top_genes_32472 <- diff_res_sig[["gse32472"]] %>%
  mutate(gene = rownames(diff_res_sig[["gse32472"]])) %>%
  arrange(desc(abs(logFC))) %>%
  group_by(sign(logFC)) %>%
  slice_head(n = 100) %>%
  pull(gene)

pcd_genes <- intersect(top_genes_32472, unique(unlist(death_genes)))

# Create matrix with only PCD genes
mat_32472 <- expr_df[["gse32472"]][pcd_genes, ] %>% as.matrix()
mat_32472 <- t(scale(t(mat_32472)))

# Assign pathways
gene_pathways <- sapply(rownames(mat_32472), function(gene) {
  for(pathway in names(death_genes)) {
    if(gene %in% death_genes[[pathway]]) return(pathway)
  }
  return("other")
})

pathway_colors <- c(
  apoptosis = "#E41A1C",
  pyroptosis = "#377EB8", 
  ferroptosis = "#4DAF4A",
  autophagy = "#984EA3",
  necroptosis = "#FF7F00",
  cuproptosis = "#FFFF33",
  parthanatos = "#F781BF",
  entotic_cell_death = "#999999",
  netotic_cell_death = "#66C2A5",
  lysosome_dependent_cell_death = "#FC8D62",
  alkaliptosis = "#8DA0CB",
  oxeiptosis = "#E78AC3",
  netosis = "#1B9E77", 
  immunogenic_cell_death = "#D95F02", 
  anoikis = "#7570B3", 
  paraptosis = "#E7298A",
  methuosis = "#66A61E", 
  entosis = "#E6AB02", 
  disufidptosis = "#A65628",
  zinc_dependent_cell_death = "#A6D854",
  other = "grey80"
)

anno_list <- list(
  BPD = c("non_BPD" = "#00BFC4", "BPD" = "#F8766D"),
  BPD_type = c(
    "non_BPD" = "#00BFC4", 
    "mild_BPD" = "#7CAE00",
    "moderate_BPD" = "#C77CFF",
    "severe_BPD" = "#F8766D"
  )
)

ha_top <- HeatmapAnnotation(
  df = list(
    BPD = meta_df[["gse32472"]]$BPD_status,
    BPD_type = meta_df[["gse32472"]]$BPD_type
  ),
  col = anno_list,
  annotation_name_rot = 45,
  annotation_name_gp = gpar(fontsize = 8)
)

ha_right <- rowAnnotation(
  PCD_pathway = gene_pathways,  # Changed from Death_Pathway to PCD_pathway
  col = list(PCD_pathway = pathway_colors),  # Update the col list key to match
  annotation_name_side = "top",
  annotation_name_rot = 45,
  annotation_name_gp = gpar(fontsize = 8)
)

gse32472_death_heatmap <- Heatmap(mat_32472,
                                  name = "z-score",
                                  top_annotation = ha_top,
                                  right_annotation = ha_right,
                                  column_split = meta_df[["gse32472"]]$BPD_status,
                                  show_column_names = FALSE,
                                  show_row_names = TRUE,
                                  cluster_columns = FALSE,
                                  cluster_rows = TRUE,
                                  row_names_gp = gpar(fontsize = 7),
                                  column_title_gp = gpar(fontsize = 10),
                                  column_title = "",
                                  col = colorRamp2(c(-6, 0, 6), c("steelblue", "white", "red"))
)

pdf("figure_reorgnized/Fig2/Fig2C_DEGs_PCD_32472.pdf", width = 6, height = 4)
draw(gse32472_death_heatmap)
dev.off()

#### functional enrichment ####

diff_gene <- diff_res_sel_sig[["gse32472_moderate"]] %>% 
  rbind(diff_res_sel_sig[["gse32472_severe"]]) %>% 
  filter(abs(logFC) >= 0.5) %>% 
  rownames() %>% 
  unique()

res_go <- enrichGO(
  gene = diff_gene,
  keyType = "SYMBOL", OrgDb = "org.Hs.eg.db",
  ont = "All", 
  readable = FALSE
) %>% 
  as.data.frame()
  
##### Fig2D GO #####
dat_go <- res_go %>%
  na.omit() %>% 
  mutate(
    GeneRatio = as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("^\\d+/", "", GeneRatio)),
    BgRatio = as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("^\\d+/", "", BgRatio)),
    enrichment_factor = GeneRatio / BgRatio
  ) %>%
  group_by(ONTOLOGY) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()

ggplot(dat_go, aes(x = enrichment_factor, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(2, 8)) +
  scale_x_continuous(limits = c(0, 14)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold"),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(
    x = "Enrichment Factor",
    y = NULL,
    color = "-log10(adj.P)",
    size = "Gene Count",
    title = "GO Enrichment Analysis"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))

ggsave("figure_reorgnized/Fig2/Fig2D_GO_common_top25.pdf", height = 9, width = 8)

##### Fig2E KEGG ####
diff_gene_ent <- bitr(diff_gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
res_kegg <- enrichKEGG(
  gene = diff_gene_ent,
  keyType = "kegg",
  organism = "hsa",
  pvalueCutoff = 0.05
)

dat_kegg <- res_kegg %>% 
  as.data.frame() %>% 
  mutate(
    GeneRatio = as.numeric(sub("/\\d+", "", GeneRatio)) / as.numeric(sub("^\\d+/", "", GeneRatio)),
    BgRatio = as.numeric(sub("/\\d+", "", BgRatio)) / as.numeric(sub("^\\d+/", "", BgRatio)),
    enrichment_factor = GeneRatio / BgRatio
  ) %>%
  slice_max(order_by = p.adjust, n = 20) %>% 
  mutate(Description = factor(Description))

ggplot(dat_kegg, aes(x = enrichment_factor, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "blue", high = "red") +
  scale_size_continuous(range = c(3, 8)) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 14, hjust = 0.5)
  ) +
  labs(
    x = "Enrichment Factor",
    y = NULL,
    color = "-log10(adj.P)",
    size = "Gene Count",
    title = "KEGG Pathway Enrichment Analysis"
  ) +
  scale_x_continuous(position = "bottom", limits = c(0, 7.5)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))

ggsave("figure_reorgnized/Fig2/Fig2E_KEGG.pdf", height = 9, width = 8)





