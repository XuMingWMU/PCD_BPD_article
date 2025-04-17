require(doParallel)
require(caret)
require(randomForest)
require(gbm)
require(glmnet)
require(xgboost)
require(rpart)
require(VennDiagram)
require(pheatmap)
require(neuralnet)
require(NeuralNetTools)
require(e1071)
require(kernlab)
require(nnet)
require(Mime1)
require(patchwork)
require(UpSetR)
require(ComplexHeatmap)
require(svglite)
require(pROC)
require(RColorBrewer)
require(ggpubr)
require(rstatix)
require(enrichplot)
require(Rcpp, lib.loc = "/usr/local/lib/R/site-require")
require(tidyverse)
set.seed(123)

#### split data ####
gene_interest <- Reduce(intersect, gene_up_list_append) %>% unique()

ml_test_pre <- expr_df[["gse32472"]] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]], by = "gse_id") %>% 
  mutate(Var = ifelse(BPD_status == "BPD", "Y", "N"), 
         Var = factor(Var, levels = c("N", "Y"))) %>% 
  select(!c(BPD_status, BPD_type)) %>% 
  rename(ID = "gse_id") %>% 
  select(ID, Var, everything())


prop_table <- prop.table(table(ml_test_pre$Var))
idx_Y <- which(ml_test_pre$Var == "Y")
idx_N <- which(ml_test_pre$Var == "N")

train_idx_Y <- sample(idx_Y, size = floor(length(idx_Y) * 0.4))
train_idx_N <- sample(idx_N, size = floor(length(idx_N) * 0.4))

train_idx <- c(train_idx_Y, train_idx_N)
ml_test_train <- ml_test_pre[train_idx, ]
ml_test_valid <- ml_test_pre[-train_idx, ]

#### train test ####
list_train_vali_Data <- list(
  training = ml_test_train, 
  validation = ml_test_valid
)

ml_test <- ML.Dev.Pred.Category.Sig(ml_test_train, 
                                    list_train_vali_Data, 
                                    candidate_genes = gene_interest,
                                    methods = c('nb', 'svmRadialWeights', 'rf', 'kknn',
                                                'adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 56)

model_names <- c(
  'nb' = 'Naive Bayes',
  'svmRadialWeights' = 'Support Vector Machine',
  'rf' = 'Random Forest',
  'kknn' = 'K-Nearest Neighbors',
  'adaboost' = 'AdaBoost',
  'LogitBoost' = 'Logistic Boosting'
)

auc_vis_category_all(ml_test,
                     dataset = c("training", "validation"),
                     order= c("training", "validation"))

#### Fig4A ML AUC using all candidate genes ####

methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost')

roc_vis_category <- function(object, model_name, dataset, order = NULL, 
                             dataset_col = NULL, anno_position = NULL) {
  require(ggplot2)

  if (is.null(anno_position)) {
    anno_position <- c(0.53, 0.35)
  }

  if (is.null(dataset_col)) {
    dataset_col <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", 
                     "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF", 
                     "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF", 
                     "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF", 
                     "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")
  }

  tmp <- data.frame()
  for (i in 1:length(dataset)) {
    if (model_name != "cancerclass") {
      roc <- cbind(
        object[["roc"]][[dataset[i]]][[model_name]][["TPR"]],
        object[["roc"]][[dataset[i]]][[model_name]][["FPR"]],
        object[["roc"]][[dataset[i]]][[model_name]][["AUC"]]
      )
    } else {
      roc <- cbind(
        object[["roc"]][[dataset[i]]][[model_name]][["sensitivities"]],
        1 - object[["roc"]][[dataset[i]]][[model_name]][["specificities"]],
        object[["roc"]][[dataset[i]]][[model_name]][["auc"]]
      )
    }
    roc <- as.data.frame(roc)
    colnames(roc) <- c("TP", "FP", "AUC")
    roc$dataset <- dataset[i]
    tmp <- rbind(tmp, roc)
  }

  if (is.null(order)) {
    order <- dataset
  }
  tmp$dataset <- factor(tmp$dataset, levels = order)

  p1 <- ggplot(tmp, aes(FP, TP, color = dataset)) +
    geom_line(linewidth = 0.5, alpha = 1) + 
    scale_color_manual(values = dataset_col, name = "Cohort") +
    labs(
      title = paste("AUC predicted by", model_name),
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), 
      plot.title = element_text(hjust = 0.5),
      legend.position = "",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "grey", linewidth = 0.5)

  for (i in 1:length(dataset)) {
    ano_y <- anno_position[2] - (i - 1) * 0.05
    p1 <- p1 +
      annotate(
        "segment",
        x = anno_position[1], xend = anno_position[1] + 0.05,
        y = ano_y, yend = ano_y,
        colour = dataset_col[i]
      ) +
      annotate(
        "text",
        x = anno_position[1] + 0.07, y = ano_y, hjust = 0,
        label = paste(
          "AUC in", order[i], ":", 
          sprintf("%.3f", unique(tmp[tmp$dataset == order[i], "AUC"]))
        )
      )
  }
  
  print(p1)
}

plot_list <- list()

for (method in methods) {
  plot <- roc_vis_category(ml_test,
                           model_name = method,
                           dataset = c("training", "validation"),
                           order = c("training", "validation"),
                           anno_position = c(0.4, 0.25))
  
  plot <- plot + 
    ggtitle(model_names[method])+
    theme(strip.text = element_text(size = 10), 
          strip.background = element_rect(fill = "grey80", colour = NA))
  
  plot_list[[method]] <- plot
}

combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "AUC predicted by machine learning models",
    subtitle = "Modeling with PCD-related DEGs",
    theme = theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
  )

combined_plot

dir.create("figure_reorgnized/Fig4")
ggsave("figure_reorgnized/Fig4/Fig4A_ML_AUC_all_candidate_genes.pdf", combined_plot, height = 9, width = 12)

#### ML variable selection process ####

n_cores <- 20
registerDoParallel(cores = n_cores)
cat("Using", n_cores, "cores for parallel processing\n")

##### pre-screen #####
expr_data <- expr_df[["gse32472"]]
meta_data <- meta_df[["gse32472"]]

design <- model.matrix(~group)
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)
top_genes <- topTable(fit, number=5000, adjust.method="BH")

selected_genes <- unique(c(rownames(top_genes),
                           common_genes %>% unlist() %>% as.vector()))

expr_matrix <- t(expr_data) %>% 
  as.data.frame() %>% 
  select(any_of(ml_test[["sig.gene"]])) %>% 
  as.matrix()

expr_matrix <- expr_matrix[meta_data$gse_id,]
group <- meta_data$BPD_status

##### GBM #####

fitControl <- trainControl(method = "repeatedcv", 
                           number = 5, 
                           repeats = 4)

gbm_fit <- train(x = expr_matrix,
                 y = group,
                 method = "gbm",
                 trControl = fitControl,
                 verbose = FALSE)

gbm_importance <- varImp(gbm_fit)
gbm_genes <- rownames(gbm_importance$importance)[gbm_importance$importance$Overall > quantile(gbm_importance$importance$Overall, 0.5)]

##### decision tree #####
dt_model <- rpart(group ~ ., 
                  data = data.frame(group = group, expr_matrix))

dt_importance <- varImp(dt_model)

dt_genes <- rownames(dt_importance)[dt_importance$Overall > quantile(dt_importance$Overall, 0.5)]

##### lasso #####
x_matrix <- as.matrix(expr_matrix)
cv_fit <- cv.glmnet(x_matrix, 
                    group, 
                    family = "binomial",
                    alpha = 1,
                    type.measure = 'deviance',
                    nfolds = 10,
                    parallel = TRUE)

lasso_coef <- coef(cv_fit, s = cv_fit$lambda.min)
lasso_genes <- rownames(lasso_coef)[which(lasso_coef != 0)][-1]

##### random forest #####
rf_model <- randomForest(x = expr_matrix,
                         y = group,
                         ntree = 500,
                         parallel = TRUE)

rf_importance <- importance(rf_model) %>% 
  as.data.frame()
rf_genes <- rf_importance %>% 
  mutate(MeanDecreaseGini = as.numeric(MeanDecreaseGini)) %>% 
  filter(MeanDecreaseGini > quantile(rf_importance$MeanDecreaseGini, 0.5)) %>% 
  rownames()

##### XGBoost #####
dtrain <- xgb.DMatrix(data = as.matrix(expr_matrix), 
                      label = as.numeric(group) - 1)

params <- list(
  objective = "binary:logistic",
  max_depth = 6,
  eta = 0.3,
  nthread = 20, 
  eval_metric = "error"
)

xgb_direct <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  verbose = 1
)

importance_matrix <- xgb.importance(
  feature_names = colnames(expr_matrix),
  model = xgb_direct
)

xgb_genes <- importance_matrix$Feature[importance_matrix$Gain > quantile(importance_matrix$Gain, 0.5)]

##### SVM-RFE #####
ctrl <- rfeControl(functions = caretFuncs,
                   method = "cv",
                   number = 10,
                   returnResamp = "all",
                   allowParallel = TRUE)
Profile <- rfe(x = expr_matrix,
               y = group, 
               sizes = c(2, 5, 10, 20, 30),
               rfeControl = ctrl,
               method = "svmRadial",
               metric = "Accuracy", 
               preProcess = c("center", "scale"))

stopImplicitCluster()

svm_genes <- Profile$optVariables

# pdf(file="figure_reorgnized/Fig4/Fig4_append_SVM-RFE.pdf", width=6, height=5.5)
# par(las=1)
# x = Profile$results$Variables
# y = Profile$results$RMSE
# plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
# lines(x, y, col="darkgreen")
# wmin=which.min(y)
# wmin.x=x[wmin]
# wmin.y=y[wmin]
# points(wmin.x, wmin.y, col="blue", pch=16)
# text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
# dev.off()

##### Fig4B common selection of ML models #####
common_genes_ai <- Reduce(intersect, list(
  gbm_genes,
  rf_genes,
  dt_genes,
  lasso_genes,
  svm_genes, 
  xgb_genes
))

common_genes_ai

count_gene_appearances <- function(test_list) {
  # Get all unique genes
  all_genes <- unique(unlist(test_list))
  
  # Create a counting dataframe
  gene_counts <- data.frame(
    gene = all_genes,
    frequency = sapply(all_genes, function(g) {
      sum(sapply(test_list, function(x) g %in% x))
    })
  )
  
  # Sort by frequency in descending order
  gene_counts <- gene_counts[order(-gene_counts$frequency), ]
  
  # Add presence/absence in each model
  for (model_name in names(test_list)) {
    gene_counts[[model_name]] <- sapply(gene_counts$gene, function(g) g %in% test_list[[model_name]])
  }
  
  # Convert boolean to more readable format
  gene_counts[names(test_list)] <- lapply(gene_counts[names(test_list)], function(x) ifelse(x, "✓", ""))
  
  return(gene_counts)
}

results <- count_gene_appearances(
  list(
    'GBM' = gbm_genes,
    'Random Forest' = rf_genes,
    'Decision Tree' = dt_genes,
    'LASSO' = lasso_genes,
    'SVM-RFE' = svm_genes,
    'XGBoost' = xgb_genes
  )
)

print(results)

plot_data_left <- results %>%
  select(gene, frequency)

p1 <- ggplot(plot_data_left, aes(y = reorder(gene, frequency))) +
  geom_segment(aes(x = 0, xend = frequency, 
                   yend = reorder(gene, frequency)),
               color = "grey70", size = 0.8) +
  geom_point(aes(x = frequency), size = 3, color = "#4169E1") +
  geom_text(aes(x = frequency, label = frequency),
            hjust = -1.2) + 
  scale_x_continuous(limits = c(0, 7), 
                     breaks = seq(0, 6, 1)) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(vjust = -1) 
  ) +
  labs(x = NULL)

plot_data_right <- results %>%
  select(-frequency) %>%
  gather(key = "Model", value = "Selected", -gene) %>%
  mutate(
    Selected = ifelse(Selected == "✓", "✓", ""),
    gene = factor(gene, levels = levels(reorder(plot_data_left$gene, plot_data_left$frequency)))
  )

p2 <- ggplot(plot_data_right, aes(x = Model, y = gene)) +
  geom_tile(fill = "white", color = "grey90") +
  geom_text(aes(label = Selected), 
            size = 4,
            color = "#2E8B57") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.title = element_blank(), 
    axis.text.y = element_blank()
  )

combined_plot <- p1 + p2 +
  plot_layout(widths = c(2, 3))

cairo_pdf("figure_reorgnized/Fig4/Fig4B_common_gene_selection_ML.pdf", width = 6, height = 5)
combined_plot
dev.off()

common_gene_expr <- expr_data[common_genes_ai, ]
write.table(common_gene_expr,
            file = "figure_reorgnized/Fig4/BPD_feature_genes_expression.txt",
            sep = "\t",
            quote = FALSE)

#### remodelling with ML selected genes ####
ml_post_sel <- expr_df[["gse32472"]] %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]], by = "gse_id") %>% 
  mutate(Var = ifelse(BPD_status == "BPD", "Y", "N"), 
         Var = factor(Var, levels = c("N", "Y"))) %>% 
  select(!c(BPD_status, BPD_type)) %>% 
  rename(ID = "gse_id") %>% 
  select(ID, Var, any_of(common_genes_ai))

# Calculate proportion table for post-selection data
prop_table_post <- prop.table(table(ml_post_sel$Var))

# Get indices for stratified sampling
idx_Y_post <- which(ml_post_sel$Var == "Y")
idx_N_post <- which(ml_post_sel$Var == "N")

# Sample 40% for training
train_idx_Y_post <- sample(idx_Y_post, size = floor(length(idx_Y_post) * 0.4))
train_idx_N_post <- sample(idx_N_post, size = floor(length(idx_N_post) * 0.4))

# Combine indices
train_idx_post <- c(train_idx_Y_post, train_idx_N_post)
ml_post_train <- ml_post_sel[train_idx_post, ]
ml_post_valid <- ml_post_sel[-train_idx_post, ]

# Create list of training and validation data
list_train_vali_Data_post <- list(
  training = ml_post_train, 
  validation = ml_post_valid
)

# Run machine learning models
ml_post <- ML.Dev.Pred.Category.Sig(ml_post_train, 
                                    list_train_vali_Data_post, 
                                    candidate_genes = common_genes_ai,
                                    methods = c('nb', 'svmRadialWeights', 'rf', 'kknn',
                                                'adaboost', 'LogitBoost', 'cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 56)

# Visualize overall AUC results
auc_vis_category_all(ml_post,
                     dataset = c("training", "validation"),
                     order = c("training", "validation"))

# Create ROC curves for each method
methods <- c('nb', 'svmRadialWeights', 'rf', 'kknn', 'adaboost', 'LogitBoost')
plot_list_post <- list()

for (method in methods) {
  plot <- roc_vis_category(ml_post,
                           model_name = method,
                           dataset = c("training", "validation"),
                           order = c("training", "validation"),
                           anno_position = c(0.4, 0.25))
  
  plot <- plot + 
    ggtitle(model_names[method]) +
    theme(strip.text = element_text(size = 10), 
          strip.background = element_rect(fill = "grey80", colour = NA))
  
  plot_list_post[[method]] <- plot
}

# Combine all ROC plots
combined_plot_post <- wrap_plots(plot_list_post, ncol = 3) +
  plot_annotation(
    title = "AUC predicted by machine learning models",
    subtitle = "Modeling with 3 ML selected genes",
    theme = theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold"), 
      plot.subtitle = element_text(size = 16, hjust = 0.5)
    )
  )

combined_plot_post

##### Fig4C gene of interest AUC #####
ggsave("figure_reorgnized/Fig4/Fig4C_AUC_post_selection.pdf", combined_plot_post, height = 9, width = 12)

#### Single gene enrichment ####
genes_of_interest <- c("THBS1", "ALOX5", "ACSL1")
meta <- meta_df[["gse32472"]]
expr <- expr_df[["gse32472"]]

group <- factor(meta$BPD_status, levels = c("non_BPD", "BPD"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr, design)

contrast_matrix <- makeContrasts(BPD_vs_non_BPD = BPD - non_BPD, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

top_table <- topTable(fit2, adjust.method = "fdr", number = Inf)
ranked_list <- sort(setNames(top_table$logFC, rownames(top_table)), decreasing = TRUE)

kegg_gmt <- "gsea/c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt"
go_gmt <- "gsea/c5.go.v2024.1.Hs.symbols.gmt"

gsea_kegg <- GSEA(
  geneList = ranked_list,
  TERM2GENE = read.gmt(kegg_gmt),
  pvalueCutoff = 0.05,
  verbose = TRUE
)
gsea_go <- GSEA(
  geneList = ranked_list,
  TERM2GENE = read.gmt(go_gmt),
  pvalueCutoff = 0.05,
  verbose = TRUE
)

plot_gsea_multiple_paths <- function(gsea_result, gene, title_prefix, n_top_paths = 10) {
  gene_set_ids <- c()
  for (i in 1:nrow(gsea_result@result)) {
    leading_edge <- unlist(strsplit(gsea_result@result$core_enrichment[i], "/"))
    if (gene %in% leading_edge) {
      gene_set_ids <- c(gene_set_ids, i)
    }
  }
  
  if (length(gene_set_ids) > n_top_paths) {
    gene_set_ids <- gene_set_ids[1:n_top_paths]
  }

  if (length(gene_set_ids) > 0) {
    gsea_plot <- gseaplot2(gsea_result, geneSetID = gene_set_ids, title = paste(title_prefix, gene), base_size = 12)
    print(gsea_plot)
  } else {
    message(paste("Gene", gene, "not found in the leading edge of any enriched gene set"))
  }
}

go_results <- gsea_go
kegg_results <- gsea_kegg

##### FigD ALOX5 #####
plot_gsea_multiple_paths(go_results, "ALOX5", "GO: Running Enrichment Score for")
ggsave("figure_reorgnized/Fig4/Fig4D_alox5_go.pdf", height = 9, width = 12)
# plot_gsea_multiple_paths(kegg_results, "ALOX5", "KEGG: Running Enrichment Score for")
# ggsave("figure_reorgnized/Fig4/Fig4E_alox5_kegg.pdf", height = 9, width = 12)

##### Fig4E-F THBS1 #####
plot_gsea_multiple_paths(go_results, "THBS1", "GO: Running Enrichment Score for")
ggsave("figure_reorgnized/Fig4/Fig4E_thbs1_go.pdf", height = 9, width = 12)
plot_gsea_multiple_paths(kegg_results, "THBS1", "KEGG: Running Enrichment Score for")
ggsave("figure_reorgnized/Fig4/Fig4F_thbs1_kegg.pdf", height = 9, width = 12)

##### FigG-H ACSL1 #####
plot_gsea_multiple_paths(go_results, "ACSL1", "GO: Running Enrichment Score for")
ggsave("figure_reorgnized/Fig4/Fig4G_acsl1_go.pdf", height = 9, width = 12)
plot_gsea_multiple_paths(kegg_results, "ACSL1", "KEGG: Running Enrichment Score for")
ggsave("figure_reorgnized/Fig4/Fig4H_acsl1_kegg.pdf", height = 9, width = 12)


