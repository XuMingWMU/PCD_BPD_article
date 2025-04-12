set.seed(5201314)
require(gtsummary)
require(rms)
require(nomogramFormula)
require(pROC)
require(car)
require(dcurves)
require(tidyverse)
require(readxl)
require(mice)
require(VIM)
require(zoo)
require(Amelia)
require(DescTools)
require(glmnet)
require(writexl)
require(dcurves)
require(glue)
require(forestplot)

#### logistic regression ####
common_genes_ai
lapply(death_genes, function(genes) intersect(common_genes_ai, genes))

expr_dca <- expr_df[["gse32472"]] %>% 
  t() %>% 
  as.data.frame() %>% 
  select(any_of(common_genes_ai)) %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]]) %>% 
  select(-BPD_type, -gse_id)
glimpse(expr_dca)

logi_tb_u <- expr_dca %>% 
  tbl_uvregression(
    method = glm,
    y = BPD_status,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  bold_p(t = 0.05) %>%
  bold_labels()

logi_multi <- expr_dca %>% 
  glm(BPD_status ~ ., family = binomial, data = .)
summary(logi_multi)

logi_tb_m <- tbl_regression(logi_multi, exponentiate = TRUE) %>%
  bold_p(t = .05)

logi_merge <- tbl_merge(
  list(logi_tb_u, logi_tb_m),
  tab_spanner = c("Univariable", "Multivariable")
)

logi_merge

# For univariable results
uni_data <- logi_tb_u$table_body %>%
  select(label, estimate, conf.low, conf.high, p.value)

# For multivariable results
multi_data <- logi_tb_m$table_body %>%
  select(label, estimate, conf.low, conf.high, p.value)

# Combine the data
forest_data <- data.frame(
  label = uni_data$label,
  uni_est = uni_data$estimate,
  uni_lower = uni_data$conf.low,
  uni_upper = uni_data$conf.high,
  uni_p = uni_data$p.value,
  multi_est = multi_data$estimate,
  multi_lower = multi_data$conf.low,
  multi_upper = multi_data$conf.high,
  multi_p = multi_data$p.value
)

##### Fig5E forestplot of logistic regression #####
dir.create("figure_reorgnized/Fig5")
pdf("figure_reorgnized/Fig5/Fig5E_hub_gene_logistic_forestplot.pdf", width = 8, height = 4)
forestplot(
  labeltext = forest_data$label,
  mean = cbind(forest_data$uni_est, forest_data$multi_est),
  lower = cbind(forest_data$uni_lower, forest_data$multi_lower),
  upper = cbind(forest_data$uni_upper, forest_data$multi_upper),
  
  title = "Forest Plot of Univariable and Multivariable Logistic Regression",
  xlab = "Odds Ratio (95% CI)",
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.9),
    xlab = gpar(cex = 1)
  ),
  col = fpColors(
    box = c("#3366CC", "#CC3366"),
    lines = c("#3366CC", "#CC3366"),
    zero = "gray50"
  ),
  legend = c("Univariable", "Multivariable"),
  legend_args = fpLegend(
    pos = list(x = 0.15, y = 0.5),
    gp = gpar(col = "#333333", fill = "#FFFFFF")
  ),
  xticks = c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  zero = 1,
  xlog = TRUE,
  grid = TRUE,
  boxsize = 0.2,
  lineheight = "auto",
  graphwidth = unit(4, "inches"),
  ci.vertices = TRUE,
  ci.vertices.height = 0.1
)
dev.off()

#### single gene logistic ####
develop_logistic_model <- function(train_set, test_set, genes = NULL) {
  if (!is.null(genes)) {
    train_set <- train_set[, c("BPD_status", genes), drop = FALSE]
    test_set <- test_set[, c("BPD_status", genes), drop = FALSE]
  }
  
  model <- glm(BPD_status ~ ., data = train_set, family = binomial())
  
  train_prob <- predict(model, train_set, type = "response")
  train_roc <- roc(train_set$BPD_status, train_prob)
  
  if (ncol(test_set) > 1) {
    test_prob <- predict(model, test_set, type = "response")
    test_roc <- roc(test_set$BPD_status, test_prob)
    
    return(list(
      Method = ifelse(is.null(genes), "All Genes", paste("Single Gene:", genes)),
      Train_ROC = train_roc,
      Test_ROC = test_roc,
      Train_AUC = as.numeric(train_roc$auc),
      Test_AUC = as.numeric(test_roc$auc)
    ))
  } else {
    return(list(
      Method = paste("Single Gene:", genes),
      Train_ROC = train_roc,
      Test_ROC = NULL,
      Train_AUC = as.numeric(train_roc$auc),
      Test_AUC = NA
    ))
  }
}

candidate_genes <- colnames(expr_dca)[!colnames(expr_dca) %in% c("BPD_status")]

single_gene_results <- lapply(candidate_genes, function(gene) {
  develop_logistic_model(train_set, test_set, genes = gene)
})

all_genes_result <- develop_logistic_model(train_set, test_set)

model_results <- c(single_gene_results, list(all_genes_result))

auc_results <- lapply(model_results, function(res) {
  data.frame(Method = res$Method, Train_AUC = res$Train_AUC, Test_AUC = res$Test_AUC)
}) %>% bind_rows()
print(auc_results)

roc_data <- lapply(model_results, function(res) {
  list(Method = res$Method, Train_ROC = res$Train_ROC, Test_ROC = res$Test_ROC)
})

tmp <- data.frame()
for (i in 1:length(roc_data)) {
  method_data <- roc_data[[i]]
  if (!is.null(method_data$Train_ROC)) {
    train_roc <- data.frame(
      TP = method_data$Train_ROC$sensitivities,
      FP = 1 - method_data$Train_ROC$specificities,
      dataset = 'training',
      Method = method_data$Method
    )
    tmp <- rbind(tmp, train_roc)
  }
  if (!is.null(method_data$Test_ROC)) {
    test_roc <- data.frame(
      TP = method_data$Test_ROC$sensitivities,
      FP = 1 - method_data$Test_ROC$specificities,
      dataset = 'validation',
      Method = method_data$Method
    )
    tmp <- rbind(tmp, test_roc)
  }
}

##### Fig5A single gene AUC #####
plot_list <- list()
methods <- unique(tmp$Method)
for (method in methods) {
  method_data <- subset(tmp, Method == method)
  train_auc <- sprintf("%.3f", as.numeric(roc_data[[which(methods == method)]]$Train_ROC$auc))
  test_auc <- sprintf("%.3f", as.numeric(roc_data[[which(methods == method)]]$Test_ROC$auc))
  
  p <- ggplot(method_data, aes(FP, TP, color = dataset)) +
    geom_path(linewidth = 0.5, alpha = 1) +
    scale_color_manual(values = c("#3182BDFF", "#E6550DFF")) +
    labs(
      title = method,
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3), 
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white")
    ) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0.01, 0)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                color = "grey", linewidth = 0.5) +
    annotate("segment", x = 0.4, xend = 0.45, y = 0.25, yend = 0.25, colour = "#3182BDFF", size = 1) +
    annotate("text", x = 0.46, y = 0.25, hjust = 0, label = paste("AUC in training: ", train_auc)) +
    annotate("segment", x = 0.4, xend = 0.45, y = 0.2, yend = 0.2, colour = "#E6550DFF", size = 1) +
    annotate("text", x = 0.46, y = 0.2, hjust = 0, label = paste("AUC in validation: ", test_auc))
  
  plot_list[[method]] <- p
}

combined_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = "AUC predicted by Logistic Regression Models",
    theme = theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold")
    )
  )

print(combined_plot)

ggsave("figure_reorgnized/Fig5/Fig5A_single_gene_AUC.pdf", combined_plot, height = 9, width = 9)

#### nomogram and DCA ####

ddist <- datadist(expr_dca)
options(datadist = 'ddist')

all_genes_model <- glm(BPD_status ~ ., data = expr_dca, family = binomial())
lrm_model <- lrm(BPD_status ~ ., data = expr_dca, x = T, y = T)

nom <- nomogram(lrm_model,
                fun = plogis, 
                fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9),
                lp = FALSE,
                funlabel = "Risk of BPD")

##### Fig5B nomogram #####
pdf("figure_reorgnized/Fig5/Fig5B_nomogram.pdf", height = 5, width = 9)
plot(nom)
dev.off()

##### Fig5D calibration curve #####
cal <- calibrate(lrm_model, method = "boot", B = 100)
pdf("figure_reorgnized/Fig5/Fig5D_calibration_curves.pdf", height = 6, width = 6)
plot(cal,
     xlab = "Predicted Probability",
     ylab = "Observed Probability",
     main = "Calibration Curve", 
     subtitles = T,
     las = 1)
dev.off()

##### Fig5F DCA ##### 

dca_pred <- predict(all_genes_model, type = "response")
dca_data <- cbind(expr_dca, `All genes model` = dca_pred)
results <- dca(BPD_status ~ `All genes model`, data = dca_data)

dca_df <- results$dca
dca_df$variable <- factor(
  dca_df$variable,
  levels = c("All genes model", "all", "none"),
  labels = c("nomogram", "All", "None")
)

ggplot(dca_df, aes(x = threshold, y = net_benefit, color = variable)) +
  geom_path(linewidth = 0.4) +
  scale_color_manual(
    values = c(
      "nomogram" = "#E41A1C",
      "All" = "gray70",
      "None" = "black"
    )
  ) +
  labs(
    x = "Threshold probability",
    y = "Net Benefit"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.margin = unit(c(5, 15, 5, 5), "pt"),
    legend.position = "none"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(-0.05, 0.65),
    breaks = seq(0, 0.6, 0.1),
    expand = c(0.01, 0)
  ) +
  annotate("segment", x = 0.6, xend = 0.68, y = 0.55, yend = 0.55, 
           color = "#E41A1C", linewidth = 0.4) +
  annotate("text", x = 0.7, y = 0.55, label = "Nomogram", 
           color = "#E41A1C", hjust = 0, size = 3.5) +
  annotate("segment", x = 0.6, xend = 0.68, y = 0.50, yend = 0.50, 
           color = "gray70", linewidth = 0.4) +
  annotate("text", x = 0.7, y = 0.50, label = "Treat all", 
           color = "gray70", hjust = 0, size = 3.5) +
  annotate("segment", x = 0.6, xend = 0.68, y = 0.45, yend = 0.45, 
           color = "black", linewidth = 0.4) +
  annotate("text", x = 0.7, y = 0.45, label = "Treat none", 
           color = "black", hjust = 0, size = 3.5)
ggsave("figure_reorgnized/Fig5/Fig5F_dca_curve.pdf", width = 4, height = 4)

##### Fig5C ROC of nomogram #####
roc_data <- list()

train_pred <- predict(all_genes_model, type = "response")
train_roc <- roc(expr_dca$BPD_status, train_pred)

cv_pred <- numeric(nrow(expr_dca))
folds <- createFolds(expr_dca$BPD_status, k = 5)

for(i in seq_along(folds)) {
  train_data <- expr_dca[-folds[[i]], ]
  valid_data <- expr_dca[folds[[i]], ]
  
  cv_model <- glm(BPD_status ~ ., data = train_data, family = binomial())
  
  cv_pred[folds[[i]]] <- predict(cv_model, newdata = valid_data, type = "response")
}

valid_roc <- roc(expr_dca$BPD_status, cv_pred)

roc_data[["Nomogram"]] <- list(
  Method = "Nomogram",
  Train_ROC = train_roc,
  Test_ROC = valid_roc
)

tmp <- data.frame(
  TP = c(train_roc$sensitivities, valid_roc$sensitivities),
  FP = c(1 - train_roc$specificities, 1 - valid_roc$specificities),
  dataset = rep(c("training", "validation"), 
                c(length(train_roc$sensitivities), length(valid_roc$sensitivities))),
  Method = "Nomogram"
)

nomogram_roc_plot <- ggplot(tmp, aes(FP, TP, color = dataset)) +
  geom_path(linewidth = 0.5, alpha = 1) +
  scale_color_manual(values = c("#3182BDFF", "#E6550DFF")) +
  labs(
    title = "Nomogram",
    x = "False Positive Rate (1-Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white")
  ) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
              color = "grey", linewidth = 0.5) +
  annotate("segment", x = 0.4, xend = 0.45, y = 0.25, yend = 0.25, 
           colour = "#3182BDFF", size = 1) +
  annotate("text", x = 0.46, y = 0.25, hjust = 0, 
           label = paste("AUC in training: ", sprintf("%.3f", auc(train_roc)))) +
  annotate("segment", x = 0.4, xend = 0.45, y = 0.2, yend = 0.2, 
           colour = "#E6550DFF", size = 1) +
  annotate("text", x = 0.46, y = 0.2, hjust = 0, 
           label = paste("AUC in validation: ", sprintf("%.3f", auc(valid_roc))))

ggsave("figure_reorgnized/Fig5/Fig5C_nomogram_roc_split.pdf", nomogram_roc_plot, height = 4, width = 4)

#### Fig5G hub gene expression boxplot ####
boxdata <- expr_df[["gse32472"]] %>% 
  t() %>% 
  as.data.frame() %>% 
  select(any_of(common_genes_ai)) %>% 
  rownames_to_column("gse_id") %>% 
  left_join(meta_df[["gse32472"]]) %>% 
  select(-BPD_status, -gse_id)

stat_test <- boxdata %>%
  pivot_longer(-BPD_type, names_to = "genes", values_to = "Expression") %>% 
  group_by(genes) %>% 
  wilcox_test(Expression ~ BPD_type) %>%
  add_significance() %>%
  add_xy_position(x = "genes", group = "BPD_type", dodge = 0.8)

boxdata %>% 
  pivot_longer(-BPD_type, names_to = "genes", values_to = "Expression") %>% 
  ggplot(aes(genes, Expression, color = BPD_type))+
  geom_boxplot(
    aes(color = BPD_type),
    lwd = 1.2,
    fill = "white",
    position = position_dodge(0.8),
    width = 0.6,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = BPD_type, group = BPD_type),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
    size = 0.6,
    alpha = 0.8
  ) +
  scale_color_nejm() +
  stat_pvalue_manual(
    stat_test,
    tip.length = 0.01,
    bracket.length = 0.07, 
    hide.ns = TRUE
  ) +
  labs(
    x = NULL,
    y = "Gene Expression",
    color = "BPD Type"
  ) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top"
  )

ggsave("figure_reorgnized/Fig5/Fig5G_hub_gene_expression_statistics.pdf", width = 6, height = 6)
