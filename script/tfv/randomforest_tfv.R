
library(tidyverse)
library(here)
library(mrggsave)
library(Tjazi)
library(randomForest)
library(ggplot2)
library(knitr)
library(patchwork)
library(fastshap)
library(ggbeeswarm)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/tfv/randomforest")
if(!dir.exists(figDir)) dir.create(figDir)

options(mrggsave.dev="pdf", 
        mrggsave.dir=figDir, 
        mrg.script="randomforest_tfv.R")

theme_set(theme_bw())

set.seed(77)

# Load data ---------------------------------------------------------------

df <- read.csv(here(dataDir, "data.csv"), header = TRUE, 
               stringsAsFactors = FALSE) %>% # Prevent automatic conversion of character as factor
  mutate(across(contains("phylum"), ~.x/100)) %>% 
  mutate(across(contains("genus"), ~.x/100))

# Select top 20 genus based on abundance
top_genus <- df %>% dplyr::select(contains("genus")) %>% 
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = everything(), names_to = "genus", values_to = "mean") %>%
  arrange(desc(mean)) %>% 
  slice_head(n=20) %>% 
  pull(genus)

# CLR transformation phylum
xx <- df %>% dplyr::select(contains("phylum"))
xx %>% reframe(total=rowSums(across(everything()))) # Make sure sum up to 100 (compositional data)
xx <- clr_c(counts=xx, samples_are = "rows")
xx %>% reframe(total=rowSums(across(everything()))) # Make sure sum up to 0

# CLR transformation genus
yy <- df %>% dplyr::select(contains("genus"))
yy %>% reframe(total=rowSums(across(everything())))
yy <- clr_c(counts=yy, samples_are = "rows")
yy %>% reframe(total=rowSums(across(everything())))

# Only keep the CLR transformed phylum and genus with top 20 abundance
df <- df %>% dplyr::select(!contains("phylum")) %>%
  dplyr::select(!contains("genus")) %>%
  bind_cols(xx) %>%
  bind_cols(yy[, top_genus])

# Functions ---------------------------------------------------------------

# Prepare data
prep_data <- function(dv, idv){
  
  df2 <- df %>% dplyr::select(all_of(dv), contains(idv), 
                              age, weight, tfv_auc_p_d1)
  
  X <- df2 %>% dplyr::select(contains(idv), 
                             age, weight, tfv_auc_p_d1) %>% as.matrix()
  X <- scale(X) # scaling predictors
  # Y <- df2[, dv]
  Y <- scale(df2[, dv])  # standardize outcome
  
  return(list(X = X, Y = Y)) 
}

# 5-fold out-of-fold CV for null model (intercept-only)
null_oof <- function(y, folds = 5, idx = NULL) {
  y <- as.numeric(y)
  n <- length(y)
  
  if (is.null(idx)) {
    idx <- sample(rep(1:folds, length.out = n))
  }
  
  oof_pred <- rep(NA_real_, n)
  
  for (k in 1:folds) {
    tr <- which(idx != k)
    te <- which(idx == k)
    
    mu_tr <- mean(y[tr])  
    oof_pred[te] <- mu_tr 
  }
  
  rmse <- sqrt(mean((y - oof_pred)^2))
  mae  <- mean(abs(y - oof_pred))
  r2   <- 1 - sum((y - oof_pred)^2) / sum((y - mean(y))^2)
  
  metrics <- tibble(RMSE = rmse, MAE = mae, R2 = r2)
  preds   <- tibble(id = seq_len(n), fold = idx, obs = y, pred = oof_pred, resid = y - oof_pred)
  list(metrics = metrics, preds = preds)
}

# 5-fold out-of-fold CV (no tuning)
rf_oof <- function(X, y, ntree = 1000, folds = 5) {
  y <- as.numeric(y)
  n <- nrow(X)
  idx      <- sample(rep(1:folds, length.out = n))
  oof_pred <- rep(NA_real_, n)
  
  for (k in 1:folds) {
    tr <- which(idx != k)
    te <- which(idx == k)
    rf_fit <- randomForest(
      x         = X[tr, , drop = FALSE],
      y         = y[tr],
      ntree     = ntree,
      importance = FALSE
    )
    oof_pred[te] <- predict(rf_fit, newdata = X[te, , drop = FALSE])
  }
  
  rmse <- sqrt(mean((y - oof_pred)^2))
  mae  <- mean(abs(y - oof_pred))
  r2   <- 1 - sum((y - oof_pred)^2) / sum((y - mean(y))^2)
  
  metrics <- tibble(RMSE = rmse, MAE = mae, R2 = r2)
  preds   <- tibble(id = seq_len(n), fold = idx, obs = y, pred = oof_pred, resid = y - oof_pred)
  list(metrics = metrics, preds = preds)
}

# Final model (fit on all data, with importance)
rf_final_mod <- function(X, y, ntree = 1000) {
  y <- as.numeric(y)
  randomForest(
    x          = X,
    y          = y,
    ntree      = ntree,
    importance = TRUE
  )
}

# Plot function
plot_obs_pred <- function(y, yhat, title, stem) {
  p <- tibble(pred = yhat, obs = y) %>%
    ggplot(aes(pred, obs)) +
    geom_point(alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Prediction", y = "Observation", title = title)
  ggsave(file.path(figDir, paste0(stem, ".pdf")), p, width  = 4.8, height = 4.8)
}

shap_beeswarm <- function(model, X, stem = "rf", nsim = 200) {
  X_df <- as.data.frame(X)
  shap_vals <- fastshap::explain(
    object        = model,
    feature_names = colnames(X_df),
    X             = X_df,
    pred_wrapper  = function(object, newdata) predict(object, newdata),
    nsim          = nsim
  )
  shap_mat <- as.matrix(shap_vals)
  df_long <- shap_mat %>%
    as.data.frame() %>%
    mutate(id = row_number()) %>%
    pivot_longer(-id, names_to = "Feature", values_to = "SHAP") %>%
    left_join(
      X_df %>% mutate(id = row_number()) %>%
        pivot_longer(-id, names_to = "Feature", values_to = "Value"),
      by = c("id", "Feature")
    )
  feat_order <- df_long %>%
    group_by(Feature) %>%
    summarise(m = mean(abs(SHAP))) %>%
    arrange(desc(m)) %>%
    pull(Feature)
  df_long$Feature <- factor(df_long$Feature, levels = rev(feat_order))
  p <- ggplot(df_long, aes(SHAP, Feature, color = Value)) +
    geom_quasirandom(alpha = 0.6, size = 1) +
    scale_color_viridis_c() +
    theme_bw() +
    labs(title = paste0("SHAP Beeswarm (", stem, ")"),
         x = "SHAP value", y = NULL)
  ggsave(file.path(figDir, paste0("shap_beeswarm_", stem, ".pdf")), p, width = 7, height = 6)
}

shap_dependence <- function(model, X, stem = "rf", top_k = 6, nsim = 200) {
  X_df <- as.data.frame(X)
  shap_vals <- fastshap::explain(
    object        = model,
    feature_names = colnames(X_df),
    X             = X_df,
    pred_wrapper  = function(object, newdata) predict(object, newdata),
    nsim          = nsim
  )
  shap_mat <- as.matrix(shap_vals)
  colnames(shap_mat) <- colnames(X_df)
  
  mean_abs <- colMeans(abs(shap_mat), na.rm = TRUE)
  top_vars <- names(sort(mean_abs, decreasing = TRUE))[seq_len(min(top_k, ncol(shap_mat)))]
  
  plots <- lapply(top_vars, function(v) {
    df_plot <- tibble(feature_value = X_df[[v]], shap_value = shap_mat[, v])
    ggplot(df_plot, aes(feature_value, shap_value)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE) +
      labs(x = paste0(v, " (scaled)"), y = "SHAP value") +
      theme_bw(base_size = 10)
  })
  p_all <- wrap_plots(plots, ncol = 2)
  ggsave(
    filename = file.path(figDir, paste0("shap_dependence_", stem, ".pdf")),
    plot     = p_all, width    = 7, height   = 7)
}


# Phylum ---------------------------------------------------------
data_phy <- prep_data("tfv_auc_t_d1", "phylum")

# Out-of-fold performance
oof_phy <- rf_oof(data_phy$X, data_phy$Y, ntree = 1000)

knitr::kable(
  oof_phy$metrics %>%
    mutate(across(everything(), ~ round(., 3))),
  caption = "OOF Performance Metrics (Phylum)"
)

# Final model
mod_phy <- rf_final_mod(data_phy$X, data_phy$Y, ntree = 1000)

# Observation vs Prediction plot
plot_obs_pred(oof_phy$preds$obs, oof_phy$preds$pred,
              "Observation vs Prediction (Phylum)", "obs_pred_phylum")

# Variable importance plot
pdf(file.path(figDir, "vip_phylum.pdf"), width = 7, height = 8)
varImpPlot(mod_phy, type  = 1, main  = "Variable Importance (Phylum)")
dev.off()

## SHAP summary
shap_beeswarm(mod_phy, data_phy$X, stem = "phylum")

# Genus ----------------------------------------------------------

data_gen <- prep_data("tfv_auc_t_d1", "genus")

# Out-of-fold performance
oof_gen <- rf_oof(data_gen$X, data_gen$Y, ntree = 1000)

knitr::kable(
  oof_gen$metrics %>%
    mutate(across(everything(), ~ round(., 3))),
  caption = "OOF Performance Metrics (Genus)"
)

# Null model out-of-fold performance
null_gen <- null_oof(
  y    = data_gen$Y,
  folds = 5,
  idx   = oof_gen$preds$fold
)

knitr::kable(
  null_gen$metrics %>%
    mutate(across(everything(), ~ round(., 3))),
  caption = "OOF Performance Metrics (Genus, Null Model)"
)

# Final model
mod_gen <- rf_final_mod(data_gen$X, data_gen$Y, ntree = 1000)

# Observation vs Prediction plot
plot_obs_pred(oof_gen$preds$obs, oof_gen$preds$pred,
              "Observation vs Prediction (Genus)", "obs_pred_genus")

# Variable importance plot
pdf(file.path(figDir, "vip_genus.pdf"), width = 7, height = 8)
varImpPlot(mod_gen, type = 1, main = "Variable Importance (Genus)")
dev.off()

## SHAP summary + dependence
shap_beeswarm(mod_gen, data_gen$X, stem = "genus")
shap_dependence(mod_gen, data_gen$X, stem = "genus", top_k = 6)
