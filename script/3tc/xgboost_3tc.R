
library(tidyverse)
library(here)
library(mrggsave)
library(Tjazi)
library(tidymodels)
library(xgboost)
library(matrixStats)
library(SHAPforxgboost)
library(ggplot2)
library(knitr)
library(patchwork)
library(ParBayesianOptimization)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/3tc/xgboost")
if(!dir.exists(figDir)) dir.create(figDir)

options(mrggsave.dev="pdf", 
        mrggsave.dir=figDir, 
        mrg.script="xgboost_3tc.R")

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
                              age, weight, X3tc_auc_p_d1)
  
  X <- df2 %>% dplyr::select(contains(idv), 
                             age, weight, X3tc_auc_p_d1) %>% as.matrix()
  X <- scale(X) # scaling predictors
  # Y <- df2[, dv]
  Y <- scale(df2[, dv])  # standardize outcome
  
  return(list(X = X, Y = Y)) 
}

# Parameter tuning with CV
xgb_tune_cv <- function(X, y, nrounds_max = 1000, early_stopping_rounds = 50, grid_size = 120) {
  dtrain <- xgb.DMatrix(data = X, label = y)
  param_grid <- expand.grid(
    max_depth        = c(2L, 3L, 4L, 5L),
    eta              = c(0.03, 0.05, 0.1, 0.2),
    subsample        = c(0.7, 0.8, 0.9),
    colsample_bytree = c(0.7, 0.8, 0.9),
    min_child_weight = c(1, 3, 5),
    lambda           = c(0, 1, 2),
    alpha            = c(0, 0.1, 0.5),
    KEEP.OUT.ATTRS   = FALSE
  )
  if (nrow(param_grid) > grid_size) {
    param_grid <- param_grid[sample(1:nrow(param_grid), grid_size), , drop = FALSE]
  }
  
  best_rmse <- Inf; best_iter <- NA; best_params <- NULL
  
  for (i in seq_len(nrow(param_grid))) {
    params <- as.list(param_grid[i, ])
    params$booster <- "gbtree"
    params$objective <- "reg:squarederror"
    params$eval_metric <- "rmse"
    
    cv <- xgb.cv(
      params = params,
      data = dtrain,
      nrounds = nrounds_max,
      nfold = 5,
      early_stopping_rounds = early_stopping_rounds,
      verbose = 0
    )
    
    rmse <- min(cv$evaluation_log$test_rmse_mean)
    iter <- which.min(cv$evaluation_log$test_rmse_mean)
    
    if (rmse < best_rmse) {
      best_rmse   <- rmse
      best_iter   <- iter
      best_params <- params
    }
  }
  
  best_df <- as_tibble(best_params) %>%
    mutate(best_iter = best_iter, best_rmse = best_rmse)
  
  list(params = best_params, iter = best_iter, rmse = best_rmse, df = best_df)
}

# 5-fold out-of-fold CV
xgb_oof <- function(X, y, params, nrounds, folds = 5) {
  n <- nrow(X)
  oof_pred <- rep(NA_real_, n)
  idx <- sample(rep(1:folds, length.out = n))
  
  for (k in 1:folds) {
    tr <- which(idx != k)
    te <- which(idx == k)
    dtr <- xgb.DMatrix(X[tr, ], label = y[tr])
    dte <- xgb.DMatrix(X[te, ], label = y[te])
    mdl <- xgb.train(params, dtr, nrounds = nrounds, verbose = 0)
    oof_pred[te] <- predict(mdl, dte)
  }
  
  rmse <- sqrt(mean((y - oof_pred)^2))
  mae  <- mean(abs(y - oof_pred))
  r2   <- 1 - sum((y - oof_pred)^2) / sum((y - mean(y))^2)
  
  metrics <- tibble(RMSE = rmse, MAE = mae, R2 = r2)
  preds   <- tibble(id = seq_len(n), fold = idx, obs = y, pred = oof_pred, resid = y - oof_pred)
  list(metrics = metrics, preds = preds)
}

# Final model fit
xgb_final_mod <- function(X, y, best) {
  dtrain <- xgb.DMatrix(X, label = y)
  xgboost(data = dtrain, params = best$params, nrounds = best$iter, verbose = 0)
}

# Plot functions
plot_obs_pred <- function(y, yhat, title, stem) {
  p <- tibble(pred = yhat, obs = y) %>%
    ggplot(aes(pred, obs)) +
    geom_point(alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    labs(x = "Prediction", y = "Observation", title = title)
  ggsave(file.path(figDir, paste0(stem, ".pdf")), p, width = 4.8, height = 4.8)
}

plot_importance <- function(model, X, top_n = NULL, stem = "vip") {
  imp <- xgb.importance(feature_names = colnames(X), model = model) %>% as_tibble()
  imp <- imp %>% slice_max(Gain, n = min(top_n, nrow(imp)))
  p <- ggplot(imp, aes(x = reorder(Feature, Gain), y = Gain)) +
    geom_col() + coord_flip() +
    labs(x = NULL, y = "Gain", title = "XGBoost Feature Importance")
  ggsave(file.path(figDir, paste0(stem, ".pdf")), p, width = 6, height = 7.5)
  list(imp = imp, plot = p)
}

plot_shap <- function(model, X, top_n = NULL, stem = "shap_summary") {
  Xmat <- as.matrix(X)
  sv   <- SHAPforxgboost::shap.values(xgb_model = model, X_train = Xmat)
  slng <- SHAPforxgboost::shap.prep(shap_contrib = sv$shap_score, X_train = as.data.frame(Xmat))
  p <- SHAPforxgboost::shap.plot.summary(slng)
  ggsave(file.path(figDir, paste0(stem, ".pdf")), plot = p, width = 7, height = 6)
  list(slong = slng, plot = p)
}

shap_dependence <- function(model, X, stem = "genus", top_k = 6) {
  dmat <- xgb.DMatrix(as.matrix(X))
  shap <- predict(model, dmat, predcontrib = TRUE)  # n x (p+1)
  shap <- shap[, -ncol(shap), drop = FALSE]
  colnames(shap) <- colnames(X)
  
  mean_abs <- colMeans(abs(shap), na.rm = TRUE)
  top_vars <- names(sort(mean_abs, decreasing = TRUE))[seq_len(min(top_k, ncol(shap)))]
  
  plist <- lapply(top_vars, function(v) {
    df <- tibble(feature_value = X[, v], shap_value = shap[, v])
    ggplot(df, aes(feature_value, shap_value)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = FALSE, color = "blue") + #method = gam/span, capture less noise
      labs(x = paste0(v, " (scaled)"), y = "SHAP value") +
      theme_bw(base_size = 10)
  })
  
  p_all <- wrap_plots(plist, ncol = 2)
  
  ggsave(file.path(figDir, paste0("shap_depend_", stem, ".pdf")),
         p_all, width = 7, height = 7)
}


# Phylum ------------------------------------------------------------
data_phy        <- prep_data("X3tc_auc_t_d1", "phylum")

# Parameter tuning
best_params_phy <- xgb_tune_cv(data_phy$X, data_phy$Y)
knitr::kable(
  best_params_phy$df %>%
    dplyr::select(max_depth, eta, subsample, colsample_bytree,
                  min_child_weight, lambda, alpha, best_iter, best_rmse) %>%
    mutate(across(where(is.numeric), ~round(., 3))),
  caption = "Best XGBoost Parameters (Phylum)"
)

# Out-of-fold performance
oof_phy         <- xgb_oof(data_phy$X, data_phy$Y, best_params_phy$params, best_params_phy$iter)
knitr::kable(
  oof_phy$metrics %>% mutate(across(everything(), ~round(., 3))),
  caption = "OOF Performance Metrics (Phylum)"
)

# Final model
mod_phy         <- xgb_final_mod(data_phy$X, data_phy$Y, best_params_phy)

# Plots
plot_obs_pred(oof_phy$preds$obs, oof_phy$preds$pred,
              "Observation vs Prediction (Phylum)", "obs_pred_phylum")
vip_phy  <- plot_importance(mod_phy, data_phy$X, stem = "vip_phylum")
shap_phy <- plot_shap(mod_phy, data_phy$X, stem = "shap_phylum")

# Genus ---------------------------------------------------------
data_gen <- prep_data("X3tc_auc_t_d1", "genus")

# Parameter tuning
best_params_gen <- xgb_tune_cv(data_gen$X, data_gen$Y)
knitr::kable(
  best_params_gen$df %>%
    dplyr::select(max_depth, eta, subsample, colsample_bytree,
                  min_child_weight, lambda, alpha, best_iter, best_rmse) %>%
    mutate(across(where(is.numeric), ~round(., 3))),
  caption = "Best XGBoost Parameters (Genus)"
)

# Out-of-fold performance
oof_gen <- xgb_oof(data_gen$X, data_gen$Y,
                   best_params_gen$params, best_params_gen$iter)
knitr::kable(
  oof_gen$metrics %>% mutate(across(everything(), ~round(., 3))),
  caption = "OOF Performance Metrics (Genus)"
)

# Final model
mod_gen <- xgb_final_mod(data_gen$X, data_gen$Y, best_params_gen)

# Plots
plot_obs_pred(oof_gen$preds$obs, oof_gen$preds$pred,
              "Observation vs Prediction (Genus)", "obs_pred_genus")
vip_gen  <- plot_importance(mod_gen, data_gen$X, stem = "vip_genus")
shap_gen <- plot_shap(mod_gen, data_gen$X, stem = "shap_genus")

shap_dependence(mod_gen, data_gen$X, stem = "genus",  top_k = 6)

