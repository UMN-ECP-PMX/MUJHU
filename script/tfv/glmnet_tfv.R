
library(tidyverse)
library(here)
library(patchwork)
library(mrggsave)
library(ggpubr)
library(ggExtra)
library(glmnet)
library(Tjazi)

dataDir <- here("data/derived")

figDir <- here("deliv/figure/tfv/glmnet")
if(!dir.exists(figDir)) dir.create(figDir)

options(mrggsave.dev="pdf", 
        mrggsave.dir=figDir, 
        mrg.script="glmnet_tfv.R")

theme_set(theme_bw())

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

# Do a few plots (CLR transformed phylum/genus vs original)
p1 <- map(names(xx), function(.x){
  p <- ggplot()+geom_point(aes(df[[.x]], xx[[.x]]))+
    xlab(paste0("Original\n", .x))+
    ylab(paste0("CLR transformed\n", .x))
  return(p)
}) %>% pmplots::pm_grid(ncol=6)
p1
mrggsave(p1, stem="phylum-clr", width=18, height=10)

p2 <- map(names(yy[, top_genus]), function(.x){
  p <- ggplot()+geom_point(aes(df[[.x]], yy[[.x]]))+
    xlab(paste0("Original\n", .x))+
    ylab(paste0("CLR transformed\n", .x))
  return(p)
}) %>% pmplots::pm_grid(ncol=5)
p2
mrggsave(p2, stem="genus-clr", width=18, height=10)

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
  
  return(list(X, Y))
}

# Perform LOOCV
loocv <- function(dat, alpha, lambda.1se=FALSE){
  
  # Set up LOOCV (n-fold CV where n = number of observations)
  n <- nrow(dat[[1]])
  
  # alpha = 0 for Ridge regression
  # alpha = 1 for LASSO regression
  cvfit <- cv.glmnet(dat[[1]], dat[[2]], alpha=alpha, nfolds=n) 
  
  # Plot CV result (optional)
  plot(cvfit)
  
  # Extract optimal lambda (based on LOOCV)
  optimal_lambda <- cvfit$lambda.min
  if(lambda.1se) optimal_lambda <- cvfit$lambda.1se # Use `lambda.1se` instead of `lambda.min`
  
  mean_value <- mean(dat[[2]])
  cv_mse <- cvfit$cvm[which(cvfit$lambda == optimal_lambda)]  # Extract MSE at best lambda
  
  cat("Optimal lambda:", optimal_lambda, "\n")
  cat("Mean:", mean_value, "\n")
  cat("RMSE:", sqrt(cv_mse), "\n")
  cat("RMSE / Mean:", sqrt(cv_mse) / mean_value, "\n")
  
  return(list(optimal_lambda, cvfit))
}

# Evaluate performance
mape <- function(model, dat){
  y_pred <- predict(model, newx = dat[[1]])
  mape <- mean(abs((dat[[2]] - y_pred) / dat[[2]])) * 100
  cat("Mean Absolute Percentage Error (MAPE):", mape, "%\n")
  return(mape)
}

# Extract coefficients
extract_coef <- function(model){
  out <- coef(model) %>% as.matrix() %>% as.data.frame() %>% 
    rownames_to_column(var="Variable")
  return(out)
}

# Fit a Ridge and a LASSO model, aggregating outputs
model_fit <- function(.data, .lambda.1se=FALSE){
  
  # Determine Lambda
  ridge_loocv  <- loocv(dat=.data, alpha=0, lambda.1se=.lambda.1se) # Returns lambda and elbow plot
  lasso_loocv  <- loocv(dat=.data, alpha=1, lambda.1se=.lambda.1se) # Returns lambda and elbow plot
  
  # Fit final model using the optimal lambda
  ridge_model  <- glmnet(.data[[1]], .data[[2]], alpha=0, lambda=ridge_loocv[[1]])
  lasso_model  <- glmnet(.data[[1]], .data[[2]], alpha=1, lambda=lasso_loocv[[1]])
  
  # Model performance
  ridge_mape  <- mape(ridge_model , dat=.data)
  lasso_mape  <- mape(lasso_model , dat=.data)
  
  # Extract coefficients
  ridge_coef  <- extract_coef(ridge_model)
  lasso_coef  <- extract_coef(lasso_model)
  
  # Assemble output
  ridge  <- list(ridge_loocv , ridge_model , ridge_mape , ridge_coef)
  names(ridge) <- c("loocv", "model", "mape", "coef")
  lasso  <- list(lasso_loocv , lasso_model , lasso_mape , lasso_coef)
  names(lasso) <- c("loocv", "model", "mape", "coef")
  
  list <- list(ridge, lasso)
  names(list) <- c("ridge", "lasso")
  return(list)
}

# Plot coefficients
plot_coeff <- function(.x){
  p <- ggplot(filter(.x, Variable != "(Intercept)"), 
              aes(x = reorder(Variable, s0), y = s0)) +
    geom_bar(stat = "identity") +
    coord_flip() + theme_minimal() +
    labs(title = "Standardized Coefficients \nfrom Ridge Regression",
         x = "Predictor",
         y = "Coefficient (Change in Tenofovir Vaginal Tissue AUC SD \n per 1 SD Increase)")
  return(p) 
}

# Phylum tfv_auc_t_d1 -----------------------------------------------------

# Prep data
data <- prep_data("tfv_auc_t_d1", "phylum")

# Fit model
fit1 <- model_fit(data)

# Plot ridge coefficients
p1 <- plot_coeff(.x=fit1$ridge$coef)
p1
mrggsave(p1, stem="ridge_coeff_phylum", width=6, height=7.5)

# Plot lasso coefficients
p2 <- plot_coeff(.x=fit1$lasso$coef)
p2
mrggsave(p2, stem="lasso_coeff_phylum", width=6, height=7.5)

# Genus tfv_auc_t_d1 ------------------------------------------------------

# Prep data
data <- prep_data("tfv_auc_t_d1", "genus")

# Fit model
fit2 <- model_fit(data)

# Plot coefficients
p1 <- plot_coeff(.x=fit2$ridge$coef)
p1
mrggsave(p1, stem="ridge_coeff_genus", width=6, height=6)

# Plot coefficients
p2 <- plot_coeff(.x=fit2$lasso$coef)
p2
mrggsave(p2, stem="lasso_coeff_genus", width=6, height=6)

# Scatter plot of top genera ----------------------------------------------

scat_plot <- function(x){
  var <- sym(as.character(x))
  p <- df %>% ggplot(aes(x=!!var, y=tfv_auc_t_d1))+
    geom_point()+geom_smooth(method = "lm")+
    ylab("TFVdp cervial exposure\n(fmol*hour/mL)")
  if(x=="age"){
      p<-p+xlab(paste0(str_to_sentence(x), " (years)"))
    }else if(x=="weight"){
      p<-p+xlab(paste0(str_to_sentence(x), " (kg)"))
    }else if(grepl("_p", x)){
      p<-p+xlab("TFV plasma exposure\n(ng*hour/mL)")
    }else{
      p<-p+xlab(paste0("CLR transformed\n", x))
    }
  return(p)}

# Pull top 6 genera/variables with most positive coefficients
top_6_ridge_coeff <- arrange(fit2$ridge$coef, desc(s0)) %>% 
  filter(Variable!="(Intercept)") %>% 
  slice_head(n=6) %>% 
  mutate(Variable=fct_reorder(Variable, s0, .desc=TRUE)) %>% 
  arrange(Variable) %>% pull(Variable)

plist1 <- map(top_6_ridge_coeff, scat_plot)
pmplots::pm_grid(plist1, ncol=2)
mrggsave_last(stem="ridge_postive_genus", width=7.5, height=8)

# Pull bottom 6 genera/variable with most negative coefficients
bottom_6_ridge_coeff <- arrange(fit2$ridge$coef, desc(s0)) %>% 
  filter(Variable!="(Intercept)") %>% 
  slice_tail(n=6) %>% 
  mutate(Variable=fct_reorder(Variable, s0, .desc=TRUE)) %>% 
  arrange(Variable) %>% pull(Variable)

plist2 <- map(bottom_6_ridge_coeff, scat_plot)
pmplots::pm_grid(plist2, ncol=2)
mrggsave_last(stem="ridge_negative_genus", width=7.5, height=8)

# Pull lasso screened genera/variable
lasso_coeff <- fit2$lasso$coef %>% 
  filter(s0!=0) %>% filter(Variable!="(Intercept)") %>% 
  mutate(Variable=fct_reorder(Variable, s0, .desc=TRUE)) %>% 
  arrange(Variable) %>% pull(Variable)

plist3 <- map(lasso_coeff, scat_plot)
pmplots::pm_grid(plist3, ncol=2)
mrggsave_last(stem="lasso_genus", width=7, height=7)


