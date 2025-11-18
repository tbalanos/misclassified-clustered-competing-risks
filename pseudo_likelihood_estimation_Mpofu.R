################################################################################
# Author: Adapted by Theofanis Balanos from Philani Brian Mpofu (2020)
# Purpose: Pseudo-likelihood estimation of misclassification probabilities
#          under clustered competing-risks data with informative cluster size (ICS).
################################################################################

library(sandwich)
library(Hmisc)
library(lmtest)
library(geepack)

# ------------------------------------------------------------------------------
# Function: data_preper (no logical changes)
# ------------------------------------------------------------------------------
data_preper <- function(dat, true_out, surr_out, out_interest, ds_var) {
  dat$r <- dat[, names(dat) %in% c(ds_var)]
  dat$c <- dat[, names(dat) %in% c(true_out)]
  dat$c_obs <- dat[, names(dat) %in% c(surr_out)]
  dat$d1 <- as.numeric(dat$c == out_interest)
  dat$d1_obs <- as.numeric(dat$c_obs == out_interest)
  return(dat)
}

# ------------------------------------------------------------------------------
# Function: prob.compute (only adds id=cluster_id and /Mi weighting)
# ------------------------------------------------------------------------------
prob.compute <- function(dat, formula_pred_val) {
  data_list <- list()
  dat$w1 <- as.numeric((1 - dat$d1_obs) == 1 & dat$r == 1 & dat$c_obs > 0)
  dat$w2 <- as.numeric(dat$d1_obs == 1 & dat$r == 1 & dat$c_obs > 0)
  
  # Convert formula to character for flexible construction
  formula_pred_val <- as.character(formula_pred_val)
  formula_pred_val <- paste(formula_pred_val[1], formula_pred_val[2])
  
  # Explicitly define both true-cause indicators
  dat$d2 <- 1 - dat$d1
  
  pred_val1 <- paste("d1", formula_pred_val)  # P(true=1 | observed=2)
  pred_val2 <- paste("d2", formula_pred_val)  # P(true=2 | observed=1)
  
  # ---- GEE models (with clustering + ICS correction) ----
  model_cause1 <- suppressWarnings(
    geeglm(as.formula(pred_val1), family = binomial(link = "logit"),
           data = dat, id = dat$cluster_id,
           weights = dat$w1 / dat$Mi, corstr = "independence")
  )
  
  model_cause2 <- suppressWarnings(
    geeglm(as.formula(pred_val2), family = binomial(link = "logit"),
           data = dat, id = dat$cluster_id,
           weights = dat$w2 / dat$Mi, corstr = "independence")
  )
  
  
  # Predict probabilities if both models converge
  if (typeof(model_cause1) != "double" && typeof(model_cause2) != "double") {
    dat$p_t1_obs2 <- predict(model_cause1, newdata = dat, type = "response")
    dat$p_t2_obs1 <- predict(model_cause2, newdata = dat, type = "response")
  }
  
  # Augment dataset with predictive probabilities
  if (typeof(model_cause1) != "double" &&
      typeof(model_cause2) != "double" &&
      typeof(dat$p_t1_obs2) != "character" &&
      typeof(dat$p_t2_obs1) != "character") {
    data <- dat
    data$y1 <- (data$r * data$d1) +
      (1 - data$r) * ((dat$d1_obs) * (1 - data$p_t2_obs1) +
                        (1 - dat$d1_obs) * data$p_t1_obs2)
    
    score1 <- estfun(model_cause1)
    score2 <- estfun(model_cause2)
    v1 <- vcov(model_cause1)
    v2 <- vcov(model_cause2)
    
    data_list <- list(data = data, score1 = score1, score2 = score2, v1 = v1, v2 = v2)
  }
  return(data_list)
}

# ------------------------------------------------------------------------------
# Function: coef.mat (adds GEE + ICS correction)
# ------------------------------------------------------------------------------
coef.mat <- function(data, formula_mis) {
  formula_mis <- as.character(formula_mis)
  formula_mis <- paste(formula_mis[1], formula_mis[2])
  
  # Explicitly define both obs-cause indicators
  data$d2_obs <- 1 - data$d1_obs
  
  formula_misclass <- paste("d2_obs", formula_mis)
  
  mod <- suppressWarnings(
    geeglm(as.formula(formula_misclass),
           family = binomial(link = "logit"),
           data = data, id = data$cluster_id,
           weights = data$y1 / data$Mi, corstr = "independence")
  )
  
  mat <- mod$coefficients
  s <- estfun(mod)
  V <- vcov(mod)
  p <- predict(mod, newdata = data, type = "response")
  names(mat) <- names(mod$coefficients)
  
  list(pseudo.coef = mat,
       pseudo.score = s,
       pseudo.V = V,
       pseudo.p = p,
       model_miscl = mod)
}

# ------------------------------------------------------------------------------
# Function: R_create (add /Mi term)
# ------------------------------------------------------------------------------
R_create <- function(dat, formula_pred_val, formula_mis, ps_list) {
  Z <- model.matrix(formula_pred_val, data = dat)
  X <- model.matrix(formula_mis, data = dat)
  w11 <- with(dat, (1 - r) * (1 - d1_obs) * (1 - p_t1_obs2) * p_t1_obs2 / Mi)
  w12 <- with(dat, (1 - r) * d1_obs * (-p_t2_obs1) * (1 - p_t2_obs1) / Mi)
  w2 <- with(dat, d1_obs - ps_list$pseudo.p)
  
  R1 <- (t(Z * w11 * w2) %*% X)
  R2 <- (t(Z * w12 * w2) %*% X)
  list(R1 = R1, R2 = R2)
}

# ------------------------------------------------------------------------------
# Function: sand_compute (adds cluster-sum)
# ------------------------------------------------------------------------------
sand_compute <- function(dat, plug_list, ps_list, R) {
  score_by_cluster <- rowsum(ps_list$pseudo.score, group = dat$cluster_id)
  H <- ps_list$pseudo.V
  v1 <- plug_list$v1
  v2 <- plug_list$v2
  R1 <- R$R1
  R2 <- R$R2
  meat_cluster <- crossprod(score_by_cluster)
  prop_term <- t(R1) %*% v1 %*% R1 + t(R2) %*% v2 %*% R2
  total_meat <- meat_cluster + prop_term
  var1 <- H %*% total_meat %*% t(H)
  var1
}

# ------------------------------------------------------------------------------
# Function: misclass_ps_est (unchanged wrapper)
# ------------------------------------------------------------------------------
misclass_ps_est <- function(dat, true_out, surr_out, out_interest,
                            formula_pred_val, formula_mis, ds_var) {
  dat_set <- data_preper(dat, true_out, surr_out, out_interest, ds_var)
  pred_val_mod <- prob.compute(dat_set, formula_pred_val)
  misMod <- coef.mat(pred_val_mod$data, formula_mis)
  jacob <- R_create(dat = pred_val_mod$data,
                    formula_pred_val = formula_pred_val,
                    formula_mis = formula_mis,
                    ps_list = misMod)
  varEst <- sand_compute(dat = dat_set,
                         plug_list = pred_val_mod,
                         ps_list = misMod,
                         R = jacob)
  parEst <- misMod$pseudo.coef
  test_summary <- coeftest(misMod$model_miscl, varEst)
  
  list(parEst = parEst,
       varEst = varEst,
       test_summary = test_summary,
       model_miscl = misMod$model_miscl)
}
