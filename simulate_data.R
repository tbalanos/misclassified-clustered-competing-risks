################################################################################
# Author: Theofanis Balanos
# Purpose: Simulate clustered synthetic datasets to illustrate the proposed 
#          marginal semiparametric regression method for misclassified competing 
#          risks data with clustering and informative cluster size (ICS).
#
# Description:
# The data-generation mechanism extends the framework of Mpofu et al. (2020) by 
# introducing a hierarchical (clustered) structure and informative cluster size. 
# Cluster-level frailty terms are used only to induce within-cluster dependence 
# and informative cluster size in the simulated data, while the proposed method 
# itself is a marginal model that does not condition on the frailty.
#
# Two dataset types can be generated:
#   (1) External validation dataset (ds = 1): true causes observed; used to 
#       estimate misclassification model parameters.
#   (2) Main analysis dataset (ds = 0): only misclassified causes observed; 
#       analyzed using the marginal semiparametric regression method that 
#       accounts for clustering and ICS through inverse cluster-size weighting.
#
# Each simulated dataset includes:
#   - Cluster and subject identifiers
#   - Covariates (z1, z2)
#   - Cluster-level frailty term (for simulation only)
#   - True and observed causes of failure
#   - Potential misclassification of causes
#   - Optional double-sampling indicators
#
# The simulated data are used solely for replication and illustration purposes 
# in the supplementary R Markdown file. No real IeDEA data are included.
#
# References:
# Mpofu, P. B., Bakoyannis, G., Yiannoutsos, C. T., Mwangi, A., & Mburu, M. 
#   (2020). Semiparametric regression analysis of misclassified competing risks 
#   data. Statistics in Medicine, 39(23–24), 3199–3216.
#
# Zhou, W., Bakoyannis, G., Zhang, Y., Musick, B., Kimaina, A., & Yiannoutsos, C. T. 
#   (2023). Semiparametric marginal regression for clustered competing risks data 
#   with missing cause of failure. Biostatistics, 24(3), 795–813.
################################################################################

# Load required packages
library(TeachingDemos)
library(sandwich)
library(Hmisc)
library(survival)
library(splines)
library(stabledist)

# Set deterministic seed for reproducibility
seed <- char2seed("DoubleSampling", set = FALSE)
set.seed(seed)

# ------------------------------------------------------------------------------
# Simulation settings
# ------------------------------------------------------------------------------

# Global storage for simulation constants
double_sample <- 0     # proportion for double sampling
truth <- NULL          # vector/list of true misclassification parameters
cause_num <- 0         # scenario indicator
#   1: uni-directional misclassification
#   2: bi-directional (equal)
#   3: bi-directional (unequal)

# ------------------------------------------------------------------------------
# Function: true_value()
# Purpose : Return true parameter values for the chosen misclassification scenario
# ------------------------------------------------------------------------------

true_value <- function(scenario) {
  if (scenario == 1) {
    # uni-directional misclassification
    truth <- c(-1.2, -0.4, 0.5, -0.5)
  } else if (scenario == 2) {
    # bi-directional misclassification with equal probabilities
    truth <- c(-0.4, -0.4, 0.5, -0.5)
  } else if (scenario == 3) {
    # bi-directional misclassification with unequal probabilities
    truth <- list(
      c(-0.4, -0.4, 0.5, -0.5),
      c(-0.3, 0.2, 0.5, 0.5)
    )
  } else {
    stop("Scenario provided is not valid. Choose 1, 2, or 3.")
  }
  return(truth)
}

# ------------------------------------------------------------------------------
# Function: simulate()
# ------------------------------------------------------------------------------
# Purpose:
#   Generate a clustered synthetic dataset under the competing risks setting
#   with cause-of-failure misclassification, optional double-sampling,
#   within-cluster dependence, and informative cluster size (ICS).
#
#   This extends the Mpofu et al. (2020) framework by introducing cluster-level
#   frailty terms to induce correlation and ICS, following Zhou et al. (2023).
#   The method itself remains *marginal*, i.e., inference does not condition on
#   the frailty but adjusts for clustering through weighting in estimation.
#
# Input:
#   n_clusters - number of clusters to generate
#   alpha      - frailty variance parameter (positive stable)
#   ds         - double-sampling proportion (0–1)
#   scen       - misclassification scenario (1, 2, or 3)
#
# Output:
#   A data.frame containing event times, true and observed causes, covariates,
#   cluster and subject identifiers, cluster size (Mi), double-sampling
#   indicators, and labeled columns ready for model fitting.
#
# References:
#   Mpofu, P. B., Bakoyannis, G., Yiannoutsos, C. T., Mwangi, A., & Mburu, M. 
#     (2020). Semiparametric regression analysis of misclassified competing risks data. 
#     Statistics in Medicine, 39(23–24), 3199–3216.
#   Zhou, W., Bakoyannis, G., Zhang, Y., Musick, B., Kimaina, A., & Yiannoutsos, C. T. 
#     (2023). Semiparametric marginal regression for clustered competing risks data 
#     with missing cause of failure. Biostatistics, 24(3), 795–813.
# ------------------------------------------------------------------------------

simulate <- function(n_clusters = 100, alpha = 0.4, ds = 1, scen = 1) {
  if (ds > 1) stop("The double-sampling proportion cannot exceed 1")
  
  ## Update global values
  double_sample <<- ds
  truth <<- true_value(scen)
  cause_num <<- scen
  
  ## Covariates
  z1 <- runif(1, min = 0, max = 1)
  z2 <- rnorm(1)
  
  ## Generate cluster-level frailties (positive stable)
  w1_all <- rstable(n_clusters, alpha = alpha, beta = 1,
                    gamma = cos(pi * alpha / 2)^(1 / alpha),
                    delta = 0, pm = 1)
  w2_all <- rstable(n_clusters, alpha = alpha, beta = 1,
                    gamma = cos(pi * alpha / 2)^(1 / alpha),
                    delta = 0, pm = 1)
  median_w1 <- median(w1_all)
  median_w2 <- median(w2_all)
  
  data_list <- list()
  id <- 1
  
  for (i in 1:n_clusters) {
    wi1 <- w1_all[i]
    wi2 <- w2_all[i]
    
    # Stronger ICS based on medians
    if (wi1 < median_w1 && wi2 < median_w2) {
      Mi <- sample(5:10, 1)   # Small clusters
    } else if (wi1 >= median_w1 && wi2 >= median_w2) {
      Mi <- sample(100:120, 1) # Large clusters
    } else {
      Mi <- sample(30:50, 1)   # Medium clusters
    }
    
    ## Loop over subjects within cluster
    for (j in 1:Mi) {
      z1 <- runif(1, 0, 1)
      z2 <- rnorm(1)
      
      ## Weibull cause-specific hazards (proportional hazards)
      lambda1 <- 0.75; alpha1 <- 2
      b11 <- 0.5; b12 <- 1
      zb1 <- b11 * z1 + b12 * z2
      
      lambda2 <- 1; alpha2 <- 2
      b21 <- -0.5; b22 <- 0.5
      zb2 <- b21 * z1 + b22 * z2
      
      ## Censoring rate
      gamma <- 0.6
      
      ## Scenario-specific parameters
      theta <- truth
      
      ## Generate true event times via inversion method
      U <- runif(1)
      t <- sqrt(-log(1 - U) / (wi1 * (lambda1^alpha1) * exp(zb1) + wi2 * (lambda2^alpha2) * exp(zb2)))
      
      ## Cause 1 probability
      p <- (wi1 * alpha1 * (lambda1^alpha1) * (t^(alpha1 - 1)) * exp(zb1)) /
           (wi1 * alpha1 * (lambda1^alpha1) * (t^(alpha1 - 1)) * exp(zb1) +
           wi2 * alpha2 * (lambda2^alpha2) * (t^(alpha2 - 1)) * exp(zb2))
      
      ## True cause of failure (before censoring)
      c <- (2 - rbinom(1, 1, p))
      
      ## Censoring time
      cens_tm <- rexp(1, rate = gamma)
      
      ## Apply censoring
      c <- c * (t <= cens_tm)
      x <- pmin(t, cens_tm)
      
      # -----------------------------------------------------------------------
      # Misclassification mechanism
      # -----------------------------------------------------------------------
      # Uni-directional misclassification (cause 1 may be misclassified as 2)
      if (scen == 1) {
        pm <- plogis(theta[1] + theta[2]*t + theta[3]*z1 + theta[4]*z2)
        u <- runif(1)
        r <- 1 - (u <= pm) * (c == 1)
        c_obs <- r * c + (1 - r) * (3 - c)
        s <- (runif(1) <= ds) * (c_obs == 2)
        r <- as.numeric((c_obs == 0) + (s == 1) + (c_obs == 1) >= 1)
      } else if (scen == 2) {
        # Equal bidirectional misclassification
        pm1 <- plogis(theta[1] + theta[2]*t + theta[3]*z1 + theta[4]*z2)
        pm2 <- pm1
        u <- runif(1)
        r <- 1 - (u <= pm1)*(c == 1) - (u <= pm2)*(c == 2)
        c_obs <- r * c + (1 - r) * (3 - c)
        s <- (c_obs > 0) * rbinom(1, 1, prob = ds)
        r <- as.numeric((c_obs == 0) + (s == 1) >= 1)
      } else if (scen == 3) {
        # Unequal bidirectional misclassification
        pm1 <- plogis(theta[[1]][1] + theta[[1]][2]*t + theta[[1]][3]*z1 + theta[[1]][4]*z2)
        pm2 <- plogis(theta[[2]][1] + theta[[2]][2]*t + theta[[2]][3]*z1 + theta[[2]][4]*z2)
        u <- runif(1)
        r <- 1 - (u <= pm1)*(c == 1) - (u <= pm2)*(c == 2)
        c_obs <- r * c + (1 - r) * (3 - c)
        s <- (runif(1) <= ds) * (c_obs > 0)
        r <- as.numeric((c_obs == 0) + (s == 1) >= 1)
      }
      
      # -----------------------------------------------------------------------
      # Store subject-level record
      # -----------------------------------------------------------------------
      d1 <- as.numeric(c == 1); d2 <- as.numeric(c == 2)
      d1_obs <- as.numeric(c_obs == 1); d2_obs <- as.numeric(c_obs == 2)
      
      data <- data.frame(d1, d2, t, x, z1, z2, s, r, d1_obs, d2_obs, c_obs, c)
      data$id <- 1:nrow(data)
      # Add cluster-level identifiers
      data$cluster_id <- i
      data$Mi <- Mi
      # Save to list
      data_list[[id]] <- data
      id <- id + 1
    }
  }
  
  # Combine all clusters
  data <- do.call(rbind, data_list)
  data$id <- 1:nrow(data)
  
  # Add descriptive labels
  label(data$c) <- "True Cause of Failure"
  label(data$c_obs) <- "Observed Cause of Failure"
  label(data$d1) <- "True cause 1 indicator"
  label(data$d2) <- "True cause 2 indicator"
  label(data$d1_obs) <- "Observed cause 1 indicator"
  label(data$d2_obs) <- "Observed cause 2 indicator"
  label(data$s) <- "Double-sampling indicator"
  label(data$r) <- "Indicator if true cause observed (1=yes)"
  label(data$t) <- "Failure time"
  label(data$x) <- "Observed time (min of event/censoring)"
  label(data$z1) <- "Covariate z1"
  label(data$z2) <- "Covariate z2"
  label(data$cluster_id) <- "Cluster identifier"
  label(data$Mi) <- "Cluster size"
  
  return(data)
}

# ------------------------------------------------------------------------------
# Function: prob.compute()
# ------------------------------------------------------------------------------
# Purpose:
#   Compute predictive values for the true cause of failure given the observed
#   cause, following the double-sampling pseudo-likelihood framework of
#   Mpofu et al. (2020). These predictive probabilities are estimated using
#   logistic regression models fitted on the double-sampled (validation) subset.
#
# Input:
#   list_obj - a list of simulated datasets created by data_list()
#
# Output:
#   A list of elements, one per dataset, each containing:
#     1. The augmented dataset with predicted probabilities (p_t1_obs2, p_t2_obs1)
#     2. Score contributions from the logistic models
#     3. Covariance matrices of model parameters
#
# Reference:
#   Mpofu, P. B., Bakoyannis, G., Yiannoutsos, C. T., Mwangi, A., & Mburu, M. (2020).
#   *Semiparametric regression analysis of misclassified competing risks data.*
#   Statistics in Medicine, 39(23–24), 3199–3216.
# ------------------------------------------------------------------------------

prob.compute <- function(list_obj) {
  n <- length(list_obj)
  data_list <- list()
  
  for (i in 1:n) {
    dat <- list_obj[[i]]
    
    # Logistic model weights identify double-sampled observations
    dat$w1 <- as.numeric(dat$d2_obs == 1 & dat$r == 1)  # model for true cause 1 | observed 2
    dat$w2 <- as.numeric(dat$d1_obs == 1 & dat$r == 1)  # model for true cause 2 | observed 1
    
    # --------------------------------------------------------------------------
    # Fit logistic models for predictive values
    # --------------------------------------------------------------------------
    model_cause1 <- tryCatch(
      glm(d1 ~ t + z1 + z2, family = binomial(link = "logit"),
          weights = w1, data = dat),
      warning = function(w) return(99)
    )
    
    model_cause2 <- tryCatch(
      glm(d2 ~ t + z1 + z2, family = binomial(link = "logit"),
          weights = w2, data = dat),
      warning = function(w) return(99)
    )
    
    # Skip datasets with quasi-separation or estimation failure
    if (typeof(model_cause1) != "double" && typeof(model_cause2) != "double") {
      dat$p_t1_obs2 <- tryCatch(
        predict(model_cause1, newdata = dat, type = "response"),
        warning = function(w) return("Stop")
      )
      dat$p_t2_obs1 <- tryCatch(
        predict(model_cause2, newdata = dat, type = "response"),
        warning = function(w) return("Stop")
      )
    }
    
    # --------------------------------------------------------------------------
    # Combine model results if both models are valid
    # --------------------------------------------------------------------------
    if (typeof(model_cause1) != "double" &&
        typeof(model_cause2) != "double" &&
        typeof(dat$p_t1_obs2) != "character" &&
        typeof(dat$p_t2_obs1) != "character") {
      
      data <- dat
      
      # Imputed probabilities for true causes (y1, y2)
      data$y1 <- as.numeric(data$r * data$c == 1) +
        (1 - as.numeric(data$r)) *
        ((dat$d1_obs == 1) * (1 - data$p_t2_obs1) +
           (dat$d2_obs == 1) * data$p_t1_obs2)
      
      data$y2 <- as.numeric(data$r * data$c == 2) +
        (1 - as.numeric(data$r)) *
        ((dat$d1_obs == 1) * data$p_t2_obs1 +
           (dat$d2_obs == 1) * (1 - data$p_t1_obs2))
      
      # Score functions and covariance matrices
      score1 <- estfun(model_cause1)
      score2 <- estfun(model_cause2)
      v1 <- vcov(model_cause1)
      v2 <- vcov(model_cause2)
      
      data_list[[i]] <- list(
        data = data,
        score1 = score1,
        score2 = score2,
        v1 = v1,
        v2 = v2
      )
    }
  }
  
  # Exclude invalid datasets (e.g., models that failed to converge)
  data_list <- data_list[lapply(data_list, typeof) == "list"]
  return(data_list)
}

# ------------------------------------------------------------------------------
# Function: data_list()
# ------------------------------------------------------------------------------
# Purpose:
#   Generate multiple simulated clustered datasets for replication or 
#   bootstrap studies.
#
# Input:
#   n           - number of datasets to generate
#   n_clusters - number of clusters per dataset
#   cause_num   - misclassification scenario indicator (1, 2, or 3)
#   ds          - double-sampling proportion (0–1)
#
# Output:
#   A list of simulated clustered datasets.
# ------------------------------------------------------------------------------

data_list <- function(n, n_clusters, cause_num, ds) {
  dat_list <- vector("list", n)
  for (i in 1:n) {
    dat_list[[i]] <- simulate(n_clusters = n_clusters, ds = ds, scen = cause_num)
  }
  return(dat_list)
}

