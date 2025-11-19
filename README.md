# **Semiparametric Marginal Regression for Clustered Competing Risks with Misclassified Causes of Failure**

R code accompanying my dissertation:
**“Semiparametric Marginal Regression for Clustered Competing Risks Data with Misclassified Causes of Failure.”**

This repository contains all functions and reproducible example code needed to:

1. Simulate clustered competing-risks data with informative cluster size (ICS)
2. Estimate time- and covariate-dependent misclassification probabilities using an external validation sample
3. Fit the proposed semiparametric marginal regression model that corrects for misclassification, clustering, and ICS
4. Conduct sensitivity analysis for transportability of the misclassification model
5. Perform cluster-bootstrap inference
6. Compute and plot cumulative incidence functions (CIFs) under different sensitivity assumptions

---

## **Repository Contents**

### **Core Estimation Functions**

* **`bssmle_clustered.R`**
  Implements the proposed semiparametric marginal regression estimator using B-splines, inverse cluster-size weighting, and externally estimated misclassification probabilities.

* **`pseudo_likelihood_estimation_Mpofu.R`**
  Pseudo-likelihood estimation of time- and covariate-dependent misclassification probabilities, extending Mpofu et al. (2020) to clustered designs.

### **Simulation and Example**

* **`simulate_data.R`**
  Generates clustered competing-risks data with frailty-based dependence, ICS, and uni- or bi-directional misclassification.

* **`example_analysis_clustered.Rmd`**
  A full end-to-end example demonstrating:

  * External validation simulation
  * Misclassification model estimation
  * Main analysis dataset generation
  * Application of the semiparametric estimator
  * Sensitivity analysis (η-grid)
  * Cluster-bootstrap inference
  * CIF reconstruction and plotting

* **`example_analysis_clustered.html`**
  Rendered output for easy viewing.

---

## **Method Overview**

### **1. External Validation & Misclassification Modeling**

Misclassification probabilities are estimated as:
$p_{21} = P(C^* = 2 \mid C = 1, T, Z)$, and $p_{12} = P(C^* = 1 \mid C = 2, T, Z)$ using a pseudo-likelihood framework with:

* logistic regression for predictive values
* logistic regression for misclassification
* cluster-level frailty and ICS incorporated in the sampling design

### **2. Semiparametric Marginal Regression Model**

The proposed estimator:

* uses B-splines to model baseline hazards
* adjusts for misclassification using ($p_{12}, p_{21}$)
* handles clustering and ICS via inverse cluster-size weighting
* yields marginal cause-specific regression effects
* provides interpretable hazard ratios and CIFs

### **3. Sensitivity Analysis (η-Shift)**

To assess violations of transportability from the validation sample, we use:

$logit(p_{jh}(η)) = logit(p_{jh}) + η$

with

$η ∈ [−0.5, −0.25, 0, 0.25, 0.5]$.


### **4. Cluster Bootstrap**

Entire clusters are resampled (with replacement) to compute:

* standard errors
* Wald statistics
* confidence intervals

### **5. CIF Computation**

Cumulative incidence functions are reconstructed from:

* spline-based baseline hazards
* estimated regression coefficients
* marginal survival and hazard components

---

## **How to Run the Example**

### **1. Clone or download the repository**

```bash
git clone https://github.com/fanismpal/misclassified-clustered-competing-risks.git
```

### **2. Open R or RStudio and source the scripts**

```r
source("simulate_data.R")
source("pseudo_likelihood_estimation_Mpofu.R")
source("bssmle_clustered.R")
```

### **3. Run the full analysis**

You can run the entire reproducible workflow by opening:

```
example_analysis_clustered.Rmd
```

and clicking **Knit**, or by viewing the pre-rendered:

```
example_analysis_clustered.html
```

---

## **Citation**

If you use this code in your research, please cite:

**Balanos, T. (2026).**
*Semiparametric Marginal Regression for Clustered Competing Risks Data with Misclassified Causes of Failure.*
Indiana University Indianapolis.

(Full journal citation will be added once published.)

---

## Contact

**Theofanis Balanos**
Department of Biostatistics & Health Data Science
Indiana University Indianapolis

Email: tbalanos@iu.edu

---
