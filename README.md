# **Semiparametric Marginal Regression for Clustered Competing Risks Data with Misclassified Causes of Failure**

R code accompanying the paper:

*"Semiparametric Marginal Regression for Clustered Competing Risks Data with Misclassified Causes of Failure."*

This repository contains all functions and reproducible example code needed to:

1. Simulate clustered competing-risks data with within-cluster dependence and informative cluster size (ICS).
2. Estimate time- and covariate-dependent misclassification probabilities from an external validation sample, with clustering and ICS accounted for via inverse cluster-size weighting and GEE with a working independence correlation structure and a sandwich variance estimator.
3. Fit the proposed semiparametric marginal regression model that simultaneously corrects for misclassification of the cause of failure, within-cluster dependence, and ICS.
4. Conduct sensitivity analyses to assess robustness against violations of the transportability assumption.
5. Perform cluster-bootstrap inference that also accounts for the uncertainty in the externally estimated misclassification probabilities.
6. Compute and plot cumulative incidence functions (CIFs) under each sensitivity-analysis setting.

---

## **Repository contents**

### **Core estimation functions**

* **`bssmle_clustered.R`**
  Implements the proposed semiparametric marginal regression estimator using a B-spline-based sieve pseudo-likelihood with inverse cluster-size weighting and externally estimated misclassification probabilities. Returns baseline cause-specific hazard spline coefficients and marginal cause-specific regression effects.

* **`pseudo_likelihood_estimation_Mpofu.R`**
  Estimates time- and covariate-dependent misclassification probabilities using the pseudo-likelihood approach of Mpofu et al. (2020), extended here to clustered designs through inverse cluster-size weighting and GEE with a working independence correlation structure and a sandwich variance estimator.

### **Simulation and example**

* **`simulate_data.R`**
  Generates clustered competing-risks data with frailty-based within-cluster dependence, informative cluster size, and uni- or bi-directional misclassification. Supports both external-validation and main-study simulation.

* **`example_analysis_clustered.Rmd`**
  A fully reproducible end-to-end example demonstrating:

  * External validation data simulation
  * Misclassification model estimation under clustering and ICS
  * Main analysis dataset generation
  * Application of the proposed semiparametric marginal estimator
  * Sensitivity analysis over an η-grid
  * Cluster-bootstrap inference that incorporates uncertainty in the estimated misclassification parameters
  * Reconstruction and plotting of CIFs

* **`example_analysis_clustered.html`**
  Rendered output of the full example for easy viewing.

---

## **Method overview**

### **1. External validation and misclassification modeling**

Misclassification probabilities are estimated as:

$p_{21} = P(C^* = 2 \mid C = 1, T, Z)$, and $p_{12} = P(C^* = 1 \mid C = 2, T, Z)$,

using a pseudo-likelihood framework with:

* logistic regression for predictive values,
* logistic regression for misclassification,
* inverse cluster-size weighting to address ICS,
* GEE with a working independence correlation structure and a sandwich variance estimator to account for within-cluster dependence.

The estimated coefficients $`\hat{\gamma}_{n'}`$ and their variance matrix $`\hat{\Omega}_{n'}`$ are then carried into the main analysis to compute predicted misclassification probabilities for each subject.

---

### **2. Semiparametric marginal regression model**

The proposed estimator:

* uses B-splines to model baseline cause-specific cumulative hazards,
* adjusts for misclassification using the externally estimated probabilities $(p_{12}, p_{21})$,
* handles within-cluster dependence and ICS via inverse cluster-size weighting,
* fits a marginal proportional cause-specific hazards model,
* simultaneously models both causes through a unified pseudo-likelihood,
* produces interpretable population-averaged hazard ratios.

---

### **3. Sensitivity analysis (η-shift)**

To examine violations of the transportability assumption, misclassification probabilities in the main dataset are adjusted via

$\text{logit}(p_{jh}(\eta)) = \text{logit}(p_{jh}) + \eta$, with $`\eta \in \{-0.5, -0.25, 0, 0.25, 0.5\}`$.

Each modified misclassification scenario yields a new set of regression estimates, allowing assessment of robustness when the misclassification mechanism in the main study differs from that in the validation study.

The conceptual framework, motivation, and detailed derivation of this η-shift sensitivity analysis are described in the companion paper:

> Balanos T, Yiannoutsos CT, Pabon-Rodriguez FM, Nan H, Bakoyannis G. *Semiparametric Regression for Misclassified Competing Risks Data.* arXiv preprint arXiv:2605.16652. 2026. <https://doi.org/10.48550/arXiv.2605.16652>
>
> Code repository: <https://github.com/tbalanos/misclassified-competing-risks-bssmle>

Readers interested in the underlying rationale for this sensitivity analysis are encouraged to consult that paper.

---

### **4. Cluster bootstrap with uncertainty in misclassification parameters**

Inference for the main regression model is carried out using a clustered bootstrap that accounts for two sources of variability:

* **Sampling variability** in the main study, addressed by resampling entire clusters (with replacement).
* **Uncertainty in the externally estimated misclassification parameters** $\hat{\gamma}_{n'}$, addressed by drawing a new realization $\tilde{\gamma}_{n'}^{(b)} \sim N(\hat{\gamma}_{n'}, \hat{\Omega}_{n'})$ at each bootstrap iteration, recomputing the misclassification probabilities, and then re-fitting the model.

This combined procedure yields:

* bootstrap standard errors,
* Wald statistics,
* confidence intervals.

Importantly, the procedure does **not** require access to the original external validation data: only the estimate $\hat{\gamma}_{n'}$ and its variance matrix $\hat{\Omega}_{n'}$ are needed.

---

### **5. CIF computation**

Cumulative incidence functions are reconstructed using:

* the spline-based baseline cumulative hazards,
* the estimated marginal cause-specific regression coefficients,
* the resulting cause-specific hazards and overall survival function.

---

## **How to run the example**

### **1. Clone or download the repository**

```bash
git clone https://github.com/tbalanos/misclassified-clustered-competing-risks.git
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

Balanos T, Yiannoutsos CT, Bakoyannis G. *Semiparametric Marginal Regression for Clustered Competing Risks Data with Misclassified Causes of Failure.* 2026.

(Full journal citation will be added once published.)

### **Related work**

Balanos T, Yiannoutsos CT, Pabon-Rodriguez FM, Nan H, Bakoyannis G. *Semiparametric Regression for Misclassified Competing Risks Data.* arXiv preprint arXiv:2605.16652. 2026. <https://doi.org/10.48550/arXiv.2605.16652>

Code repository: <https://github.com/tbalanos/misclassified-competing-risks-bssmle>

---

## Contact

**Theofanis Balanos, Ph.D.**
Department of Biostatistics and Health Data Science
Richard M. Fairbanks School of Public Health
Indiana University Indianapolis

Email: **[tbalanos@iu.edu](mailto:tbalanos@iu.edu)**

---
