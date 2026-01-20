# zibFit

<!-- badges: start -->

<!-- badges: end -->

The goal of zibFit is to make it easier to model with zero-inflated beta regression models.

## Introduction

The Zero-Inflated Beta Regression (ZIBR) model, as proposed by Eric Z. Chen and Hongzhe Li in 2016, is a proposed solution for modelling longitudinal microbiome compositional data. This data tends to be highly skewed, is bounded between [0, 1) and is often sparse, with many zeroes. Additionally, observations from repeated measurements are correlated. The ZIBR model proposes a two part zero-inflated beta regression model with random effects for testing the association between microbial abundance and clinical covariates for longitudinal microbiome data. The model includes a logistic component to model presence/absence of the microbe in samples and a Beta component to model non-zero microbial abundance. Both components include a random effect to take into account the correlation among repeated measurements on the same subject.

## Statistical Model

The details of the statistical model are as follows:

ZIBR assumes $y_{it}$, the relative abundance of bacterial taxon in individual $i$ at time $t$, for $1 \leq i \leq N, 1 \leq t \leq{T_i}$, follows the following distribution; 
$\begin{align*}
    y_{it} \sim 
    \begin{cases}
        0 & \text{with probability }1-p_{it}\\
        Beta(\mu_{it}\phi, (1-\mu_{it})\phi) & \text{with probability }p_{it}
    \end{cases}
\end{align*}$
Where $\phi > 0$ and $0 < \mu_{it}, p_{it} < 1$.

Accordingly, we can test three biologically relevant null hypotheses:\
- H0: α_j = 0. This is to test the coefficients in the logistic component, if the covariates are associated with the bacterial taxon by affecting its presence or absence;\
- H0: β_j = 0. This is to test the coefficients in the Beta component, if the taxon is associated with the covariates by showing different abundances;\
- H0: α_j = 0 and β_j = 0 for each covariate X_j and Z_j. This is to joinly test the coefficients in both logistic and Beta components, if the covariates affect the taxon both in terms of presence/absence and its non-zero abundance.

The zibFit package features 4 functions to assist in modelling longitudinal microbiome composition data; zibSim (to simulate data), zibPlot (plot method for plotting relative abundance), zibClean (to aid in data cleaning prior to implementing the zibFitter function) and zibFitter (which fits the ZIBR model to a microbiome dataset).

## Installation

You can install the development version of zibFit from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("joe-tanner/zibFit")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(zibFit)
# simulate zibr data and plot

# initialise covariate matrix, simulating a clinical study where half recieve the treatment (1) and half do not (0)
n = 1000
t = 3
X <- as.matrix(c(rep(0, n/2 * t), rep(1, n/2 * t))) 

y1 <- zibSim(n = n, t = t, a = 0.5, b = 0.8, sigma1 = 2.5, sigma2 = 0.5, phi = 7.2, alpha = 0.5, beta = -0.5, X = X, Z = X, seed = 6874)

# plot the data
plot(y1, microbe_name = "Simulated Microbe")

## fit a zibr model to longitudinal microbiome composition data
# simulate five different microbiomes
y2 <- zibSim(n = n, t = t, a = 0.5, b = -0.5, sigma1 = 0.8, sigma2 = 1.2, phi = 5, alpha = -0.5, beta = -0.5, X = X, Z = X, seed = 6874)
y3 <- zibSim(n = n, t = t, a = 0.2, b = 0.3, sigma1 = 3, sigma2 = 2, phi = 4, alpha = -0.5, beta = -0.5, X = X, Z = X, seed = 6874)
y4 <- zibSim(n = n, t = t, a = 0.1, b = 0.9, sigma1 = 0.4, sigma2 = 0.5, phi = 8.1, alpha = 0.5, beta = 0.5, X = X, Z = X, seed = 6874)
y5 <- zibSim(n = n, t = t, a = 0.2, b = 1, sigma1 = 2, sigma2 = 0.5, phi = 3.4, alpha = -0.5, beta = 0.5, X = X, Z = X, seed = 6874)

ra <- cbind(y1$rel_abundance, y2$rel_abundance, y3$rel_abundance, y4$rel_abundance, y5$rel_abundance)
rownames(ra) <- rep(1:3000) # rownames and column names are required
colnames(ra) <- rep(1:5)

# create covariate data frame
sample_id <- rep(1:3000)
covariates <- data.frame(subject = y1$subject_ind, time = y1$time_ind, treat = y1$log_covariates$X, sample = sample_id)

zibData = zibClean(data = ra, metadata = covariates, log_covs = c("subject", "time", "treat"), id_column = "sample", subject_column = "subject", time_column = "time")
# if working on a real world problem, just use zibClean with your relative abundance and metadata files in the correct format!

# results list now includes p-values for all estimates, as well as lists of each bacteria with statistically significant treatment, subject, time and baseline abundance effects
results <- zibFitter(zibData)


```

## Citation

Eric Z. Chen and Hongzhe Li (2016). A two-part mixed effect model for analyzing longitudinal microbiome data. Bioinformatics. [Link](https://academic.oup.com/bioinformatics/article-abstract/32/17/2611/2450750?rss=1)
