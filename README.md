# QBL

**Quantitative Bayesian Lasso for Bivariate Haplotype Association**

QBL provides `QBL.CC()` for haplotype association analysis with two correlated continuous phenotypes. It uses Markov chain Monte Carlo (MCMC) sampling to obtain regression coefficient and haplotype frequency summaries, credible intervals, and Bayes factors for individual haplotypes and the haplotype block jointly.

The current model does not include environmental covariates.

## Installation

Install QBL from GitHub:

```r
install.packages("remotes")
remotes::install_github("IHSajal/QBL", upgrade = "never")
```

The package declares R >= 3.5.0 and depends on `hapassoc` and `stats`. It contains C source code, so installation requires suitable compilation tools. On macOS, these include the Command Line Tools for Xcode.

## Example

`QBL.ex1` contains two continuous phenotypes and allelic genotype data for five SNPs. Each row represents a subject. The first two columns are the phenotypes, followed by two allele columns per SNP.

```r
library(QBL)
data("QBL.ex1", package = "QBL")

# Short demonstration to check execution; not a run for inference.
fit <- QBL.CC(
  dat = QBL.ex1,
  numSNPs = 5,
  allelic = TRUE,
  seed = 123,
  burn.in = 100,
  num.it = 500
)

fit$beta_b
fit$beta_c
fit$BF_bivariate_hap
fit$BF_bivariate_allHap
```

The labels `b` and `c` refer to the first and second continuous phenotypes, respectively. Regression summaries are on the transformed outcome scale described in `?QBL.CC`.

The default chain has 100,000 iterations, including 50,000 burn-in iterations. These defaults do not by themselves establish convergence. The current function returns posterior summaries rather than raw MCMC chains. Short-run results, including zero or infinite estimated Bayes factors, should not be used for scientific conclusions.

## Optional credible intervals

Set the following arguments to `TRUE` to include additional 95% credible intervals in the returned list:

- `print.freq.ci`: haplotype frequencies.
- `print.lambda.ci`: both shrinkage parameters.
- `print.D.ci`: the inbreeding parameter.
- `print.sigma_sq_u.ci`: the latent-effect **standard deviation**, returned as `CI.sigma_u`.
- `print.sigma_sq_b.ci`: the first phenotype's residual variance.
- `print.sigma_sq_c.ci`: the second phenotype's residual variance.

## Documentation

```r
help("QBL", package = "QBL")
help("QBL.CC", package = "QBL")
help("QBL.ex1", package = "QBL")
```

## Package checks

An initial local check on R 4.5.1 and Apple silicon macOS completed with **0 errors, 1 warning, and 0 notes**. The warning originated in R's `R_ext/Boolean.h` header when compiling with Apple clang 21. Documentation and example checks passed. This records the local packaging check; validation of the sampler and Bayes-factor calculations is a separate task.

## Authors and maintainer

Authors: Ibrahim Hossain Sajal, Swati Biswas, and Shili Lin.

Maintainer: Ibrahim Hossain Sajal (<ihossain2@isrt.ac.bd>).

## References

- Sajal IH and Biswas S (2023). Bivariate quantitative Bayesian LASSO for detecting association of rare haplotypes with two correlated continuous phenotypes. *Frontiers in Genetics*, 14:1104727. [doi:10.3389/fgene.2023.1104727](https://doi.org/10.3389/fgene.2023.1104727).
- Yuan X and Biswas S (2021). Detecting rare haplotype association with two correlated phenotypes of binary and continuous types. *Statistics in Medicine*, 40(8):1877–1900. [doi:10.1002/sim.8877](https://doi.org/10.1002/sim.8877).
- Yuan X and Biswas S (2019). Bivariate logistic Bayesian LASSO for detecting rare haplotype association with two correlated phenotypes. *Genetic Epidemiology*, 43(8):996–1017. [doi:10.1002/gepi.22258](https://doi.org/10.1002/gepi.22258).
- Biswas S and Lin S (2012). Logistic Bayesian LASSO for identifying association with rare haplotypes and application to age-related macular degeneration. *Biometrics*, 68(2):587–597. [doi:10.1111/j.1541-0420.2011.01680.x](https://doi.org/10.1111/j.1541-0420.2011.01680.x).
- Burkett K, Graham J and McNeney B (2006). hapassoc: Software for Likelihood Inference of Trait Associations with SNP Haplotypes and Other Attributes. *Journal of Statistical Software*, 16(2):1–19. [doi:10.18637/jss.v016.i02](https://doi.org/10.18637/jss.v016.i02).

## License

GPL-3.
