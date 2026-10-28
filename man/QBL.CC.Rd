% Manually maintained help page following the layout of LBLGXE/man/LBL.Rd.
\name{QBL.CC}
\alias{QBL.CC}
\title{Quantitative Bayesian Lasso for Bivariate Haplotype Association}
\usage{
QBL.CC(dat, numSNPs, maxMissingGenos = 1, allelic = TRUE, haplo.baseline = "missing",
  cov.baseline = "missing", a = 20, b = 20, start.beta = 0.01, gamma = 0.01,
  lambda = 1, D = 0, e = 0.1, seed = NULL, burn.in = NULL, num.it = NULL,
  start.u = 0.01, sigma_sq_u = 1, sigma_sq_b = NULL, sigma_sq_c = NULL,
  e_allHap = 0.4, print.freq.ci = FALSE, print.lambda.ci = FALSE,
  print.D.ci = FALSE, print.sigma_sq_u.ci = FALSE,
  print.sigma_sq_b.ci = FALSE, print.sigma_sq_c.ci = FALSE)
}
\arguments{
\item{dat}{A data frame with one row per subject. The first two columns
contain numeric continuous phenotypes. The remaining columns contain SNP
genotypes: two allele columns per SNP when \code{allelic = TRUE}, or one
genotype column per SNP when \code{allelic = FALSE}. Environmental covariates
are not supported by this implementation.}

\item{numSNPs}{Number of SNPs per haplotype.}

\item{maxMissingGenos}{Maximum number of SNP genotypes with missing data
allowed per subject during preprocessing. Subjects exceeding this limit or
having missing non-SNP data are removed by \code{pre.hapassoc}.}

\item{allelic}{Logical. Whether SNPs are supplied in allelic format.
Missing alleles may be coded as \code{NA} or \code{""}. In genotypic format,
one observed allele can represent a partially missing genotype and
\code{""} a completely missing genotype.}

\item{haplo.baseline}{Haplotype label used as the reference. The default
\code{"missing"} selects the haplotype with the largest initial frequency
estimated by \code{pre.hapassoc}.}

\item{cov.baseline}{retained for compatibility with the existing interface. This option is ignored because environmental covariates are not included in this model.}

\item{a}{first hyperparameter of the gamma update for the shrinkage parameters, \code{lambda_b} and \code{lambda_c}; default is 20.}

\item{b}{rate hyperparameter used with \code{a} in the shrinkage updates; default is 20.}

\item{start.beta}{Initial value of the regression coefficients, including
the intercepts. Use a nonzero value because proposal scales depend on the
current coefficient magnitudes.}

\item{gamma}{retained for compatibility with the existing interface. This option is ignored by the current implementation.}

\item{lambda}{starting value of the \code{lambda_b} and \code{lambda_c} shrinkage parameters; default is 1.}

\item{D}{Initial value of the within-population inbreeding parameter.}

\item{e}{a small positive number epsilon defining the null hypothesis for a single haplotype: both phenotype coefficients have absolute values no greater than epsilon. The default is 0.1. The alternative is that at least one coefficient is outside this interval.}

\item{seed}{Optional seed passed to \code{set.seed} before MCMC sampling.
\code{NULL} leaves the current R random-number state unchanged at this step.}

\item{burn.in}{burn-in period of the MCMC sampling scheme. The default \code{NULL} uses 50000 iterations. Must be a nonnegative integer smaller than \code{num.it}.}

\item{num.it}{total number of MCMC iterations including burn-in. The default \code{NULL} uses 100000 iterations. Must be an integer greater than \code{burn.in}.}

\item{start.u}{starting value of the subject-specific latent variables, which induce correlation between the two phenotypes for each individual; default is 0.01. Use a nonzero value because proposal scales depend on the current magnitude.}

\item{sigma_sq_u}{Positive initial variance of the subject-specific latent
effects. The sampler uses its square root, \code{sigma_u}.}

\item{sigma_sq_b}{Positive initial residual variance for the first
phenotype. \code{NULL} uses the sample variance of the first column of
\code{dat}.}

\item{sigma_sq_c}{Positive initial residual variance for the second
phenotype. \code{NULL} uses the sample variance of the second column of
\code{dat}.}

\item{e_allHap}{epsilon defining the joint null hypothesis that all nonbaseline haplotype coefficients for both phenotypes have absolute values no greater than epsilon. The default is 0.4.}

\item{print.freq.ci}{Whether to calculate and return 95\% credible
intervals for haplotype frequencies.}

\item{print.lambda.ci}{Whether to calculate and return 95\% credible
intervals for both shrinkage parameters.}

\item{print.D.ci}{Whether to calculate and return a 95\% credible
interval for \code{D}.}

\item{print.sigma_sq_u.ci}{Whether to calculate and return a 95\% credible
interval for the latent-effect standard deviation \code{sigma_u}.
Despite the argument name, this interval is for the standard deviation,
not the variance.}

\item{print.sigma_sq_b.ci}{Whether to calculate and return a 95\% credible
interval for the first phenotype's residual variance.}

\item{print.sigma_sq_c.ci}{Whether to calculate and return a 95\% credible
interval for the second phenotype's residual variance.}
}
\value{
\item{BF_bivariate_hap}{A vector of Bayes factors for testing each nonbaseline haplotype with both phenotypes jointly. The alternative is association with at least one of the two phenotypes.}

\item{BF_bivariate_allHap}{The joint Bayes factor for testing all nonbaseline haplotypes in a block with both phenotypes jointly. The alternative is association of at least one haplotype with at least one phenotype.}

\item{beta_b}{A named vector of estimated posterior means of the regression coefficients, including the intercept, for the first continuous phenotype.}

\item{beta_c}{A named vector of estimated posterior means of the regression coefficients, including the intercept, for the second continuous phenotype.}

\item{CI.beta_b}{95\% credible sets for the first phenotype regression coefficients. The data frame has columns \code{Hap}, \code{Lower}, and \code{Upper}. These marginal intervals do not constitute a joint test of association with both phenotypes.}

\item{CI.beta_c}{95\% credible sets for the second phenotype regression coefficients, in the same format as \code{CI.beta_b}. These marginal intervals do not constitute a joint test of association with both phenotypes.}

\item{freq}{A named vector of posterior means of the haplotype frequencies, with the reference haplotype last.}

\item{CI.freq}{95\% credible sets for each haplotype frequency. Optional, returned only if \code{print.freq.ci = TRUE}. The data frame has columns \code{Hap}, \code{Lower}, and \code{Upper}. Otherwise, \code{""} is returned.}

\item{CI.lambda_b}{95\% credible set for \code{lambda_b}. Optional, returned only if \code{print.lambda.ci = TRUE}; otherwise, \code{""}.}

\item{CI.lambda_c}{95\% credible set for \code{lambda_c}. Optional, returned only if \code{print.lambda.ci = TRUE}; otherwise, \code{""}.}

\item{CI.D}{95\% credible set for \code{D}. Optional, returned only if \code{print.D.ci = TRUE}; otherwise, \code{""}.}

\item{CI.sigma_u}{95\% credible set for \code{sigma_u}, the square root of \code{sigma_sq_u}. Optional, returned only if \code{print.sigma_sq_u.ci = TRUE}; otherwise, \code{""}.}

\item{CI.sigma_sq_b}{95\% credible set for the residual variance of the first continuous phenotype. Optional, returned only if \code{print.sigma_sq_b.ci = TRUE}; otherwise, \code{""}.}

\item{CI.sigma_sq_c}{95\% credible set for the residual variance of the second continuous phenotype. Optional, returned only if \code{print.sigma_sq_c.ci = TRUE}; otherwise, \code{""}.}
}
\description{
Bayesian lasso is used to obtain posterior distributions of regression
coefficients for two correlated continuous phenotypes. The posterior
samples are used to calculate Bayes factors and credible sets for
haplotype association. Subject-specific latent variables model correlation
between the two phenotypes. Environmental covariates are not included.

The function first calls \code{pre.hapassoc} from the \pkg{hapassoc}
package. Arguments including \code{dat}, \code{numSNPs},
\code{maxMissingGenos}, and \code{allelic} are used in this preprocessing.
Rows of the input data frame correspond to subjects. The first two columns
are continuous phenotypes, followed by SNP genotypes. Subjects with too
many missing genotypes or with missing non-SNP data are removed.

Continuous phenotypes are centered and scaled after haplotype expansion.
Subjects with more compatible haplotype pairs therefore contribute more
rows to those scaling calculations. Regression summaries refer to this
transformed scale. Initial residual variances default to sample variances
from the original input columns. Internally, labels \code{b} and \code{c}
refer to the first and second continuous phenotypes, respectively.

Credible sets are equal-tailed 95\% posterior intervals. Coefficient means
are rounded to three decimal places and coefficient interval bounds to
two decimal places. Arguments beginning with \code{print.} control which
optional intervals are included in the returned list.

Short runs check execution but do not establish MCMC convergence or reliable
posterior inference. Raw draws and convergence diagnostics are not included
in the current return value. When the prior null probability is strictly
between zero and one, no saved draws inside the null region give an infinite
estimated Bayes factor. Such boundary estimates are not exact posterior
probabilities. The prior probabilities used for Bayes factors must match
the sampler's joint prior.
}
\examples{
# Load the example dataset: two continuous phenotypes and five SNPs.
data(QBL.ex1)

# Run QBL.CC to check that the package is loaded and runs correctly.
# Defaults use a longer chain. This short example is not for inference.
out.QBL <- QBL.CC(QBL.ex1, numSNPs = 5, allelic = TRUE,
                  seed = 123, burn.in = 100, num.it = 500)

out.QBL$beta_b
out.QBL$beta_c
}
\references{

Sajal IH and Biswas S (2023). Bivariate quantitative Bayesian LASSO for
 detecting association of rare haplotypes with two correlated continuous
 phenotypes. Frontiers in Genetics, 14:1104727.
 \url{https://doi.org/10.3389/fgene.2023.1104727}.

Yuan X and Biswas S (2021). Detecting rare haplotype association with two
 correlated phenotypes of binary and continuous types. Statistics in
 Medicine, 40(8):1877--1900. \url{https://doi.org/10.1002/sim.8877}.

Yuan X and Biswas S (2019). Bivariate logistic Bayesian LASSO for detecting
 rare haplotype association with two correlated phenotypes. Genetic
 Epidemiology, 43(8):996--1017. \url{https://doi.org/10.1002/gepi.22258}.

Biswas S and Lin S (2012). Logistic Bayesian LASSO for identifying
 association with rare haplotypes and application to age-related macular
 degeneration. Biometrics, 68(2):587--597.
 \url{https://doi.org/10.1111/j.1541-0420.2011.01680.x}.

Burkett K, Graham J and McNeney B (2006). hapassoc: Software for Likelihood
 Inference of Trait Associations with SNP Haplotypes and Other Attributes.
 Journal of Statistical Software, 16(2):1--19.
 \url{https://doi.org/10.18637/jss.v016.i02}.

}
\seealso{
\code{\link[hapassoc]{pre.hapassoc}}
}
\author{
Ibrahim Hossain Sajal, Swati Biswas, Shili Lin
}
