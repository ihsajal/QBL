\name{QBL-package}
\alias{QBL-package}
\alias{QBL}
\docType{package}
\title{
Quantitative Bayesian Lasso for Bivariate Haplotype Association
}
\description{
The main function of this package is \code{QBL.CC}. For details, see
\code{\link{QBL.CC}}.
}
\details{
\tabular{ll}{
Package: \tab QBL\cr
Type: \tab Package\cr
Version: \tab 1.0.0\cr
License: \tab GPL-3\cr
NeedsCompilation: \tab yes\cr
}

Quantitative Bayesian lasso is used for haplotype association analysis
with two correlated continuous phenotypes. The package takes a dataset
containing the two phenotypes and SNP genotypes, with no environmental
covariates. Markov chain Monte Carlo sampling is used to estimate posterior
distributions of regression coefficients and haplotype frequencies.
Outputs include posterior means, credible sets, and Bayes factors for
individual haplotype and overall haplotype association.

Currently available analysis function: \code{QBL.CC}. Type
\code{?QBL.CC} for details.

Example dataset: \code{\link{QBL.ex1}}. It contains two continuous
phenotypes and allelic genotype data for five SNPs.
}
\author{
Ibrahim Hossain Sajal, Swati Biswas, Shili Lin

Maintainer: Ibrahim Hossain Sajal \email{ihossain2@isrt.ac.bd}
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
\code{\link{QBL.CC}}, \code{\link{QBL.ex1}},
\code{\link[hapassoc]{pre.hapassoc}}
}
\keyword{package}
\examples{
# See ?QBL.CC for examples of haplotype association analysis.
}
