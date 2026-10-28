\docType{data}
\name{QBL.ex1}
\alias{QBL.ex1}
\title{Two Continuous Phenotypes and SNP Genotypes}
\format{
A data frame with 12 variables:
\describe{
  \item{Y1}{The first continuous phenotype.}
  \item{Y2}{The second continuous phenotype.}
  \item{X1, X2}{The two alleles at SNP 1.}
  \item{X3, X4}{The two alleles at SNP 2.}
  \item{X5, X6}{The two alleles at SNP 3.}
  \item{X7, X8}{The two alleles at SNP 4.}
  \item{X9, X10}{The two alleles at SNP 5.}
}
Rows correspond to subjects. The first two columns contain the phenotypes,
followed by ten allele columns representing five SNPs per haplotype.
}
\usage{
QBL.ex1
}
\description{
An example dataset containing two continuous phenotypes, no environmental
covariates, and SNP genotypes in allelic format. This dataset can be used
with \code{QBL.CC}, setting \code{numSNPs = 5} and \code{allelic = TRUE}.
}
\examples{
data(QBL.ex1)
str(QBL.ex1)
}
\seealso{
\code{\link{QBL.CC}}
}
\keyword{datasets}
