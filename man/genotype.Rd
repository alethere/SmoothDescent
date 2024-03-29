% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geno_predict.R
\name{genotype}
\alias{genotype}
\alias{genotype.matrix}
\alias{genotype.list}
\title{Genotype calculator}
\usage{
genotype(IBD, homologue, threshold = 0.8, ploidy = 2)

\method{genotype}{matrix}(IBD, homologue, threshold = 0.8, ploidy = 2)

\method{genotype}{list}(IBD, homologue, threshold = 0.8, ploidy = 2)
}
\arguments{
\item{IBD}{numeric matrix or list of matrices.}

\item{homologue}{homologue matrix for both parents. The first `ploidy` columns
from parent 1 and the second `ploidy` columns from parent 2.}

\item{threshold}{probability threshold to assign a homologue to an individual (
if IBD probability for hom 1 > threshold, it is considered that that individual has
inherited hom 1).}

\item{ploidy}{numeric}
}
\value{
A vector or matrix of imputed genotypes. If not enough IBD probabilities
are above the threshold for one individual, NA is returned.
}
\description{
Given a set of IBD probabilities, genotypes are calculated. Can take as input a single
matrix containing all homologue probabilities of a single individual; or a list of
matrices, where each element contains the probabilities of all individuals for one homologue
(same as output of `calc_IBD.list`)
}
\details{
Genotypes can be though of as a multiplication between a set of homologue probabilites
and homologue genotypes. For example, in a diploid we have the following IBD probabilities:
p(h1) = 0.9, p(h2) = 0.1, p(h3) = 1, p(h4) = 0; and genotypes g(h1) = 1, g(h2) = 0, g(h3) = 0, g(h4) = 1.
The resulting genotype using a simple multiplication would be \eqn{\sum{p(hi)*g(hi)}}, 0.9 in this case.\cr
We are more interested in imputing integer genotypes, so probabilities must be turned into
integers. For that reason we use a threshold, default 0.8. All probabilities above 0.8 are turned to 1,
and all probabilities below 0.8 are turned to 0. If not enough homologues are assigned this way,
that is, an individual cannot be assigned at least `ploidy` homologues, no genotype can be
calculated and NA is returned.
}
\section{Methods (by class)}{
\itemize{
\item \code{matrix}: Method for a single individual

\item \code{list}: Method for a list of IBD matrices
}}

\examples{
data("smooth_descent")
data("homologue")
IBD <- sdescent$predIBD

#For all individuals at once
gen <- genotype(IBD, hom)
#Or one single individual
#Note this is a matrix where each column is one homologue
ind <- sapply(IBD,'[',,1)
gen <- genotype(ind, hom)
}
