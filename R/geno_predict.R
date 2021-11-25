#3 Genotype prediction -----------------

#Function to calculate genotypes from a single individual
#ibd matrix and a homologue matrix

#' Genotype calculator
#'
#' Given a set of IBD probabilities, genotypes are calculated. Can take as input a single
#' matrix containing all homologue probabilities of a single individual; or a list of
#' matrices, where each element contains the probabilities of all individuals for one homologue
#' (same as output of `calc_IBD.list`)
#'
#' @param IBD numeric matrix or list of matrices.
#' @param homologue homologue matrix for both parents. The first `ploidy` columns
#' from parent 1 and the second `ploidy` columns from parent 2.
#' @param threshold probability threshold to assign a homologue to an individual (
#' if IBD probability for hom 1 > threshold, it is considered that that individual has
#' inherited hom 1).
#' @param ploidy numeric
#'
#' @details Genotypes can be though of as a multiplication between a set of homologue probabilites
#' and homologue genotypes. For example, in a diploid we have the following IBD probabilities:
#' p(h1) = 0.9, p(h2) = 0.1, p(h3) = 1, p(h4) = 0; and genotypes g(h1) = 1, g(h2) = 0, g(h3) = 0, g(h4) = 1.
#' The resulting genotype using a simple multiplication would be \eqn{\sum{p(hi)*g(hi)}}, 0.9 in this case.\cr
#' We are more interested in imputing integer genotypes, so probabilities must be turned into
#' integers. For that reason we use a threshold, default 0.8. All probabilities above 0.8 are turned to 1,
#' and all probabilities below 0.8 are turned to 0. If not enough homologues are assigned this way,
#' that is, an individual cannot be assigned at least `ploidy` homologues, no genotype can be
#' calculated and NA is returned.
#'
#' @return A vector or matrix of imputed genotypes. If not enough IBD probabilities
#' are above the threshold for one individual, NA is returned.
#' @export
#'
#' @examples
#' data("smooth_descent")
#' data("homologue")
#' IBD <- sdescent$predIBD
#'
#' #For all individuals at once
#' gen <- genotype(IBD, hom)
#' #Or one single individual
#' #Note this is a matrix where each column is one homologue
#' ind <- sapply(IBD,'[',,1)
#' gen <- genotype(ind, hom)
genotype <- function(IBD,homologue,threshold = 0.8, ploidy = 2){
  UseMethod("genotype",IBD)
}

#' Genotype calculator for a single individual
#'
#' @describeIn genotype Method for a single individual
genotype.matrix <- function(IBD,homologue,threshold = 0.8,ploidy = 2){
  # #Above threshold are turned into 1
  # IBD[IBD >= threshold] <- 1
  # IBD[IBD <= 1-threshold] <- 0

  marks <- intersect(rownames(homologue),rownames(IBD))
  IBD <- IBD[marks,]
  homologue <- homologue[marks,]
  #Genotypes are calculated
  perhom_geno <- IBD*homologue
  geno <- rowSums(perhom_geno)
  not_certain <- rowSums(IBD)/ploidy < threshold

  #There must be two columns per individual over the threshold
  #This is somehow problematic though, what if two columns of p1 instead of 1 in p1 and one in p2?
  # not_certain <- rowSums(IBD > threshold) < ploidy
  # not_certain <- (perhom_geno != 0 & perhom_geno != 1)
  # not_certain <- apply(not_certain, 1, any) | rowSums(IBD) != ploidy
  # decimals <- geno %% 1
  # not_certain <- decimals > threshold | decimals < threshold
  geno[not_certain] <- NA
  geno <- round(geno)

  return(geno)
}

#' Genotype calculator for lists
#'
#' @describeIn genotype Method for a list of IBD matrices
genotype.list <- function(IBD,homologue,threshold = 0.8, ploidy = 2){
  if(!all(colnames(homologue) %in% names(IBD))) stop("IBD list should contain one element per column in homologue matrix")

  ind_ibd <- lapply(1:ncol(IBD[[1]]),function(i){
    ib <- sapply(IBD,'[',,i)
  })
  geno <- sapply(ind_ibd,genotype.matrix,homologue = homologue,
                 threshold = threshold, ploidy = ploidy)
  colnames(geno) <- colnames(IBD[[1]])
  return(geno)
}
