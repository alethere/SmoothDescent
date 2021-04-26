#IBD calculation -----------------
#This file contains a generic and the methods for performing naïve
#IBD calculation regardless of ploidy in an F1.

#' Calculation of IBD
#'
#' @param dos integer or matrix containing genotypes (dosages of allele A, from 0 to ploidy).
#' If `dos` is a matrix, each row must be labelled with a marker name, and each column is
#' interpreted as an individual.
#' @param p1hom numeric vector or matrix indicating which parental homologue, from parent 1, contains allele A.
#' If `dos` is an integer, `p1hom` must be a vector of length ploidy. If `dos` is a matrix,
#' `p1hom` must be a matrix of ncol = ploidy and nrow = nrow(dos) and the same rownames
#' as `dos`.
#' @param p2hom numeric vector or matrix indicating which parental homologue, from parent 2, contains allele A.
#' If `dos` is an integer, `p2hom` must be a vector of length ploidy. If `dos` is a matrix,
#' `p2hom` must be a matrix of ncol = ploidy and nrow = nrow(dos) and the same rownames
#' as `dos`.
#' @param ploidy numeric, indicating the ploidy of an organism. Only even ploidies allowed and
#' all individuals are expected to have the same ploidy.
#'
#' @return The naïve probability of having inherited each parental homologue. If the
#' offspring dosage is impossible given the parental homologues, NA is returned.\n
#' If `dos`is a vector, returns a one-row matrix, each column containing the
#' probability of having inherited that parental homologue (the first columns corresponding
#' to parent 1, and the rest to parent 2). \n If `dos` is a matrix, a list of matrices is returned.
#' Each element is the probabilities for all individuals of having inherited one parental
#' homologue. There should be ploidy*2 elements in the list.
#' @export
#'
#' @examples
calc_IBD <- function(dos,p1hom,p2hom,ploidy = 2){
  dos_p1 <- identical(class(dos),class(p1hom))
  dos_p2 <- identical(class(dos),class(p2hom))
  if(!(dos_p1 & dos_p2)){
    stop("All data provided to IBD must be of the same class (i.e. all numeric or all matrices)")
  }else{
    UseMethod("calc_IBD",dos)
  }
}



#' Calculation of IBD for a single data point
#'
#' @describeIn calc_IBD Function for a single marker of a single individual.
#' @export
calc_IBD.numeric <- function(dos,p1hom,p2hom,ploidy = 2){
  if(all(is.na(dos))) return(NA)
  assertthat::assert_that(all(dos <= ploidy))
  assertthat::assert_that(length(p1hom) == ploidy)
  assertthat::assert_that(length(p2hom) == ploidy)

  configlist <- config(dos,ploidy = ploidy)
  parenthom <- list(p1 = p1hom,
                    p2 = p2hom)
  parentcols <- list(p1 = 1:(ploidy/2),
                     p2 = (ploidy/2+1):ploidy)
  #This expresses the maximum and minimum allele dosage
  #that each parent can give. It's useful to find out
  #impossible inheritances
  parentrange <- sapply(parenthom,function(homp){
    min <- sum(sort(homp)[1:(ploidy/2)])
    max <- sum(sort(homp,decreasing = T)[1:(ploidy/2)])
    return(c(min = min,max = max))
  })
  #We will calculate the homologue probability
  #per each configuration
  config_prob <- lapply(configlist,function(con){

    #This calculates the homologue probabilites
    #per parent for the configuration at which we're at
    homprobs <- lapply(c(1,2),function(parent){
      conp <- con[parentcols[[parent]]]
      homp <- parenthom[[parent]]

      if(any(is.na(homp)) | sum(conp) > parentrange["max",parent] | sum(conp) < parentrange["min",parent]){
        impossible <- T
      }else{
        impossible <- F
      }

      homprobs <- t(sapply(conp,function(inherited){
        sapply(homp,function(hom){
          if(impossible) return(NA)
          res <- (inherited == hom)/sum(homp == inherited)
          if(is.nan(res)) res <- 0
          return(res)
        })
      }))

      rownames(homprobs) <- conp
      return(homprobs)
    })

    #Here something else needs to be added in order to make it
    #for polyploids
    homprobs <- do.call(cbind,homprobs)
    homprobs <- colSums(homprobs)

    #probability of parent 1 to give first genotypes
    probp1 <- sapply(con[parentcols$p1],function(gamete) mean(parenthom$p1 == gamete) )
    probp1 <- prod(probp1,na.rm = T)
    #probability of parent 2 to give other genotypes
    probp2 <- sapply(con[parentcols$p2],function(gamete) mean(parenthom$p2 == gamete) )
    probp2 <- prod(probp2,na.rm = T)

    return(list(conprob = probp1*probp2,
                homprob = homprobs))

  })

  #Now each configuration probability
  #First the probability of each configuration
  conprob <- extract(config_prob,"conprob")
  homprob <- extract(config_prob,"homprob")
  homprob <- matrix(homprob,nrow = length(p1hom) + length(p2hom))
  possible_configs <- !is.na(colSums(homprob))
  conprob <- conprob[possible_configs]
  homprob <- homprob[,possible_configs,drop = F]
  #Now we normalize it (is this correct?)
  conprob <- conprob/sum(conprob,na.rm = T)
  #Then we multiply the probability of the configuration times the
  #probability of each homologue. Summing the homologue probabilities
  #provides us with the final homologue probabilities
  homprob <- t(homprob %*% conprob)
  colnames(homprob) <- c(names(p1hom),names(p2hom))

  if(all(!possible_configs)) return(matrix(NA,ncol = length(p1hom) + length(p2hom)))
  return(homprob)
}

#IBD probability calculator for a matrix of markers x individuals.
#dos is a numeric matrix with dosages, rownames and columnames.
#p1hom is a homologue matrix for parent 1. NAs allowed.
#p2hom is a homologue matrix for parent 2. NAs allowed.
#Returns a list of IBD matrices for each parental homologue

#' Calculation of IBD for multiple markers and/or individuals
#'
#' @describeIn calc_IBD Function for handling matrix input to calc_IBD
#' @export
calc_IBD.matrix <- function(dos,p1hom,p2hom,ploidy = 2){
  if(is.vector(dos)){
    nams <- names(dos)
    dos <- matrix(dos,ncol = 1)
    rownames(dos) <- nams
  }
  if(is.null(rownames(dos))) stop("Dosage matrix should contain rownames")
  if(is.null(rownames(p1hom))) stop("P1 homologue matrix should contain rownames")
  if(is.null(rownames(p2hom))) stop("P2 homologue matrix should contain rownames")
  markers_p1 <- rownames(dos) %in% rownames(p1hom)
  if(!all(markers_p1)) stop("Some rownames in dosage matrix were not found in P1 homologue matrix: ",
                            paste(rownames(dos)[!markers_p1],collapse = " "))
  markers_p2 <- rownames(dos) %in% rownames(p2hom)
  if(!all(markers_p2)) stop("Some rownames in dosage matrix were not found in P2 homologue matrix:\n",
                            paste(rownames(dos)[!markers_p2],collapse = " "))
  assertthat::assert_that(ncol(p1hom) == ploidy)
  assertthat::assert_that(ncol(p2hom) == ploidy)

  #To save time, we only calculate each unique case
  #To identify unique cases we transform the data into one long vector
  #This is memory expensive, might need to rework for big individual x genotype matrices
  nmark <- nrow(dos)
  nind <- ncol(dos)
  dosvec <- as.vector(dos)
  names(dosvec) <- rep(rownames(dos),nind)
  p1hom <- p1hom[names(dosvec),]
  p2hom <- p2hom[names(dosvec),]

  #First identify the number of combinations
  data <- data.frame(dos = dosvec, p1hom, p2hom)
  parentcols <- list(p1 = 2:(ncol(p1hom)+1),
                     p2 = (ncol(p1hom)+2):ncol(data))
  comb <- apply(data,1,paste,collapse = "_")
  marks <- 1:nrow(data)
  cases <- split(marks,comb)

  #Then we fill in a result matrix
  #we only calculate one time each combination
  res <- matrix(0,ncol = ncol(data) -1,nrow = nrow(data))
  rownames(res) <- names(dosvec)
  for(i in 1:length(cases)){
    cas <- names(cases)[i]
    dat <- suppressWarnings(as.numeric(strsplit(cas,"_")[[1]]))
    prob <- calc_IBD.numeric(dos = dat[1],p1hom = dat[parentcols$p1],
                             p2hom = dat[parentcols$p2],ploidy = ploidy)
    res[cases[[i]],] <- rep(prob,each = length(cases[[i]]))
  }

  #We then reform the output into one homologue list.
  homlist <- lapply(1:ncol(res),function(i){
    x <- res[,i]
    res <- matrix(x,ncol = nind, nrow= nmark)
    rownames(res) <- rownames(dos)
    colnames(res) <- colnames(dos)
    return(res)
  })
  names(homlist) <- c(colnames(p1hom),colnames(p2hom))
  return(homlist)
}

#Small function used in calc_IBD.numeric. It returns a list of the possible
#genotype configurations given a dosage (g) and a ploidy. For instance
#config(1,ploidy = 2) returns 1,0 and 0,1
#' Configuration function
#'
#' Small internal function that returns all possible configurations of
#' parental inheritance given a genotype
#'
#' @param g integer, a genotype <= ploidy and >=0
#' @param ploidy integer, indicating ploidy.
#'
#' @return a list of vectors, where each element indicates
#' one possible parental chromosome inheritance configuration.
#' @example
#' config(1,2)
config <- function(g,ploidy = 2){
  doses <- c(rep(1,g),rep(0,ploidy-g))
  return(unique(combinat::permn(doses)))
}

#Test --------
source("R/Utils.R")
dos <- as.matrix(read.table("test/test_geno.txt",header = T))
dos <- dos[,-1:-2] #we take out parental genotypes
hom <- as.matrix(read.table("test/test_hom.txt",header = T))
#Numeric
p2_test <- calc_IBD(dos[1,1],hom[1,1:2],hom[1,3:4],ploidy = 2)

if(!is.matrix(p2_test)){
  stop("Output of calc_IBD.numeric is not matrix")
}
if(!ncol(p2_test) == 2*2){
  stop("Length of calc_IBD.numeric output is not ploidy*2")
}
if(!all(p2_test == c(0.5,0.5,1,0))){
  stop("Output of calc_IBD.numeric does not coincide with saved output")
}

#Matrix
ploidy = 2
matrix_test <- calc_IBD(dos,p1hom = hom[,1:2],p2hom = hom[,3:4],ploidy = 2)
#saveRDS(matrix_test,"test/calc_IBD.matrix.RDS")


if(!is.list(matrix_test)){
  stop("Output of calc_IBD.matrix is not list")
}
if(!length(matrix_test) == ploidy*2){
  stop("Length of calc_IBD.matrix output is not ploidy*2")
}
if(!all(sapply(matrix_test,ncol) == ncol(dos))){
  stop("Not all elements from output of calc_IBD.matrix have ncol == ncol(dos)")
}
saved_mtest <- readRDS("test/calc_IBD.matrix.RDS")
if(!identical(matrix_test,saved_mtest)){
  stop("Output of calc_IBD.matrix does not coincide with saved output")
}
