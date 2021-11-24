#IBD calculation -----------------

#' Calculation of IBD
#'
#' This function computes Identity-by-descent probabilities according to a naive model
#' based on random assortment of markers. It requires homologue matrices (`p1hom` and `p2hom`)
#' indicating which homologues of each parent carry a genoage of 1 or 0 for the tracked allele.
#' Based on this, a probability matrix is returned: for each individual how likely it is,
#' at each marker, that it has inherited each of the homologues.
#' As a generic method, it can use a single genoage of one marker (and two vectors
#' of homologue assignment); a vector of markers (and two homologue matrices) or a matrix of markers
#' (and two homologue matrices).
#'
#' @param geno integer or matrix containing genotypes (genoages of allele A, from 0 to ploidy).
#' If `geno` is a matrix, each row must be labelled with a marker name, and each column is
#' interpreted as an individual. Note that for many markers of one single individual, a
#' 1 column matrix must be provided.
#' @param p1hom numeric vector or matrix indicating which parental homologue, from parent 1, contains allele A.
#' If `geno` is an integer, `p1hom` must be a vector of length ploidy. If `geno` is a matrix, `p1hom` must be a matrix of ncol = ploidy and nrow = nrow(geno) and the same rownames
#' as `geno`.
#' @param p2hom numeric vector or matrix indicating which parental homologue, from parent 2, contains allele A.
#' If `geno` is an integer, `p2hom` must be a vector of length ploidy. If `geno` is a matrix,`p2hom` must be a matrix of ncol = ploidy and nrow = nrow(geno) and the same rownames
#' as `geno`.
#' @param ploidy numeric, indicating the ploidy of an organism. Only even ploidies allowed and
#' all individuals are expected to have the same ploidy.
#'
#' @return The naive probability of having inherited each parental homologue. If the
#' offspring genoage is impossible given the parental homologues, NA is returned.
#' If `geno`is a vector, returns a one-row matrix, each column containing the
#' probability of having inherited that parental homologue (the first columns corresponding
#' to parent 1, and the rest to parent 2). If `geno` is a matrix, a list of matrices is returned.
#' Each element is the probabilities for all individuals of having inherited one parental
#' homologue. There should be ploidy*2 elements in the list.
#' @export
#'
#' @examples
#' data("genotype")
#' geno <- geno[,-1:-2] #we take out parental genotypes
#' data("homologue")
#'
#' #For a single genotype
#' IBD <- calc_IBD(geno[1,1],hom[1,1:2],hom[1,3:4], ploidy = 2)
#' #For a single individual
#' #The genotype must be a matrix of one column!
#' IBD <- calc_IBD(geno[,1,drop = FALSE],hom[,1:2],hom[,3:4], ploidy = 2)
#' #For all individuals
#' IBD <- calc_IBD(geno,hom[,1:2],hom[,3:4], ploidy = 2)
calc_IBD <- function(geno,p1hom,p2hom,ploidy = 2){
  geno_p1 <- identical(class(geno),class(p1hom))
  geno_p2 <- identical(class(geno),class(p2hom))
  if(!(geno_p1 & geno_p2)){
    stop("All data provided to IBD must be of the same class (i.e. all numeric or all matrices)")
  }else{
    UseMethod("calc_IBD",geno)
  }
}



#' Calculation of IBD for a single data point
#'
#' @describeIn calc_IBD Function for a single marker of a single individual.
calc_IBD.numeric <- function(geno,p1hom,p2hom,ploidy = 2){
  if(all(is.na(geno))) return(NA)
  assertthat::assert_that(all(geno <= ploidy))
  assertthat::assert_that(length(p1hom) == ploidy)
  assertthat::assert_that(length(p2hom) == ploidy)

  configlist <- config(geno,ploidy = ploidy)
  parenthom <- list(p1 = p1hom,
                    p2 = p2hom)
  parentcols <- list(p1 = 1:(ploidy/2),
                     p2 = (ploidy/2+1):ploidy)
  #This expresses the maximum and minimum allele genoage
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

#' Calculation of IBD for multiple markers and/or individuals
#'
#' @describeIn calc_IBD Function for handling matrix input to calc_IBD
calc_IBD.matrix <- function(geno,p1hom,p2hom,ploidy = 2){
  if(is.vector(geno)){
    nams <- names(geno)
    geno <- matrix(geno,ncol = 1)
    rownames(geno) <- nams
  }
  if(is.null(rownames(geno))) stop("genoage matrix should contain rownames")
  if(is.null(rownames(p1hom))) stop("P1 homologue matrix should contain rownames")
  if(is.null(rownames(p2hom))) stop("P2 homologue matrix should contain rownames")
  markers_p1 <- rownames(geno) %in% rownames(p1hom)
  if(!all(markers_p1)) stop("Some rownames in genoage matrix were not found in P1 homologue matrix: ",
                            paste(rownames(geno)[!markers_p1],collapse = " "))
  markers_p2 <- rownames(geno) %in% rownames(p2hom)
  if(!all(markers_p2)) stop("Some rownames in genoage matrix were not found in P2 homologue matrix:\n",
                            paste(rownames(geno)[!markers_p2],collapse = " "))
  assertthat::assert_that(ncol(p1hom) == ploidy)
  assertthat::assert_that(ncol(p2hom) == ploidy)

  #To save time, we only calculate each unique case
  #To identify unique cases we transform the data into one long vector
  #This is memory expensive, might need to rework for big individual x genotype matrices
  nmark <- nrow(geno)
  nind <- ncol(geno)
  genovec <- as.vector(geno)
  names(genovec) <- rep(rownames(geno),nind)
  p1hom <- p1hom[names(genovec),]
  p2hom <- p2hom[names(genovec),]

  #First identify the number of combinations
  data <- data.frame(geno = genovec, p1hom, p2hom)
  parentcols <- list(p1 = 2:(ncol(p1hom)+1),
                     p2 = (ncol(p1hom)+2):ncol(data))
  comb <- apply(data,1,paste,collapse = "_")
  marks <- 1:nrow(data)
  cases <- split(marks,comb)

  #Then we fill in a result matrix
  #we only calculate one time each combination
  res <- matrix(0,ncol = ncol(data) -1,nrow = nrow(data))
  rownames(res) <- names(genovec)
  for(i in 1:length(cases)){
    cas <- names(cases)[i]
    dat <- suppressWarnings(as.numeric(strsplit(cas,"_")[[1]]))
    prob <- calc_IBD.numeric(geno = dat[1],p1hom = dat[parentcols$p1],
                             p2hom = dat[parentcols$p2],ploidy = ploidy)
    res[cases[[i]],] <- rep(prob,each = length(cases[[i]]))
  }

  #We then reform the output into one homologue list.
  homlist <- lapply(1:ncol(res),function(i){
    x <- res[,i]
    res <- matrix(x,ncol = nind, nrow= nmark)
    rownames(res) <- rownames(geno)
    colnames(res) <- colnames(geno)
    return(res)
  })
  names(homlist) <- c(colnames(p1hom),colnames(p2hom))
  return(homlist)
}


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
#'
#' @example
#' config(1,4)
config <- function(g,ploidy = 2){
  genoes <- c(rep(1,g),rep(0,ploidy-g))
  return(unique(combinat::permn(genoes)))
}
