#IBD calculation -----------------
calc_IBD <- function(geno,p1hom,p2hom,map = NULL,ploidy = 2,method = "naive",
                     p1name = "P1",p2name = "P2"){

  method <- match.arg(method, c("naive","hmm","heuristic"))
  geno_p1 <- identical(class(geno),class(p1hom))
  geno_p2 <- identical(class(geno),class(p2hom))
  if(!(geno_p1 & geno_p2)){
    stop("All data provided to IBD must be of the same class (i.e. all numeric or all matrices)")
  }else{
    switch(method,
           naive = UseMethod("naive_IBD",geno),
           hmm = UseMethod("hmm_IBD",geno),
           heuristic = UseMethod("heur_IBD",geno))
  }
}


### SD Implementation ---------

#' Calculation of IBD
#'
#' This function computes Identity-by-descent probabilities according to a naive model
#' based on random assortment of markers. It requires homologue matrices (`p1hom` and `p2hom`)
#' indicating which homologues of each parent carry a genotype of 1 or 0 for the tracked allele.
#' Based on this, a probability matrix is returned: for each individual how likely it is,
#' at each marker, that it has inherited each of the homologues.
#' As a generic method, it can use a single genotype of one marker (and two vectors
#' of homologue assignment); a vector of markers (and two homologue matrices) or a matrix of markers
#' (and two homologue matrices).
#'
#' @param geno integer or matrix containing genotypes (genotypes of allele A, from 0 to ploidy).
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
#' offspring genotype is impossible given the parental homologues, NA is returned.
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
#' IBD <- naive_IBD(geno[1,1],hom[1,1:2],hom[1,3:4], ploidy = 2)
#' #For a single individual
#' #The genotype must be a matrix of one column!
#' IBD <- naive_IBD(geno[,1,drop = FALSE],hom[,1:2],hom[,3:4], ploidy = 2)
#' #For all individuals
#' IBD <- naive_IBD(geno,hom[,1:2],hom[,3:4], ploidy = 2)
naive_IBD <- function(geno,p1hom,p2hom,ploidy = 2,p1name = "P1",p2name = "P2",...){
  geno_p1 <- identical(class(geno),class(p1hom))
  geno_p2 <- identical(class(geno),class(p2hom))
  if(!(geno_p1 & geno_p2)){
    stop("All data provided to IBD must be of the same class (i.e. all numeric or all matrices)")
  }else{
    UseMethod("naive_IBD",geno)
  }
}

#' Calculation of IBD for a single data point
#'
#' @describeIn naive_IBD Function for a single marker of a single individual.
naive_IBD.numeric <- function(geno,p1hom,p2hom,ploidy = 2,...){
  if(all(is.na(geno))) return(NA)
  assertthat::assert_that(all(geno <= ploidy))
  assertthat::assert_that(length(p1hom) == ploidy)
  assertthat::assert_that(length(p2hom) == ploidy)

  configlist <- config(geno,ploidy = ploidy)
  parenthom <- list(p1 = p1hom,
                    p2 = p2hom)
  parentcols <- list(p1 = 1:(ploidy/2),
                     p2 = (ploidy/2+1):ploidy)
  #This expresses the maximum and minimum allele genotype
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
#' @describeIn naive_IBD Function for handling matrix input to naive_IBD
naive_IBD.matrix <- function(geno,p1hom,p2hom,ploidy = 2, p1name = "P1",p2name = "P2",...){
  if(is.vector(geno)){
    nams <- names(geno)
    geno <- matrix(geno,ncol = 1)
    rownames(geno) <- nams
  }
  if(is.null(rownames(geno))) stop("genotype matrix should contain rownames")
  if(is.null(rownames(p1hom))) stop("P1 homologue matrix should contain rownames")
  if(is.null(rownames(p2hom))) stop("P2 homologue matrix should contain rownames")
  markers_p1 <- rownames(geno) %in% rownames(p1hom)
  if(!all(markers_p1)) stop("Some rownames in genotype matrix were not found in P1 homologue matrix: ",
                            paste(rownames(geno)[!markers_p1],collapse = " "))
  markers_p2 <- rownames(geno) %in% rownames(p2hom)
  if(!all(markers_p2)) stop("Some rownames in genotype matrix were not found in P2 homologue matrix:\n",
                            paste(rownames(geno)[!markers_p2],collapse = " "))
  assertthat::assert_that(ncol(p1hom) == ploidy)
  assertthat::assert_that(ncol(p2hom) == ploidy)
  pindex <- which(colnames(geno) %in% c(p1name,p2name))
  if(length(pindex) > 0)  geno <- geno[,-pindex,drop=F]

  #To save time, we only calculate each unique case
  #To identify unique cases we transform the data into one long vector
  #This is memory expensive, might need to rework for big individual x genotype matrices
  #Alterative method to identify cases
  homcases <- apply(cbind(p1hom,p2hom),1,paste0,collapse = "_")
  unhom <- unique(homcases)
  ung <- unique(as.vector(geno))

  whichom <- lapply(unhom,function(h) homcases == h)
  whichg <- lapply(ung, function(g) geno == g)

  cases <- lapply(whichg,function(g){
    lapply(whichom, function(h){
      mh <- matrix(rep(h,ncol(g)),nrow = length(h))
      which(as.logical(mh*g))
    })
  })
  cases <- unlist(cases,recursive = F)
  names(cases) <- apply(expand.grid(unhom,ung),1,paste,collapse = "_")

  #Then we fill in a result matrix
  #we only calculate one time each combination
  matfill <- rep(NA, (ncol(p1hom) + ncol(p2hom)))
  res <- lapply(matfill,matrix,
                ncol = ncol(geno),
                nrow = nrow(geno))
  parentcols <- list(p1 = 1:ncol(p1hom),
                     p2 = ncol(p1hom) + 1:ncol(p2hom))
  for(i in 1:length(cases)){
    cas <- names(cases)[i]
    dat <- suppressWarnings(as.numeric(strsplit(cas,"_")[[1]]))

    prob <- naive_IBD.numeric(geno = dat[ncol(p1hom) + ncol(p2hom) + 1],
                             p1hom = dat[parentcols$p1],
                             p2hom = dat[parentcols$p2],ploidy = ploidy)
    for(j in 1:length(res)){
      res[[j]][cases[[i]]] <- prob[j]
    }
  }
  for(i in 1:length(res)){
    dimnames(res[[i]]) <- dimnames(geno)
  }


  names(res) <- c(colnames(p1hom),colnames(p2hom))
  return(res)
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

### HMM implementation ----------
#Based on polyqtlR's hmm model
hmm_IBD <- function(geno,p1hom,p2hom,map,ploidy = 2,
                    p1name = "P1", p2name = "P2", ...){
  if(!c("marker","position") %in% colnames(map)) stop("Map must be a data.frame with marker and position columns")
  UseMethod("hmm_IBD",geno)
}


hmm_IBD.matrix <- function(geno,p1hom,p2hom,map,ploidy = 2,
                           p1name = "P1", p2name = "P2",...){

  homnames <- c(colnames(p1hom),colnames(p2hom))
  phased_map <- list(LG1 = data.frame(marker = map$marker,
                                      position = map$position,
                                      p1hom[map$marker,],p2hom[map$marker,]))
  colnames(phased_map[[1]]) <- c("marker","position",paste0("h",1:(ploidy*2)))
  IBD <- polyqtlR::estimate_IBD(phased_maplist = phased_map, genotypes = geno,
                        ploidy = ploidy, method = "hmm",parent1 = p1name,parent2 = p2name)
  homArray <- per_homologue(IBD$LG1$IBDarray)
  IBD <- lapply(1:dim(homArray)[[2]],function(i) homArray[,i,])
  names(IBD) <- homnames

  return(IBD)
}

hmm_IBD.numeric <- function(...){
  stop("Numeric method for hmm is not implemented: at least two markers required")
}

#Using polyqtlR's estimate_IBD
per_homologue <- function(IBDarray){
  #Function expects an IBD array with following dimensions:
  #[locus, genotypeclass, individuals]
  #And translates it to an array with following dimensions:
  #[locus, homologue, individuals]
  #Genotypeclass should be expressed as a combination
  #of letters where each letter represents a parental homologue

  #First create indicator matrix
  genoclass <- dimnames(IBDarray)[[2]]
  homnames <- unique(unlist(strsplit(genoclass,"")))
  indicatorMatrix <- sapply(homnames,function(x) stringr::str_count(genoclass,x) )
  indicatorMatrix <- indicatorMatrix[,sort(colnames(indicatorMatrix))]

  #Then multiply the matrices
  haploIBD <- apply(IBDarray, c(1,3), function(x) x %*% indicatorMatrix)
  haploIBD <- aperm(haploIBD,c(2,1,3))

  dimnames(haploIBD)[[2]] <- homnames
  return(haploIBD)
}

### Heuristic implementation ----------
heur_IBD <- function(geno,p1hom,p2hom,map = NULL,ploidy = 2,
                     p1name = "P1",p2name = "P2",...){
  UseMethod("heur_IBD",geno)
}

heur_IBD.numeric <- function(geno,p1hom,p2hom,ploidy = 2,
                                p1name = "P1",p2name = "P2",...){
  homnames <- c(names(p1hom),names(p2hom))
  map <- data.frame(marker = "m1",position = 1)
  phased_map <- list(LG1 = data.frame(marker = map$marker,
                                      position = map$position,
                                      p1hom,p2hom))
  colnames(phased_map[[1]]) <- c("marker","position",paste0("h",1:(ploidy*2)))
  g <- matrix(c(sum(p1hom),sum(p2hom),geno),nrow = 1)
  rownames(g) <- "m1"
  colnames(g) <- c("P1","P2","geno")

  IBD <- faster_IBD(phased_maplist = phased_map, dosage_matrix = g,
                    ploidy = ploidy, p1name = p1name, p2name = p2name,
                    method.exp = "interval")
  IBD <- IBD$LG1$IBDarray
  names(IBD) <- homnames
  return(IBD)
}

heur_IBD.matrix <- function(geno,p1hom,p2hom,map = NULL,ploidy = 2,
                            p1name = "P1",p2name = "P2",...){
  homnames <- c(colnames(p1hom),colnames(p2hom))
  phased_map <- list(LG1 = data.frame(marker = map$marker,
                                      position = map$position,
                                      p1hom[map$marker,],p2hom[map$marker,]))
  colnames(phased_map[[1]]) <- c("marker","position",paste0("h",1:(ploidy*2)))

  IBD <- faster_IBD(phased_maplist = phased_map, dosage_matrix = geno,
                    ploidy = ploidy, p1name = p1name, p2name = p2name,
                    method.exp = "interval")
  IBD <- lapply(1:dim(IBD$LG1$IBDarray)[[2]],function(i) IBD$LG1$IBDarray[,i,])
  names(IBD) <- homnames
  return(IBD)
}

#' Extremely fast estimation of identity-by-descent (IBD) probabilities.
#' @description The method of "quick-and-dirty" IBD estimation was originally developed by Bourke (2014) for tetraploid data only, and was subsequently
#' generalised by van Geest et al. (2017). Can be useful for a first quick analysis, particularly in large hexaploid datasets. However, the higher accuracy of IBD
#' probabilities generated by \code{\link{hmm_IBD}} makes that function the preferred choice.
#' @param phased_maplist A list of linkage maps calculated by \code{polymapR::create_phased_maplist}
#' @param dosage_matrix An integer matrix with markers in rows and individuals in columns
#' @param map_function The mapping function to calculate recombination frequency based on map distance (haldane or kosambi)
#' @param ploidy Ploidy level of parents or of the first parent
#' @param ploidy2 Ploidy level of the second parent. By default \code{NULL}, if parents have equal ploidy levels.
#' @param fix_threshold The threshold to fix the IBD probabilities while correcting for the sum of probabilities.
#' @param factor_dist Factor to increase or decrease the recombination frequencies as calculated from the map distances.
#' @param ncores Number of cores to use for multi-core processing.
#' @param p1name Name of parent 1 in the column names of \code{dosage_matrix}
#' @param p2name Name of parent 2 in the column names of \code{dosage_matrix}
#' @param method.exp Model to use for probability expansion. Either "simple" or "interval"
#' @return
#' A nested list (with the same length as phased_maplist). Each list element contains the following items:
#' \item{IBDtype}{Always "haplotypeIBD" for the output of this function}
#' \item{IBDarray}{An array of IBD probabilities. The dimensions of the array are: markers, homologues and individuals.}
#' \item{map}{Integrated linkage map positions of markers used in IBD calculation}
#' \item{parental_phase}{The parental marker phasing, coded in 1 and 0's}
#' \item{biv_dec}{\code{NULL}}
#' \item{gap}{The gap size used in IBD interpolation, by default \code{NULL}. See \code{\link{spline_IBD}}}
#' \item{genocodes}{\code{NULL}}
#' \item{pairing}{\code{NULL}}
#' \item{ploidy}{ploidy of parent 1}
#' \item{ploidy2}{ploidy of parent 2}
#' \item{method}{The method used, here "heur" (heuristic)}
#' \item{error}{The error prior used, not relevant here thus \code{NULL}}
#' @references
#' Bourke P.M. (2014) QTL analysis in polyploids: Model testing and power calculations. Wageningen University (MSc thesis)
#' @examples
#' data("phased_maplist.4x","SNP_dosages.4x")
#' IBD_list.4x <- fast_IBD(phased_maplist = phased_maplist.4x,
#'                         dosage_matrix = SNP_dosages.4x,
#'                         ploidy = 4)
faster_IBD <- function(phased_maplist,
                       dosage_matrix,
                       map_function = "haldane",
                       ploidy,
                       ploidy2 = NULL,
                       fix_threshold = 0.1,
                       factor_dist = 1,
                       ncores = 1,
                       p1name = "P1",
                       p2name = "P2",
                       method.exp = "simple"){

  if(is.null(ploidy2)) ploidy2 <- ploidy
  ploidy1 <- ploidy

  map_function <- match.arg(map_function, choices = c("kosambi", "haldane"))

  # check compatibility of dosages and phased maps
  p_check <- function(){
    pcols <- list(p1 = 3:(2+ploidy1) ,
                  p2 = 3:(2+ploidy2) + ploidy1)
    res <- sapply(1:2,function(p){
      sapply(seq(phased_maplist),function(n){
        pname <- switch(p, '1' = p1name, '2' = p2name)
        phom <- phased_maplist[[n]][,pcols[[p]]]
        all(dosage_matrix[rownames(phom),pname] == rowSums(phom,na.rm = T))
      })
    })
    return(all(res))
  }

  if(!p_check())
    stop("Inconsistent input detected.\nIt is likely that unconverted dosages are being combined with a converted phased maplist, or vice versa.")

  if(map_function == "kosambi"){
    mapfun <- function(d, f = 1){
      d <- d/100/f
      0.5*((exp(4*d)-1)/(exp(4*d)+1))
    }
  } else if (map_function == "haldane"){
    mapfun <- function(d, f = 1){
      d <- d/100/f
      0.5*(1-exp(-2*d))
    }
  }

    outlist <- lapply(seq(phased_maplist),function(i){
      ph <- phased_maplist[[i]]
      posdf <- ph[,c("marker", "position"),drop = F]

      ## If marker data is conflicting and co-located, can lead to problems. Jitter positions if necessary
      while(length(unique(posdf$position)) < nrow(posdf)){
        warning("Non-unique marker positions detected. Jittering positions.")
        posdf[duplicated(posdf$position),"position"] <- jitter(posdf[duplicated(posdf$position),"position"],amount=0.1)
        #Map might no longer be in order:
        ph <- ph[order(posdf$position),,drop = F]
        posdf <- posdf[order(posdf$position),,drop = F]
      }

      posvec <- posdf$position
      mark <- as.character(posdf$marker)
      names(posvec) <- mark

      #Same distmat calculation without the need of reshape2
      #or additional objects
      distmat <- abs(outer(posvec,posvec,FUN = "-"))
      k0_mat <- mapfun(distmat, f = factor_dist)

      rownames(ph) <- ph$marker
      ph <- ph[,-which(colnames(ph) %in% c("marker", "position"))]
      ph <- as.matrix(ph)

      scd <- dosage_matrix[rownames(ph),,drop = F]

      # get parental scores
      #Names should not be assumed
      par <- scd[,c(p1name,p2name),drop = F]
      progeny <- scd[,-which(colnames(scd) %in% c(p1name,p2name)),drop = F]

      # make a garray (=genotype probability array) with only 0.5
      garray <- array(data = 0.5, dim = c(nrow(ph), ploidy1 + ploidy2, ncol(scd)-2),
                      dimnames = list(rownames(scd), colnames(ph), colnames(scd)[3:ncol(scd)]))

      # get all markers with fully informative dosage
      P1_hp <- scd[,p1name,drop = F] == 0.5*ploidy1
      P2_hp <- scd[,p2name,drop = F] == 0.5*ploidy2

      # get all maximum scores per marker which are informative in offspring
      sumP <- rowSums(par)

      # fill in garray
      for(m in seq(nrow(scd))){
        # get progeny that is informative - this is the kernel of the original approach:
        # progeny score which equals (or exceeds, DR?) the sum of the parental scores have inherited
        # all the parental alleles, and so we assign all these homologues with a probability of 1
        progS <- colnames(progeny)[progeny[m,] >= sumP[m] & !is.na(progeny[m,])]

        # assign 1 if so.
        garray[m,ph[m,] == 1,progS] <- 1

        # special case where dosage is fully informative (dosage = 0.5*ploidy)
        # only if one of the P1 dosages should be 0.5*ploidy
        # and marker is assigned to 0.5*ploidy homologues
        # apply zeros to all others
        if((P1_hp[m] | P2_hp[m]) & sum(ph[m,]) == sumP[m]){
          if(P1_hp[m])
            garray[m,c(ph[m,1:ploidy1] == 0, rep(FALSE, ploidy2)),progS] <- 0
          if(P2_hp[m])
            garray[m,c(rep(FALSE, ploidy1) ,ph[m,(ploidy1+1):(ploidy1 + ploidy2)] == 0),progS] <- 0
        }
        # if allele is absent in progeny, also informative:
        prog0 <- colnames(progeny)[progeny[m,] == 0 & !is.na(progeny[m,])]
        garray[m,ph[m,] == 1,prog0] <- 0
      }

      # order garray
      garray <- garray[colnames(distmat),,]

      # find nearest informative and change genotype probability accordingly
      if(method.exp == "simple"){
        garray <- simple_expansion(garray,k0_mat,posvec)
      }else{
        garray <- interval_expansion(garray,k0_mat,posvec)
      }

      # call correct_to_sumScore to array:
      for(g in dimnames(garray)[[3]]){
        corr_P1 <- correct_to_sumScore(subg_P = garray[,1:ploidy1,g],
                                       threshold = fix_threshold, ploidy = ploidy1)
        garray[,1:ploidy1,g] <- corr_P1

        corr_P2 <- correct_to_sumScore(subg_P = garray[,(ploidy1 +1):(ploidy1 + ploidy2),g],
                                       threshold = fix_threshold, ploidy = ploidy2)
        garray[,(ploidy1 +1):(ploidy1 + ploidy2),g] <- corr_P2
      }

      return(list("IBDtype" = "haplotypeIBD",
                  "IBDarray" = garray,
                  "map" = cbind("chromosome" = i,
                                posdf),
                  "parental_phase" = ph,
                  "biv_dec" = NULL,
                  "gap" = NULL,
                  "genocodes" = NULL,
                  "pairing" = NULL,
                  "ploidy" = ploidy,
                  "ploidy2" = ploidy2,
                  "method" = "heur",
                  "error" = NULL))
    })

  if(!is.null(names(phased_maplist))){
    names(outlist) <- names(phased_maplist)
  } else{
    warning(paste("phased_maplist had no LG names assigned. Assigning names",
                  paste0("'LG",seq(length(outlist)),"'", collapse = ", "),"to output list."))
    names(outlist) <- paste0("LG",seq(length(outlist)))
  }

  return(outlist)
} #fast_IBD


#This is the reduction approach (HeurB)
interval_expansion <- function(garray,k0_mat,posvec){
  #Expand probabilities
  for(h in dimnames(garray)[[2]]){
    gp <- garray[,h,]
    inf1 <- gp == 1
    inf0 <- gp == 0
    informative <- inf1 | inf0
    intervals <- group_informative(informative)
    case <- unique(unlist(intervals))
    caseindex <- sapply(intervals,function(i) case %in% i)
    rownames(caseindex) <- case

    for(ca in case){
      bounds <- as.numeric(strsplit(ca,":")[[1]])

      #Start and end boundaries are computed differently
      #So we need to change it in case
      if(bounds[1] == -1 & is.infinite(bounds[2])){
        next()
      }

      if(bounds[1] == -1){
        marks <- 1:bounds[2]
        bounds[1] <- bounds[2]
      }else if(is.infinite(bounds[2])){
        marks <- bounds[1]:length(posvec)
        bounds[2] <- bounds[1]
      }else{
        marks <- bounds[1]:bounds[2]
      }

      #All markers that are not boundary are non-informative
      notinf <- marks[!marks %in% bounds]
      near_inf <- (posvec[notinf] > mean(posvec[bounds])) + 1

      #To which individuals should these results be allocated
      locate <- caseindex[ca,]

      newval <- k0_mat[cbind(notinf,bounds[near_inf])]
      #We check the type (repulsion/coupling)
      closest <- gp[bounds[near_inf],locate,drop = F]
      #All results are equal (but if in coupling they must be reversed)
      newval <- matrix(newval,ncol = ncol(closest),nrow = nrow(closest))
      #If in coupling the probability is reversed
      newval[closest == 1] <- 1 - newval[closest == 1]
      gp[notinf,locate] <- newval
    }
    garray[,h,] <- gp
  }
  return(garray)
}
#This is a helper function for HeurB
group_informative <- function(mat){
  #Given a logical matrix, returns a list
  #Where each element contains the intervals of
  #False values in each column. The first
  #Interval starts at -1 and the last ends at Inf
  intervals <- apply(mat,2,function(i){
    m <- which(i)
    if(length(m) == 0) return("-1:Inf")
    move <- diff(m)

    #We define boundaries of non-informative intervals
    m[move != 1]
    starts <- c(m[move != 1])
    ends <- c(m[-1][move != 1])
    #Not sure why but this is necessary for correct interval definition
    if(length(ends) != length(starts)) starts <- starts[-length(starts)]
    intervals <- paste0(starts,":",ends)
    if(all(intervals == ":")) intervals <- NULL
    intervals <- c(paste0(c(-1,m[length(m)]),":",c(m[1],Inf)),intervals)
  })

}

#This is equivalent to the previous approach (HeurA)
simple_expansion <- function(garray,k0_mat,posvec){
  for(g in dimnames(garray)[[3]]){
    for(h in dimnames(garray)[[2]]){
      gp <- garray[,h,g]
      inf1 <- gp == 1
      inf0 <- gp == 0
      informative <- inf1 | inf0

      if(any(informative)){
        mdist <- k0_mat[names(gp)[!informative], names(gp)[informative], drop = FALSE]
        near_inf <- findNearest(posvec[!informative],posvec[informative])
        #The k1 is equal to 1-k0 so this makes sense:
        if(any(inf1)) mdist[,names(inf1)[inf1]] <- 1 - mdist[,names(inf1)[inf1]]
        garray[rownames(mdist),h,g] <- mdist[cbind(seq(near_inf), near_inf)]
      }

    }
  }
  return(garray)
}
#This function simplifies finding the nearest marker
findNearest <- function(a,b){
  span <- range(c(a,b))
  intervals <- c(span[1] - 1,
                 (b[-length(b)] + b[-1])/2,
                 span[2] + 1)
  findInterval(a,intervals)
}

# function to correct total sum of genotype probabilties to 0.5*ploidy
correct_to_sumScore <- function(subg_P, threshold = 0.1, ploidy){
  sumScore <- 0.5*ploidy
  parentalDose <- rowSums(subg_P)

  #Make sure that IBDs don't exceed 0.5*ploidy
  exceed <- parentalDose > sumScore
  subg_P[exceed,] <- subg_P[exceed,]*sumScore/parentalDose[exceed]

  #Fully informative if 1/2 is 0
  inf0 <- rowSums(subg_P == 0) == 0.5*ploidy
  subg_P[inf0,][subg_P[inf0,] != 0] <- 1
  #Fully informative if 1/2 is 1
  inf1 <- rowSums(subg_P == 1) == 0.5*ploidy
  subg_P[inf1,][subg_P[inf1,] != 1] <- 0

  #Those that are not informative are threshold corrected
  #Not sure what is happening here
  fix <- (subg_P < threshold) | (subg_P > (1 - threshold))
  fixable <- rowSums(fix) == ploidy
  fix[fixable,] <- subg_P[fixable,] == 1 | subg_P[fixable,] == 0

  corval <- subg_P
  corval[!fix] <- 0
  corr <- sumScore - rowSums(corval)

  corval <- subg_P
  corval[fix] <- 0
  corval <- corval*corr/rowSums(corval)
  subg_P[!fix] <- corval[!fix]

  # Some probabilites can end up >1. Correct for this:
  #What to do with double reductions?
  if(any(subg_P > 1)){
    #Is over 1 but not close to 2 (1.5 or more is considered close )
    over <- rowSums(subg_P > 1) > 0 & rowSums(round(subg_P) == 2) != 1
    need_cor <- subg_P[over,,drop = F]

    addp <- need_cor < 1 & need_cor > 0 #where to add the probabilities
    remp <- need_cor > 1 #where to remove probabilities

    corr_good <- need_cor
    corr_good[!addp] <- 0 #This object only stores where we need to add
    corr_bad <- need_cor
    corr_bad[!remp] <- 0 #This object only stores where we need to remove

    addval <- (rowSums(corr_bad) - rowSums(remp))/rowSums(addp)
    addval <- rep(addval,times = rowSums(addp))
    need_cor[addp] <- need_cor[addp] + addval

    need_cor[remp] <- 1
    subg_P[over,] <- need_cor
  }

  return(subg_P)
}
