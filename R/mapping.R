#Mapping functions ---------
#These functions are mostly made to simplify the process of using polyploid mapping
#procedures defined in polymapR and MDSmap

#This function allows to calculate a list of arguments to pass on
#to the linkage function in order to estimate meaningful linkages between
#all marker types found in the dataset. This allows to
#create a signle data frame containing all linkages


#' Linkage parameter calculator
#'
#' This function aids in the definition of parameters to pass to the
#' `linkage` function in the package `polymapR`.
#'
#' @param genop1 Numeric vector containing genotypes of parent 1
#' @param genop2 Numeric vector containing genotypes of parent 2
#' @param pnames Character vector of length 2 containing parent names
#'
#' @return A list where each element contains the parameters `markertype1`,
#' `markertype2`, `target_parent` and `other_parent` for each type of marker
#' defined by the combination of P1 and P2 genotypes.
#' @keywords internal
#'
#' @examples
linkpar <- function(genop1,genop2,pnames){
  pgeno <- cbind(genop1,genop2)
  genoclass <- unique(apply(pgeno,1,paste,collapse = "_"))
  genoclass <- strsplit(genoclass,"_")
  pgeno <- cbind(as.numeric(extract(genoclass,1)),
                 as.numeric(extract(genoclass,2)))

  self_links <- apply(pgeno,1,function(mark){
    target <- pnames[which.max(mark)]
    other <- pnames[!pnames %in% target]
    type1 <- sort(mark,decreasing = T)
    type2 <- NULL
    return(list(target = target,other = other,markertype1 = type1,markertype2 = type2))
  })

  mark_combs <- combn(1:nrow(pgeno),2)
  other_links <- apply(mark_combs,2,function(index){
    m1index <- which.min(rowSums(pgeno[index,]))
    m2index <- which.max(rowSums(pgeno[index,]))
    m1 <- pgeno[index[m1index],]
    m2 <- pgeno[index[m2index],]
    target <- pnames[as.logical(m1)]
    other <- pnames[!as.logical(m1)]
    m1 <- sort(m1,decreasing = T)
    m2 <- sort(m2,decreasing = T)

    list(target = target,
         other = other,
         markertype1 = m1,
         markertype2 = m2)
  })


  all_events <- c(self_links,other_links)
  type1 <- extract(all_events,"markertype1",simplify = F)
  type2 <- extract(all_events,"markertype2",simplify = F)
  meaningful <- sapply(1:length(type1),function(i){
    if(is.null(type2[[i]])){
      if(sum(type1[[i]])!= 0) return(T)
      else return(F)
    }else{
      t2 <- t(rev(type2[[i]]))
    }
    as.vector(t2 %*% type1[[i]]) != 0
  })
  all_events <- all_events[meaningful]
  all_events_code <- lapply(all_events,function(r) paste(unlist(r),collapse = "_"))
  all_events <- all_events[!duplicated(all_events_code)]
  res <- lapply(names(all_events[[1]]),function(nam){
    extract(all_events,nam,simplify = F)
  })
  names(res) <- names(all_events[[1]])
  return(res)
}

#This is a shortcut for linkage estimation using polymapR
#' Linkage estimation
#'
#' The `polymapR::linkage()` allows to calculate pair-wise linkage estimates between markers
#' using the `MDSmap` package. However, it requires to be run once for each marker pair
#' combination (for 1x0 - 1x0, 1x0 - 1x1, 1x1 - 0x1 ...). To simplify this process one can
#' use this shortcut, which analyzes the present genotypes and launches `linkage()` as many
#' times as necessary to obtain all meaningful pair-wise estimates (for instance, 1x0 - 0x1
#' pairs are ignored since they are always independent). It returns a data-frame
#' with all pair-wise links.
#'
#' @param geno Genotype matrix with individuals on columns and markers on rows. Colum names
#' are expected.
#' @param p1name Column name of the first parent.
#' @param p2name Column name of the second parent.
#' @param ploidy Ploidy of the population (only even ploidies identical for the whole population
#' are allowed)
#' @param ncores Number of cores to use in the linkage calculation. It is recommended to
#' set this number at most to the number of available cores - 1.
#'
#' @return Data.frame containing recombination frequency, LOD estimate and LOD for independence
#' for each pair of markers. See `polymapR::linkage()` for more information.
#' @keywords internal
#'
#' @examples
linkdf_shortcut <- function(geno,p1name,p2name,ploidy,ncores = 1){
  linkparlist <- linkpar(geno[,p1name],geno[,p2name],c(p1name,p2name))
  links <- lapply(1:length(linkparlist[[1]]),function(i){
    suppressWarnings(
      suppressMessages(
        polymapR::linkage(dosage_matrix = as.matrix(geno),
                          markertype1 = linkparlist$markertype1[[i]],
                          markertype2 = linkparlist$markertype2[[i]],
                          target_parent = linkparlist$target[[i]],
                          other_parent = linkparlist$other[[i]],
                          ploidy = 2,G2_test = T,ncores = ncores)
      )
    )
  })
  links <- do.call(rbind,links)
  return(links)
}

#' MDSmap caller
#'
#' Shortcut function to launch MDSmap
#'
#' @param linkdf Linkage data.frame as defined in the `MDSMap::estimate.map()` function.
#' @param ndim 2 or 3, number of dimensions to use in MDS-mapping.
#'
#' @return Estimated linkage map, including position and other parameters.
#' @keywords internal
#'
#' @examples
mdsmap <- function(linkdf,ndim = 2){
  linkdf <- linkdf[order(linkdf$marker_a,linkdf$marker_b),]
  text <- paste(c(length(unique(unlist(linkdf[,1:2]))),
                  apply(linkdf[,c(1,2,3,4)],1,paste,collapse = " ")),
                collapse = "\n")
  write(x = text,file = "tmp.map")
  if(ndim == 2){
    newmap <- MDSMap::calc.maps.pc("tmp.map",weightfn = "lod2")
  }else if(ndim == 3){
    newmap <- custom_calc.maps.sphere("tmp.map",weightfn = "lod2",p= 200)
  }
  rem <- file.remove("tmp.map")
  return(newmap)
}


#' 3-dimensional MDS-mapping
#'
#' Custom function identical to `MDSMap::calc.maps.sphere()` except for a reduction
#' in the number of iterations used from 1e06 to 2e04, which greately reduces computational
#' time with little change in the final outcome.
#'
#' @param fname Character string specifying the base name of the file fname.txt
#' which contains the data to be analysed.
#' The file should be white space or tab separated.
#' @param p Integer - the penalty for deviations from the sphere - higher p
#' forces points more closely onto a sphere.
#' @param n Vector of integers or strings containing markers to be omitted from the analysis.
#' @param weightfn Character string specifying the values to use for
#' the weight matrix in the MDS 'lod2' or 'lod'.
#' @param mapfn Character string specifying the map function to use on the recombination
#' fractions 'haldane' is default, 'kosambi' or 'none'.
#'
#' @return Estimated linkage map. See `MDSMap::calc.maps.sphere()` for further details.
#' @keywords internal
#'
#' @examples
custom_calc.maps.sphere <- function (fname, p = 100, n = NULL, weightfn = "lod2", mapfn = "haldane") {
  lodrf <- MDSMap::calc.pair.rf.lod(fname, weightfn)
  confplotno <- 1:lodrf$nloci
  if (!is.null(n)) {
    if (!is.numeric(n))
      n <- which(lodrf$locinames %in% n)
    r <- lodrf$rf[-n, -n]
    lod <- lodrf$lod[-n, -n]
    confplotno <- confplotno[-n]
  }
  else {
    r <- lodrf$rf
    lod <- lodrf$lod
  }
  M <- MDSMap::dmap(r, mapfn)
  nloci = length(confplotno)
  smacofsym <- smacof::smacofSym(M, ndim = 2, weightmat = lod,
                                 itmax = 1e+04)
  smacofsphere <- smacof::smacofSphere(M, ndim = 2, algorithm = "dual",
                                       weightmat = lod,
                                       penalty = p,
                                       itmax = 2e+04, mod = 10,
                                       verbose = FALSE)
  mapsphere <- MDSMap::map.to.interval(smacofsphere, nloci)
  length <- mapsphere$chromlength[nloci]
  distmap <- outer(mapsphere$maporder, mapsphere$maporder,
                   Vectorize(function(i, j) M[i, j]))
  lodmap <- outer(mapsphere$maporder, mapsphere$maporder,
                  Vectorize(function(i, j) lod[i, j]))
  if (!is.null(n)) {
    locikey <- data.frame(locus = lodrf$locinames[-n], confplotno = confplotno)
  }
  else {
    locikey <- data.frame(locus = lodrf$locinames, confplotno = confplotno)
  }
  sr = smacofsphere$stress/smacofsym$stress
  ssphere = smacofsphere$stress
  ssym = smacofsym$stress
  nnfit <- MDSMap::calc.nnfit(distmap, lodmap, mapsphere$chromlength)
  locimap <- data.frame(confplotno = confplotno[mapsphere$maporder],
                        locus = locikey$locus[mapsphere$maporder], position = mapsphere$chromlength,
                        nnfit = nnfit$pointfits, row.names = 1:nloci)
  if (!is.null(n)) {
    removedloci <- data.frame(n, lodrf$locinames[n], row.names = NULL)
  }
  else {
    removedloci <- n
  }
  retlist <- list(smacofsym = smacofsym, smacofsphere = smacofsphere,
                  mapsphere = mapsphere, distmap = distmap, lodmap = lodmap,
                  locimap = locimap, length = length, removed = n, locikey = locikey,
                  stressratio = sr, ssphere = ssphere, ssym = ssym, meannnfit = nnfit$meanfit)
  class(retlist) <- "spheremap"
  retlist
}

