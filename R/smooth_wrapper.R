#Wrappers
#These functions streamline the smooth descent algorithm: from IBD prediction to error estimation
#and genotype calculation.

#Two flavours exist: only smoothing, or smoothing and remapping (using MDSMap algorithm).

#Smooth Descent -------------
#' Smooth Descent wrapper
#'
#' This function applies the IBD calculation, IBD prediction, error estimation
#' and genotype prediction functions according to the Smooth Descent algorithm.
#' Provided with a genotype matrix, a parental homologue assignment matrix and
#' a genetic map, it is able to estimate putative genotyping errors as well
#' as recombination counts and imputed genotypes with less errors.
#'
#' @param geno matrix with markers on the rows, individuals on the columns. Row names
#' are expected.
#' @param homologue matrix with markers on the rows, homologue names on the columns.
#' Rownames and columnames expected.
#' @param map data.frame with at least columns "marker" and "position". If it is
#' not specified, a map will be estimated from the uncorrected genotype data with polymapR.
#' @param ploidy numeric indicating the ploidy. Both parents must be of the same ploidy,
#' and it is assumed that "homologue" has 2*ploidy columns.
#' @param p1name character, name of the first parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the first column.
#' @param p2name character, name of the second parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the second column.
#' @param prediction_interval numeric, interval to be used during the IBD prediction
#' step. It should be specified in the same units as the "position" in the map.
#' @param error_threshold numeric, threshold over which a marker is considered erroneous.
#' Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)
#' @param verbose logical, should smooth descent report the steps it takes?
#'
#' @return A list containing the following items:
#' * obsIBD: list of observed IBD matrices (marker x individual) for each parental homologue
#' * predIBD: list of predicted IBD matrices (marker x individual) for each parental homologue
#' * oldmap: data.frame containing the original map
#' * error: list of error matrices (marker x individual) for each parental homologue
#' * newIBD: list of observed IBD matrices of the corrected genotypes.
#' * newgeno: matrix containing the new genotypes.
#' * rec: list containing the recombination counts using the observed IBD (`obs`),
#' the predicted IBD (`pred`) and the updated IBD (`new`), using the given map order.
#' For more information see `rec_count()`.
#' @export
#'
#' @examples
smooth_descent <- function(geno,
                           homologue,
                           map,
                           ploidy = 2,
                           p1name = NULL,
                           p2name = NULL,
                           prediction_interval = 10,
                           error_threshold = 0.8,
                           verbose = T){
  talk <- function(msg){
    if(verbose) cat(msg)
  }
  ## INPUT CHECK ##
  assertthat::assert_that(!is.null(rownames(geno)))
  assertthat::assert_that(!is.null(rownames(homologue)))
  if(!is.null(map)){
    assertthat::assert_that(assertthat::has_name(map,"marker"))
    assertthat::assert_that(assertthat::has_name(map,"position"))
  }

  #We eliminate non-assigned homologues
  homologue <- homologue[!is.na(rowSums(homologue)),]

  #Then we make sure that geno, homologue and map
  #have the same markers and are in the same order
  #First the matrices with rownames
  shared_marks <- intersect(rownames(geno),rownames(homologue))
  if(length(shared_marks) == 0)    stop("No shared markers found between genotype and homologue matrices")
  geno <- geno[shared_marks,]
  homologue <- homologue[shared_marks,]

  #Map and genotypes need to be synchronized
  mapped_marks <- intersect(map$marker,rownames(geno))
  if(length(mapped_marks) == 0) stop("No shared markers found between genotype matrix and map data.frame")
  map <- subset(map,marker %in% mapped_marks)
  map <- map[order(map$position),]
  geno <- geno[map$marker,]
  homologue <- homologue[map$marker,]


  #Error estimation
  parentcols <- list(p1 = 1:ploidy,
                     p2 = (ploidy+1):ncol(homologue))
  genocols <- which(colnames(geno) %in% c(p1name,p2name))

  talk("Obtaining IBD\n")
  #The parents should not be included in the IBD calculation
  full_obsIBD <- calc_IBD(dos = as.matrix(geno[,-genocols]),
                          p1hom = as.matrix(homologue[,parentcols$p1]),
                          p2hom = as.matrix(homologue[,parentcols$p2]),
                          ploidy = ploidy)
  obsIBD <- full_obsIBD
  # #With this I eliminate from the observed IBD the non-informative markers
  # obsIBD <- lapply(full_obsIBD,function(obs){
  #   #All markers for which all IBDs are >0.3 and <0.7 are considered
  #   #uninformative. They are not used for calculating predicted IBDs
  #   #and do not influence genotype estimation
  #   uninf <- obs > 0.3 & obs < 0.7
  #   inf_rows <- (rowSums(uninf | is.na(obs))) != ncol(obs)
  #   obs[inf_rows,]
  # })
  talk("Predicting IBD\n")
  predIBD <- predict_IBD(full_obsIBD, map, interval = prediction_interval)
  talk("Detecting errors\n")
  errors <- lapply(1:length(obsIBD),function(i){
    err <- abs(obsIBD[[i]] - predIBD[[i]]) > error_threshold
    err | is.na(obsIBD[[i]])
  })
  tots <- sapply(errors,sum)

  if(sum(tots) == 0){
    warning("No errors detected, no smoothing performed")

    res <- list(obsIBD = full_obsIBD,
                predIBD = predIBD,
                oldmap = map,
                error = errors,
                newIBD = full_obsIBD,
                oldgeno = geno,
                newgeno = geno)

  }else{
    newIBD <- lapply(1:length(obsIBD),function(i){
      ibd <- obsIBD[[i]]
      if(length(ibd) == 0) return(ibd)
      ibd[errors[[i]]] <- round(predIBD[[i]][errors[[i]]])
      return(ibd)
    })
    names(newIBD) <- names(obsIBD)

    talk("Updating genotypes\n")
    new_geno <- genotype(newIBD,
                         homologue = homologue,
                         threshold = 0.8,
                         ploidy = ploidy)
    new_geno <- new_geno[rownames(geno),colnames(geno[,-genocols])]
    new_geno[is.na(new_geno)] <- geno[,-genocols][is.na(new_geno)]
    new_geno <- cbind(geno[,genocols],new_geno)
    newIBD <- calc_IBD(dos = as.matrix(new_geno[,-genocols]),
                       p1hom = as.matrix(homologue[,parentcols$p1]),
                       p2hom = as.matrix(homologue[,parentcols$p2]))

    res <- list(obsIBD = full_obsIBD,
                predIBD = predIBD,
                oldmap = map,
                error = errors,
                newIBD = newIBD,
                oldgeno = geno,
                newgeno = new_geno)
  }

  #Diagnostics
  #Because here we're just recalculating genotypes, not
  #re-mapping, only the recombination calculation can be made
  talk("Counting recombinations\n")
  obs_rec <- rec_count(res$obsIBD,map)
  pred_rec <- rec_count(res$predIBD,map)
  new_rec <- rec_count(res$newIBD,map)

  res$rec <- list(obs = obs_rec,
                  pred = pred_rec,
                  new = new_rec)

  res <- res[c("obsIBD","predIBD","error","newIBD","oldgeno","newgeno","oldmap","rec")]
  return(res)
}

#Smooth Map --------------

#' Smooth Mapping
#'
#' This function allows to apply one single interation of the Smooth Descent
#' algorithm, including a re-mapping step using new genotypes. Optionally, the
#' preliminary map can also be estimated instead of being provided. Both mapping
#' procedures are performed with packages `polymapR` and `MDSMap`, which perform
#' multi-dimensional scaling mapping. For usage of the Smooth Descent algorithm
#' with a different map algorithm see `smooth_descent()`, which performs the
#' genotype correction without re-mapping.
#'
#' @param geno matrix with markers on the rows, individuals on the columns.
#' Rownames and columnames expected.
#' @param homologue matrix with markers on the rows, homologue names on the columns.
#' Rownames and columnames expected.
#' @param map optionally, data.frame with at least columns "marker" and "position". If it is
#' not specified, a map will be estimated from the uncorrected genotype data with polymapR.
#' @param ploidy numeric indicating the ploidy. Both parents must be of the same ploidy,
#' and it is assumed that "homologue" has 2*ploidy columns.
#' @param prediction_interval numeric, interval to be used during the IBD prediction
#' step. It should be specified in the same units as the "position" in the map.
#' @param error_threshold numeric, threshold over which a marker is considered erroneous.
#' Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)
#' @param ncores number of cores to use for linkage estimation.
#' @param mapping_ndim 2 or 3, number of dimensions to use for multi-dimensional mapping
#' @param p1name character, name of the first parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the first column.
#' @param p2name character, name of the second parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the second column.
#' @param estimate_premap logical, whether to use polymapR to estimate a preliminary map.
#' @param verbose logical, should smooth descent report the steps it takes?
#'
#' @return list containing the following items:
#' * obsIBD: list of observed IBD matrices (marker x individual) for each parental homologue
#' * predIBD: list of predicted IBD matrices (marker x individual) for each parental homologue
#' * oldmap: data.frame containing the original map
#' * error: list of error matrices (marker x individual) for each parental homologue
#' * newIBD: list of observed IBD matrices of the corrected genotypes.
#' * newmap: data.frame containing the new updated map.
#' * newgeno: matrix containing the new genotypes.
#' * rec: list containing the recombination counts using the observed IBD (`obs`),
#' the predicted IBD (`pred`) and the updated IBD (`new`). `obs` and `pred` use the old
#' map order while `new` uses the re-estimated map. For more information see `rec_count()`.
#' * recdist: data.frame containing pair-wise recombination and distance between markers,
#' useful to plot using `recdist_plot()`
#' * tau: reordering parameter tau (Kendall's rank correlation). Obtained with `reorder_tau()`
#' * eliminated: markers eliminated due to a large average neighbour distance, they tend
#' to be problematic when re-mapping.
#'
#' @export
#'
#' @examples
smooth_map <- function(geno,
                      homologue,
                      map = NULL,
                      ploidy = 2,
                      p1name = NULL,
                      p2name = NULL,
                      prediction_interval = 10,
                      error_threshold = 0.8,
                      ncores = 10,
                      mapping_ndim = 2,
                      estimate_premap = F,
                      verbose = T){
  talk <- function(msg){
    if(verbose) cat(msg)
  }

  #In case the preliminary map has not been given, or
  #a new map should be estimated, we calculate a new one
  if(is.null(map) | estimate_premap){
    talk("Estimating preliminary map\n")
    if(is.null(p1name)) p1name <- colnames(geno)[1]
    if(is.null(p2name)) p2name <- colnames(geno)[2]
    linkdf <- linkdf_shortcut(geno,ncores,p1name = p1name, p2name = p2name)
    map <- mdsmap(linkdf,mapping_ndim)
    map <- map$locimap
    colnames(map)[2] <- "marker"
  }

  #We should eliminate markers that are distant from everything
  #For the remapping process, distant markers should be eliminated
  closeness <- sapply(map$position,function(p){
    mean(sort(abs(p - map$position))[1:5])
  })
  far_marks <- closeness >= 3
  if(any(far_marks)){
    warning(sum(far_marks)," markers were eliminated from the map due to large neighbour distance: ",
            paste(map$marker[far_marks],collapse = " "))
  }
  map <- map[!far_marks,]

  res <- smooth_descent(geno = geno, homologue = homologue, map = map,
                 ploidy = ploidy, p1name = p1name, p2name = p2name,
                 prediction_interval = prediction_interval, error_threshold = error_threshold,
                 verbose = verbose)

  talk("Estimating new linkage\n")
  linkdf <- linkdf_shortcut(geno = res$newgeno,ploidy = 2,ncores = ncores,
                            p1name = p1name, p2name = p2name)


  talk("Mapping new genotypes\n")
  new_map <- mdsmap(linkdf,mapping_ndim)
  new_map <- new_map$locimap
  colnames(new_map)[2] <- "marker"

  #Diagnostics and output formatting
  res$oldmap = map
  res$newmap = new_map
  res$eliminated = map$marker[far_marks]
  res$recdist <- recdist_calc(res$newmap,linkdf)
  res$tau <- reorder_tau(res$oldmap,res$newmap)
  res$rec$new <- rec_count(res$newIBD,map = res$newmap)

  return(res)
}


#Test --------
for(file in list.files("R/",full.names = T)){ source(file)}
geno <- read.table("test/test_geno.txt",header = T)
map <- read.table("test/test_map.txt",header = T)
hom <- read.table("test/test_hom.txt",header = T)

smooth_test <- smooth_descent(geno,map,homologue = hom,ploidy = 2,
               p1name = "P1",p2name = "P2",verbose = T)

smoothmap_test <- smooth_map(geno,map,homologue = hom,ploidy = 2,
                             p1name = "P1",p2name = "P2",verbose = T,
                             mapping_ndim = 3)


iterplot(smoothmap_test[c("oldmap","newmap")])
