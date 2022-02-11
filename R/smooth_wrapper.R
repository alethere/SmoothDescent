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
#' @param prediction_threshold float, probability threshold for imputing new genotypes.
#' All new genotypes with a probability under this threshold will be considered uncertain.
#' Defaults to 0.8.
#' @param prediction_points numeric, number of points to use for IBD prediction. If NULL, all
#' points in map$position are used, otherwise n equally spaced points are used. Greatly
#' improves efficiency if the number of markers is very large.
#' @param error_threshold numeric, threshold over which a marker is considered erroneous.
#' Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)
#' @param non_inf numeric, lower and upper probability boundaries to consider an
#' IBD probability non-informative (if they fall within the threshold they will be
#' ignored during prediction). Defaults to 0.3 - 0.7. Symmetrical boundaries are
#' recommended but not necessary.
#' @param verbose logical, should smooth descent report the steps it takes?
#'
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
#'
#' data("genotype")
#' data("homologue")
#' data("map")
#'
#' res <- smooth_descent(geno,hom,map, ploidy = 2, p1name = "P1", p2name = "P2")
#'
#'
smooth_descent <- function(geno,
                           homologue,
                           map,
                           ploidy = 2,
                           p1name = NULL,
                           p2name = NULL,
                           prediction_interval = 10,
                           prediction_threshold = 0.8,
                           prediction_points = NULL,
                           error_threshold = 0.8,
                           non_inf = c(0.3,0.7),
                           verbose = T,
                           obs.method = "naive",
                           pred.method = "prediction"){
  talk <- function(...){
    if(verbose) cat(...)
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
  map <- subset(map,map$marker %in% mapped_marks)
  map <- map[order(map$position),]
  geno <- geno[map$marker,]
  homologue <- homologue[map$marker,]

  #Matching parameters
  obs.method <- match.arg(obs.method,c("naive","heuristic"))
  pred.method <- match.arg(pred.method,c("prediction","hmm"))
  talk("Will estimate errors with",obs.method,"-",pred.method,"combination\n")


  ## Error estimation ##
  if(is.null(p1name)) p1name <- colnames(geno)[1]
  if(is.null(p2name)) p2name <- colnames(geno)[2]
  parentcols <- list(p1 = 1:ploidy,
                     p2 = (ploidy+1):ncol(homologue))
  genocols <- which(colnames(geno) %in% c(p1name,p2name))

  talk("Obtaining IBD\n")
  obsIBD <- calc_IBD(geno = as.matrix(geno),
                      p1hom = as.matrix(homologue[,parentcols$p1]),
                      p2hom = as.matrix(homologue[,parentcols$p2]),
                      ploidy = ploidy,map = map, p1name = p1name, p2name = p2name,
                     method = obs.method)

  talk("Predicting IBD\n")
  if(pred.method == "hmm"){
    predIBD <- calc_IBD(geno = as.matrix(geno),
                        p1hom = as.matrix(homologue[,parentcols$p1]),
                        p2hom = as.matrix(homologue[,parentcols$p2]),
                        ploidy = ploidy,map = map, p1name = p1name, p2name = p2name,
                        method = "hmm")
  }else if(pred.method == "prediction"){
    predIBD <- predict_IBD(obsIBD, map, interval = prediction_interval,non_inf = non_inf,
                           pred_points = prediction_points)
  }

  talk("Detecting errors\n")
  errors <- lapply(1:length(obsIBD),function(i){
    err <- abs(obsIBD[[i]] - predIBD[[i]]) > error_threshold
    err[is.na(err)] <- T
    return(err)
  })
  errors <- Reduce('|',errors)
  errors <- lapply(names(obsIBD),function(i) errors)
  names(errors) <- names(obsIBD)
  tots <- sapply(errors,sum,na.rm = T)

  if(sum(tots) == 0){
    warning("No errors detected, no smoothing performed")

    res <- list(obsIBD = obsIBD,
                predIBD = predIBD,
                oldmap = map,
                error = errors,
                newIBD = obsIBD,
                oldgeno = geno,
                newgeno = geno,
                obs.method = obs.method,
                pred.method = pred.method)

  }else{
    newIBD <- lapply(1:length(obsIBD),function(i){
      ibd <- obsIBD[[i]]
      if(length(ibd) == 0) return(ibd)
      err <- errors[[i]]
      ibd[err] <- predIBD[[i]][err]
      return(ibd)
    })
    names(newIBD) <- names(obsIBD)

    talk("Updating genotypes\n")
    new_geno <- genotype(newIBD,
                         homologue = homologue,
                         threshold = prediction_threshold,
                         ploidy = ploidy)

    new_geno <- new_geno[rownames(geno),colnames(geno[,-genocols])]
    new_geno[is.na(new_geno)] <- geno[,-genocols][is.na(new_geno)]
    new_geno <- cbind(geno[,genocols],new_geno)
    newIBD <- calc_IBD(geno = as.matrix(new_geno[,-genocols]),
                       p1hom = as.matrix(homologue[,parentcols$p1]),
                       p2hom = as.matrix(homologue[,parentcols$p2]),
                       ploidy = ploidy)

    res <- list(obsIBD = obsIBD,
                predIBD = predIBD,
                oldmap = map,
                error = errors,
                newIBD = newIBD,
                oldgeno = geno,
                newgeno = new_geno,
                obs.method = obs.method,
                pred.method = pred.method)
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

  res <- res[c("obsIBD","predIBD","error","newIBD","oldgeno","newgeno","oldmap","rec","obs.method","pred.method")]
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
#' @param p1name character, name of the first parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the first column.
#' @param p2name character, name of the second parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the second column.
#' @param prediction_interval numeric, interval to be used during the IBD prediction
#' step. It should be specified in the same units as the "position" in the map.
#' @param prediction_threshold float, probability threshold for imputing new genotypes.
#' All new genotypes with a probability under this threshold will be considered uncertain.
#' Defaults to 0.8.
#' @param prediction_points numeric, number of points to use for IBD prediction. If NULL, all
#' points in map$position are used, otherwise n equally spaced points are used. Greatly
#' improves efficiency if the number of markers is very large.
#' @param error_threshold numeric, threshold over which a marker is considered erroneous.
#' Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)
#' @param ncores number of cores to use for linkage estimation.
#' @param mapping_ndim 2 or 3, number of dimensions to use for multi-dimensional mapping
#' @param estimate_premap logical, whether to use polymapR to estimate a preliminary map.
#' @param max_distance numeric, markers that have near neighbours will be eliminated. This
#' parameter defines the maximum neighbour distance allowed. A warning will be issued if some markers
#' are eliminated.
#' @param non_inf numeric, lower and upper probability boundaries to consider an
#' IBD probability non-informative (if they fall within the threshold they will be
#' ignored during prediction). Defaults to 0.3 - 0.7. Symmetrical boundaries are
#' recommended but not necessary.
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
#' * r2: R-squared parameter of pair-wise recombination and final map distance. A higher
#' value indicates a better newmap.
#' * tau: reordering parameter tau (Kendall's rank correlation). Obtained with `reorder_tau()`
#' * eliminated: markers eliminated due to a large average neighbour distance, they tend
#' to be problematic when re-mapping.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data("genotype")
#' data("homologue")
#' data("map")
#'
#' res <- smooth_descent(geno,hom,map, ploidy = 2, p1name = "P1", p2name = "P2",
#'    estimate_premap = F, mapping_ndim = 3, ncores = 1)
#' }
smooth_map <- function(geno,
                      homologue,
                      map = NULL,
                      ploidy = 2,
                      p1name = NULL,
                      p2name = NULL,
                      prediction_interval = 10,
                      prediction_threshold = 0.8,
                      prediction_points = NULL,
                      error_threshold = 0.8,
                      ncores = 1,
                      mapping_ndim = 2,
                      estimate_premap = F,
                      max_distance = 10,
                      non_inf = c(0.3,0.7),
                      verbose = T,
                      obs.method = "naive",
                      pred.method = "prediction"){
  talk <- function(msg){
    if(verbose) cat(msg)
  }
  if(is.null(p1name)) p1name <- colnames(geno)[1]
  if(is.null(p2name)) p2name <- colnames(geno)[2]

  #In case the preliminary map has not been given, or
  #a new map should be estimated, we calculate a new one
  if(is.null(map) | estimate_premap){
    talk("Estimating preliminary map\n")

    linkdf <- linkdf_shortcut(geno = as.matrix(geno),ploidy = ploidy,
                              p1name = p1name, p2name = p2name,
                              ncores = ncores)
    map <- mdsmap(linkdf,mapping_ndim)
    map <- map$locimap
    colnames(map)[2] <- "marker"
  }

  #We should eliminate markers that are distant from everything
  #For the remapping process, distant markers should be eliminated
  closeness <- sapply(map$position,function(p){
    mean(sort(abs(p - map$position))[1:5])
  })
  far_marks <- closeness >= max_distance
  if(any(far_marks)){
    warning(sum(far_marks)," markers were eliminated from the map due to large neighbour distance: ",
            paste(map$marker[far_marks],collapse = " "),"\n")
  }
  map <- map[!far_marks,]

  res <- smooth_descent(geno = geno, homologue = homologue, map = map,
                 ploidy = ploidy, p1name = p1name, p2name = p2name,
                 prediction_interval = prediction_interval,
                 prediction_threshold = prediction_threshold,
                 prediction_points = prediction_points,
                 error_threshold = error_threshold,
                 verbose = verbose, obs.method = obs.method,
                 pred.method = pred.method)

  talk("Estimating new linkage\n")
  geno <- as.matrix(res$newgeno)
  linkdf <- linkdf_shortcut(geno = geno,ploidy = ploidy,
                            p1name = p1name, p2name = p2name,
                            ncores = ncores)

  talk("Mapping new genotypes\n")
  new_map <- mdsmap(linkdf,mapping_ndim)
  new_map <- new_map$locimap
  colnames(new_map)[2] <- "marker"

  #Diagnostics and output formatting
  res$oldmap = map
  res$newmap = new_map
  res$eliminated = map$marker[far_marks]
  res$recdist <- recdist_calc(res$newmap,linkdf)
  res$r2 <- attr(res$recdist,"r2")
  res$tau <- reorder_tau(res$oldmap,res$newmap)

  return(res)
}

#Iterators ---------------------
#These function iterates smooth_map repeatedly

#' Iterative Smooth Mapping
#'
#' This function allows to apply multiple interations of the Smooth Descent
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
#' @param iters numeric indicating the number of iterations to perform.
#' @param map optionally, data.frame with at least columns "marker" and "position". If it is
#' not specified, a map will be estimated from the uncorrected genotype data with polymapR.
#' @param ploidy numeric indicating the ploidy. Both parents must be of the same ploidy,
#' @param p1name character, name of the first parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the first column.
#' @param p2name character, name of the second parent. Must be present in the geno columnames.
#' If it's not specified it will be taken as the name of the second column.
#' and it is assumed that "homologue" has 2*ploidy columns.
#' @param prediction_interval numeric, interval to be used during the IBD prediction
#' step. It should be specified in the same units as the "position" in the map.
#' @param prediction_threshold float, probability threshold for imputing new genotypes.
#' All new genotypes with a probability under this threshold will be considered uncertain.
#' Defaults to 0.8.
#' @param prediction_points numeric, number of points to use for IBD prediction. If NULL, all
#' points in map$position are used, otherwise n equally spaced points are used. Greatly
#' improves efficiency if the number of markers is very large.
#' @param error_threshold numeric, threshold over which a marker is considered erroneous.
#' Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)
#' @param ncores number of cores to use for linkage estimation.
#' @param mapping_ndim 2 or 3, number of dimensions to use for multi-dimensional mapping
#' @param estimate_premap logical, whether to use polymapR to estimate a preliminary map.
#' @param max_distance numeric, markers that have near neighbours will be eliminated. This
#' parameter defines the maximum neighbour distance allowed. A warning will be issued if some markers
#' are eliminated.
#' @param non_inf numeric, lower and upper probability boundaries to consider an
#' IBD probability non-informative (if they fall within the threshold they will be
#' ignored during prediction). Defaults to 0.3 - 0.7. Symmetrical boundaries are
#' recommended but not necessary.
#' @param verbose logical, should smooth descent report the steps it takes?
#'
#' @return list of lists, where each element contains the following items:
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
#' * r2: R-squared parameter of pair-wise recombination and final map distance. A higher
#' value indicates a better newmap.
#' * tau: reordering parameter tau (Kendall's rank correlation). Obtained with `reorder_tau()`
#' * eliminated: markers eliminated due to a large average neighbour distance, they tend
#' to be problematic when re-mapping.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data("genotype")
#' data("homologue")
#' data("map")
#'
#' res <- smooth_descent(geno,hom,iters = 5, ploidy = 2, p1name = "P1", p2name = "P2",
#'    estimate_premap = F, mapping_ndim = 3, ncores = 1)
#'
#' }
smooth_map_iter <- function(geno,
                            homologue,
                            iters,
                            map = NULL,
                            ploidy = 2,
                            p1name = NULL,
                            p2name = NULL,
                            prediction_interval = 10,
                            prediction_threshold = 0.8,
                            prediction_points = NULL,
                            error_threshold = 0.8,
                            ncores = 1,
                            mapping_ndim = 2,
                            estimate_premap = F,
                            max_distance = 10,
                            non_inf = c(0.3,0.7),
                            verbose = T,
                            obs.method = "naive",
                            pred.method = "prediction"){

  talk <- function(msg){
    if(verbose) cat(msg)
  }

  #Error threshold manager
  if(length(error_threshold) == 2){
    err <- seq(max(error_threshold),min(error_threshold),length.out = iters)
  }else if(length(error_threshold) == iters){
    err <- error_threshold
  }else if(length(error_threshold) == 1){
    err <- rep(error_threshold,iters)
  }else{
    err <- seq(max(error_threshold),min(error_threshold),length.out = iters)
    warning("Error threshold was not length 1, 2 or ",iters,".\n",
            "These thresholds will be used: ",paste(err,collapse = " "))
  }

  res <- list()
  for(i in 1:iters){
    talk(paste0("Iteration ",i," ----------\n"))
    if(i == 1){
      this_iter <- try(smooth_map(geno = geno, homologue = homologue,
                 map = map, ploidy = ploidy, p1name = p1name,
                 p2name = p2name, prediction_interval = prediction_interval,
                 prediction_threshold = prediction_threshold,
                 error_threshold = err[i],
                 ncores = ncores,
                 mapping_ndim = mapping_ndim, estimate_premap = estimate_premap,
                 max_distance = max_distance, verbose = verbose,
                 prediction_points = prediction_points, obs.method = obs.method,
                 pred.method = pred.method))
    }else{
      #In the next iterations, a new starting map and new genotype is used
      if(class(this_iter) == "try-error") next
      m <- this_iter$newmap
      g <- this_iter$newgeno
      this_iter <- try(smooth_map(geno = g,
                              map = m,
                              homologue = homologue,
                              ploidy = ploidy,
                              p1name = p1name,
                              p2name = p2name,
                              prediction_interval = prediction_interval,
                              prediction_threshold = prediction_threshold,
                              prediction_points = prediction_points,
                              error_threshold = err[i],
                              ncores = ncores,
                              mapping_ndim = mapping_ndim,
                              max_distance = max_distance,
                              verbose = verbose,
                              obs.method = obs.method,
                              pred.method = pred.method))
    }
    res[[i]] <- this_iter
    res[[i]]$error_threshold <- err[i]
  }

  names(res) <- paste0("iter",1:length(res))
  return(res)
}

