#IBD prediction ------------

#' Predict IBD
#'
#' IBD prediction algorithm, performs a distance-informed weighted average using
#' a genetic map. Weights are calculated as a proportion of the non-recombination
#' probability between two markers, which are estimated using mapping functions (Kosambi's,
#' Haldane's or Morgan's). Low-informative IBD probabilities are ignored, by default probabilities
#' between 0.3 and 0.7 are considered low-informative.
#'
#' @param IBD observed IBD. Can be given as a vector containing the
#' observed IBDs of one individual, a single matrix (markers x individual) or a list
#' of matrices. Matrices must contain rownames and columnames.
#' @param map data.frame containing "marker" and "position" columns at least.
#' @param interval numeric indicating the interval size to be used in prediction.
#' @param method string indicating the mapping function used to calculate weights.
#' Either "kosambi", "haldane" or "morgan". Kosambi takes into account chiasma interference,
#' Haldane is the most common mapping function and Morgan assumes a linear relationship
#' between distance and recombination (is inaccurate for large distances).
#' @param non_inf numeric vector of two digits, containing the lower and upper
#' bound for probabilities to be considered non-informative. These probabilities will be
#' ignored during prediction. By default 0.3 and 0.7. Symmetric and stringent thresholds
#' are recommended.
#' @param pred_points numeric, number of points to use for IBD prediction. If NULL, all
#' points in map$position are used, otherwise n equally spaced points are used. Greatly
#' improves efficiency if the number of markers is very large.
#'
#'
#' @return the output will have the same format as the input given. That is, an input
#' IBD vector will return an IBD vector, a matrix will return a matrix and a list of matrices
#' will return a list of matrices.
#' @export
#'
#' @examples
#' data("genotype")
#' geno <- geno[,-1:-2] #we take out parental genotypes
#' data("homologue")
#' data(map)
#' IBD <- calc_IBD(geno,hom[,1:2],hom[,3:4], ploidy = 2)
#'
#' #One homologue of one individual
#' pred <- predict_IBD(IBD[[1]][,1], map)
#' #One homologue for all individuals
#' pred <- predict_IBD(IBD[[1]], map)
#' #Or all homologues
#' pred <- predict_IBD(IBD, map)
predict_IBD <- function(IBD,map,interval = 10, method = "kosambi",
                        non_inf = c(0.3,0.7),pred_points = NULL){
  type <- class(IBD)
  if(type[1] == "list"){
    types <- sapply(IBD,function(i) is.matrix(i) | is.numeric(i))
    if(all(types)){
      UseMethod("predict_IBD",IBD)
    }else{
      stop("This type of list is not supported. All elements must be matrices or named numerics")
    }
  }else if(type[1] == "data.frame"){
    predict_IBD(IBD,map,interval, method = method,
                       non_inf = non_inf,pred_points = pred_points)
  }

  UseMethod("predict_IBD",IBD)
}

#Method to apply IBD prediction to a list of IBD matrices from the same map.
#It returns a list of predicted IBD matrices.

#' IBD prediction for a list of IBD matrices
#'
#' @describeIn predict_IBD Function for a list of matrices
#' @export
predict_IBD.list <- function(IBD,map,interval = 10,method = "kosambi",
                             non_inf = c(0.3,0.7), pred_points = NULL){
  lapply(IBD,predict_IBD, map = map, interval = interval,
         method = method, non_inf = non_inf,pred_points = pred_points)
}

#Method to apply IBD prediction to a single matrix. It calls upon the numeric method.
#It returns a predicted IBD matrix.

#' IBD prediction for an IBD matrix
#'
#' @describeIn predict_IBD Function for a single matrix
#' @export
predict_IBD.matrix <- function(IBD,map,interval = 10,method = "kosambi", non_inf = c(0.3,0.7),
                               pred_points = NULL){
  if(nrow(IBD) == 0) return(NULL)
  if(is.null(rownames(IBD))) stop("IBD matrix should contain rownames")
  if(!"marker" %in% colnames(map)) stop("Map data.frame should contain a \"marker\" column")
  if(!"position" %in% colnames(map)) stop("Map data.frame should contain a \"position\" column")
  if(length(non_inf) != 2) stop("non_inf must be a numeric vector of length 2")
  if(non_inf[1] > non_inf[2]) stop("non_inf[1] is larger than non_inf[2]")
  if(any(non_inf > 1) | any(non_inf < 0)){
    stop("non_inf must contain numbers between 0 and 1")}

  #This can be done faster in a different way
  #pred <- apply(IBD,2,predict_IBD,map = map, interval = interval)

  #First compute the markers that should go in each interval
  IBD <- IBD[map$marker,]
  IBD[IBD > non_inf[1] & IBD < non_inf[2]] <- NA

  #Here parallelism could be added
  #Or an estimation procedure based on points rather than markers
  if(!is.null(pred_points)){
    est_points <- seq(from = min(map$position),
                      to = max(map$position),
                      length.out = pred_points)
  }else{
    #By default we estimate per each point on the map
    est_points <- map$position
  }

  pred <- lapply(est_points,function(p){
    #First we compute distance of all markers to the estimation point
    distance <- abs(p - map$position)
    #Then we obtain the weights and the nearby markers
    weights <- 1 - r_from_dist(distance,method)
    close_marks <- distance < interval & distance != 0

    #Then compute each row at once
    ib <- IBD[close_marks,,drop = F]
    weights <- weights[close_marks,drop = F]
    wibd <- ib*weights
    wmat <- matrix(rep(weights,ncol(ib)),ncol = ncol(ib))
    wmat[is.na(ib)] <- NA
    colSums(wibd,na.rm = T)/colSums(wmat,na.rm = T)
  })
  pred <- do.call(rbind,pred)

  #If we use point-based estimation, we must re-dimension result
  #to fit the number of markers in the map
  if(!is.null(pred_points)){
    closest_est_points <- sapply(map$position,function(p){
      which.min(abs(est_points - p))
    })
    pred <- pred[closest_est_points,]
  }
  rownames(pred) <- rownames(IBD)

  return(pred)
}

#' IBD prediction for a a numeric vector
#'
#' @describeIn predict_IBD Function for a numeric vector
#' @export
predict_IBD.numeric <- function(IBD,map,interval = 10,
                                method = "kosambi",
                                non_inf = c(0.3,0.7),
                                pred_points = NULL){
  if(!"position" %in% colnames(map)) stop("Map data.frame should contain a \"position\" column")
  if(!"marker" %in% colnames(map)) stop("Map data.frame should contain a \"marker\" column")
  if(is.null(names(IBD))) stop("IBD vector should be named")
  marks_in_map <- names(IBD) %in% map$marker
  if(!all(marks_in_map)) stop("Some names in IBD vector were not found in map$markers: ",
                              paste(names(IBD)[!marks_in_map],collapse = ", "))
  if(length(non_inf) != 2) stop("non_inf must be a numeric vector of length 2")
  if(non_inf[1] > non_inf[2]) stop("non_inf[1] is larger than non_inf[2]")
  if(any(non_inf > 1) | any(non_inf < 0)){
    stop("non_inf must contain numbers between 0 and 1")}

  # IBD[IBD == non_inf] <- NA
  IBD[IBD > non_inf[1] & IBD < non_inf[2]] <- NA
  #To ensure position and IBD vectors are in the same order
  IBD <- IBD[map$marker]

  #Or an estimation procedure based on points rather than markers
  if(!is.null(pred_points)){
    est_points <- seq(from = min(map$position),
                      to = max(map$position),
                      length.out = pred_points)
  }else{
    #By default we estimate per each point on the map
    est_points <- map$position
  }

  #Here parallelism could be added
  pred <- sapply(est_points,function(p){
    distance <- abs(map$position - p)
    local_d <- distance < interval & distance != 0
    weights <- 1 - r_from_dist(distance[local_d],method = method)
    inf <- !is.na(IBD[local_d])
    weights <- weights[inf]/sum(weights[inf])
    sum(IBD[local_d][inf]*weights)
  })
  #If we use point-based estimation, we must re-dimension result
  #to fit the number of markers in the map
  if(!is.null(pred_points)){
    closest_est_points <- sapply(map$position,function(p){
      which.min(abs(est_points - p))
    })
    pred <- pred[closest_est_points]
  }
  names(pred) <- names(IBD)
  return(pred)
}

#Reverse mapping functions --------

#' Reverse mapping function
#'
#' Function to calculate recombination frequencies from distance based on three
#' mapping functions. Morgan's, which assumes a linear relationship between
#' distance and recombination frequency (inaccurate for large distances); Haldane's,
#' the most common function, accounting for possible double recombinations between
#' distant markers; and Kosambi's, which accounts for intereference between nearby
#' recombinations. Kosambi's is used by default, as it tends to be the most biologically
#' accurate.
#'
#' @param d numeric, distance vector in cM.
#' @param method character, one of the following: kosambi, haldane, morgan. Partial
#' matching is allowed.
#'
#' @return a numeric vector of recombination frequencies.
#' @export
#'
#' @examples
#' d <- 1:100
#' kos <- r_from_dist(d, method = "k")
#' hal <- r_from_dist(d, method = "h")
#' mor <- r_from_dist(d, method = "m")
#'
#' plot(d,mor,ylim = range(c(kos,hal,mor)),type = "l", col = 1)
#' points(d,hal,type = "l", col = 2)
#' points(d, kos, type = "l",col = 3)
#' legend("topright", legend = c("Morgan","Haldane","Kosambi"),col = 1:3,lty = 1)
r_from_dist <- function(d,method = "kosambi"){
  method <- match.arg(method,choices = c("haldane","morgan","kosambi"))

  if(method == "haldane"){
    r <- (1 - exp(-d/50))/2
  }else if(method == "morgan"){
    r <- d/100
  }else if(method == "kosambi"){
    r <- ((exp(d/25) - 1)/(exp(d/25) + 1))/2
  }

  return(r)
}


