#IBD prediction ------------

#' Predict IBD
#'
#' IBD prediction algorithm, performs a distance-informed weighted average using
#' a genetic map. IBD probabilities between 0.3 and 0.7 are ignored as non-informative.
#'
#' @param IBD observed IBD. Can be given as a vector containing the
#' observed IBDs of one individual, a single matrix (markers x individual) or a list
#' of matrices. Matrices must contain rownames and columnames.
#' @param map data.frame containing "marker" and "position" columns at least.
#' @param interval numeric indicating the interval size to be used in prediction.
#'
#'
#' @return the output will have the same format as the input given. That is, an input
#' IBD vector will return an IBD vector, a matrix will return a matrix and a list of matrices
#' will return a list of matrices.
#' @export
#'
#' @examples
predict_IBD <- function(IBD,map,interval = 10){
  type <- class(IBD)
  if(type[1] == "list"){
    types <- sapply(IBD,function(i) is.matrix(i) | is.numeric(i))
    if(all(types)){
      UseMethod("predict_IBD",IBD)
    }else{
      stop("This type of list is not supported. All elements must be matrices or named numerics")
    }
  }else if(type[1] == "data.frame"){
    predict_IBD.matrix(IBD,map,interval)
  }

  UseMethod("predict_IBD",IBD)
}

#Method to apply IBD prediction to a list of IBD matrices from the same map.
#It returns a list of predicted IBD matrices.

#' IBD prediction for a list of IBD matrices
#'
#' @describeIn predict_IBD Function for a list of matrices
#' @export
predict_IBD.list <- function(IBD,map,interval = 10){
  lapply(IBD,predict_IBD, map = map, interval = interval)
}

#Method to apply IBD prediction to a single matrix. It calls upon the numeric method.
#It returns a predicted IBD matrix.

#' IBD prediction for an IBD matrix
#'
#' @describeIn predict_IBD Function for a single matrix
#' @export
predict_IBD.matrix <- function(IBD,map,interval = 10){
  if(nrow(IBD) == 0) return(NULL)
  if(is.null(rownames(IBD))) stop("IBD matrix should contain rownames")
  if(!"marker" %in% colnames(map)) stop("Map data.frame should contain a \"marker\" column")
  if(!"position" %in% colnames(map)) stop("Map data.frame should contain a \"position\" column")

  #This can be done faster in a different way
  #pred <- apply(IBD,2,predict_IBD,map = map, interval = interval)

  #First compute the markers that should go in each interval
  IBD <- IBD[map$marker,]
  IBD[IBD > 0.3 & IBD < 0.7] <- NA

  pred <- lapply(1:nrow(map),function(i){
    #Compute weights for each marker (they are the same accross columns)
    distance <- abs(map$position[i]-map$position)
    weights <- 1 - distance/100
    close_marks <- distance < interval & distance != 0

    #The compute each row at once
    ib <- IBD[close_marks,,drop = F]
    weights <- weights[close_marks,drop = F]
    wibd <- ib*weights
    wmat <- matrix(rep(weights,ncol(ib)),ncol = ncol(ib))
    wmat[is.na(ib)] <- NA
    colSums(wibd,na.rm = T)/colSums(wmat,na.rm = T)
  })
  pred <- do.call(rbind,pred)
  rownames(pred) <- rownames(IBD)
  return(pred)
}

#' IBD prediction for a a numeric vector
#'
#' @describeIn predict_IBD Function for a numeric vector
#' @export
predict_IBD.numeric <- function(IBD,map,interval = 10){
  if(!"position" %in% colnames(map)) stop("Map data.frame should contain a \"position\" column")
  if(!"marker" %in% colnames(map)) stop("Map data.frame should contain a \"marker\" column")
  if(is.null(names(IBD))) stop("IBD vector should be named")
  marks_in_map <- names(IBD) %in% map$marker
  if(!all(marks_in_map)) stop("Some names in IBD vector were not found in map$markers: ",
                              paste(names(IBD)[!marks_in_map],collapse = ", "))

  # IBD[IBD == non_inf] <- NA
  IBD[IBD > 0.3 & IBD < 0.7] <- NA
  #To ensure position and IBD vectors are in the same order
  position <- map$position
  names(position) <- map$marker
  position <- position[names(IBD)]
  sapply(1:length(IBD),function(i){
    distance <- abs(position - position[i])
    weights <- 1-distance[distance < interval & !is.na(IBD)]/100
    weights <- weights/sum(weights)
    sum(IBD[distance < interval & !is.na(IBD)]*weights)
  })
}
