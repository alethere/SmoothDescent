#Diagnostic tools ----------
#These functions aid in evaluating the output of smooth descent steps. They allow to:
# - count the number of recombinations per individual
# - calculate a reordering parameter between map iterations
# - calculate the "distance plot" for one or more iterations
# - plot multiple orderings of the same markers in an "iteration" plot.


#Recombination count ------------

#' Recombination counter
#'
#' Using IBD matrices, one can estimate the number of recombinations in a single chromosome. This
#' function calculates the number of recombinations using a single IBD matrix
#' or a list of IBD matrices of all parental homologues. Recombinations are
#' calculated per individual and per window. A genetic map is required.
#'
#' @param ibd IBD matrix or list of matrices. Must contain row.names indicating
#' the marker name.If `ibd` is a list, it must contain matrices for all parental homologues.
#' @param map data.frame with at least columns "position" and "marker". It is advisable
#' to use a single
#' @param non_inf Interval determining what is to be considered a non informative
#' probability. Values below non_inf[1] will be treated as 0 and above non_inf[2]
#' will be treated as 1. Values between the two will be ignored.
#'
#' @return A list containing:
#' * individual: recombination counts per individual. Provided per parent if
#' `ibd` was a list.
#' * on_map: a data.frame with recombination counts across individuals per
#' window. Provided per parent if `ibd` was a list.
#' @export
#'
#' @examples
rec_count <- function(ibd,map,non_inf = c(0.3,0.7)){
  if("data.frame" %in% class(ibd)) ibd <- as.matrix(ibd)
  UseMethod("rec_count",ibd)
}

#' Genotype calculator for a single IBD matrix
#'
#' @describeIn rec_count Recombination count for a single IBD matrix.
rec_count.matrix <- function(ibd,map,non_inf = c(0.3,0.7)){
  if(is.null(rownames(ibd))) stop("Row names are required in the ibd matrix")
  markers_present <- rownames(ibd) %in% map$marker
  if(!all(markers_present)){
    stop(paste("These markers in the ibd matrix were not found in the map dataframe:\n",
         paste0(rownames(ibd)[!markers_present],collapse = " ")))
  }

  #Recombination test
  map <- subset(map,marker %in% rownames(ibd))
  map <- map[order(map$position),]
  pos <- map$position
  names(pos) <- map$marker
  ibd <- ibd[map$marker,]
  ind_recs <- apply(ibd,2,function(ib_ind){
    ib_ind <- ib_ind[ib_ind < non_inf[1] | ib_ind > non_inf[2]]
    ib_ind <- round(ib_ind)
    rec_pos <- ib_ind[-1] != ib_ind[-length(ib_ind)]
    ind = sum(rec_pos)
    rec_map <- pos[names(rec_pos)[rec_pos]]
    return(list(ind = ind,rec_map = rec_map))
  })
  ind <- extract(ind_recs,"ind")
  rec_pos <- unlist(extract(ind_recs,"rec_map"))

  h <- hist(rec_pos,plot = F,breaks = seq(0,max(map$position)+5,by = 5))
  recs_on_map <- data.frame(start = h$breaks[-length(h$breaks)],
                            end = h$breaks[-1],
                            mid = h$mids,
                            count = h$counts)
  return(list(individual = ind,
              on_map = recs_on_map))

}


#' Recombination counter for a list of IBD matrices
#'
#' @describeIn rec_count Recombination count for a list of IBD matrices
#' for all parental homologues.
rec_count.list <- function(ibd,map,non_inf = c(0.3,0.7)){

  recs <- lapply(ibd,rec_count.matrix,map = map,non_inf = non_inf)

  parentindex <- list(p1=1:(length(ibd)/2),
                      p2 = (length(ibd)/2+1):length(ibd))
  recs <- lapply(parentindex,function(p){
    ind = rowSums(extract(recs[p],"individual"))/2
    counts <- sapply(extract(recs[p],"on_map",simplify = F),`[`,,"count")
    on_map <- data.frame(recs[p][[1]]$on_map[,1:3],count = rowSums(counts)/2)
    return(list(individual = ind,
                on_map = on_map))
  })
  return(recs)


}

#Reordering ----------


#' Reorder parameter calculator
#'
#' Given two maps with the same markers, calculates the
#' reordering parameter \tau between them. This is based on the
#' rank-correlation, or Kendall correlation.
#'
#' @param oldmap a data.frame containing at least "marker" and "position" columns.
#' @param newmap a data.frame containing at least "marker" and "position" columns.
#' All marker names must be present in `oldmap`.
#'
#' @return A single parameter indicating the rank order correlation.
#' @export
#'
#' @examples
reorder_tau <- function(oldmap,newmap){
  in_oldmap <- newmap$marker %in% oldmap$marker
  if(!all(in_oldmap)){
    stop(paste0("These markers were found in newmap but not in oldmap:\n",
               paste0(newmap$marker[!in_oldmap],collapse = " ")))
  }
  nm <- newmap
  rownames(nm) <- nm$marker
  tau <- cor(oldmap$position,nm[oldmap$marker,"position"],method = "kendall",use = "pair")
  return(tau)
}


#Recdist plot -----------

#' Recombination-distance plot calculator
#'
#' Given a map and a linkage data.frame with pairwise recombination
#' estimate, it returns a data.frame with all the information needed
#' to produce a recombination-distance plot. These plots are useful to
#' assess whether a given map order matches well with recombination estimates,
#' since larger recombination estimates should be equivalent to
#' larger map distances.
#'
#' @param map a data.frame containing at least "marker" and "position" columns.
#' @param linkdf a data.frame with columns "marker_a", "marker_b", "r" for
#' recombination and "LOD_independence" for the LOD score. Can be obtained
#' using the `linkage_shortcut` function and `polymapR::linkage` function.
#'
#' @return A data.frame with columns "recombination", for the recombination
#' frequency; "distance", for the physical distance between two markers; "marker_a"
#' for the name of the first marker; "marker_b" for the name of the second marker;
#' and "LOD" for the logarithm of odds score for that recombination estimate.
#'
#' @describeIn recdist_plot
#' @export
#'
#' @examples
recdist_calc <- function(map,linkdf){
  in_linkdf <- (map$marker %in% linkdf$marker_a) | (map$marker %in% linkdf$marker_b)
  if(!all(in_linkdf)){
    stop(paste0("Some markers found in map were not found in linkdf:\n",
               paste0(map$marker[!in_linkdf],collapse = " ")))
  }
  pos <- map$position
  names(pos) <- map$marker
  dist_vec <- apply(linkdf,1,function(linkrow){
    abs(pos[linkrow["marker_a"]] - pos[linkrow["marker_b"]])
  })
  recdist <- data.frame(recombination = linkdf$r,
                         distance = dist_vec,
                         marker_a = linkdf$marker_a,
                         marker_b = linkdf$marker_b,
                         LOD = linkdf$LOD_independence)

  return(recdist)
}

#' Recombination-distance plotter
#'
#' Plots a recombination-distance plot that allows to assess whether a given map order
#' matches well with recombination estimates,
#' since larger recombination estimates should be equivalent to
#' larger map distances. It is useful to colour the dots with
#' the LOD score.
#'
#'
#' @param recdist data.frame containing at least columns "recombination" and "distance".
#' Can be obtained from `recdist_calc`.
#' @param ... additional parameters passed onto `plot`
#'
#' @describeIn recdist_calc
#' @export
#'
#' @examples
recdist_plot <- function(recdist,...){
  sp <- lm(distance ~ poly(recombination,2),data = recdist)
  new_x <- seq(0,0.5,by = 0.01)
  pred <- predict(sp,newdata = data.frame(recombination = new_x))

  plot(recdist$recombination,recdist$distance,xlab = "recombination",ylab = "distance",...)
  points(new_x,pred, type = "l",col = "red3",lwd =  3)
}

#Iteration maps ---------

#' Plotting of map comparisons
#'
#' Given a list of maps with the same markers (re-orderings of the same markers)
#' it plots all maps with segments connecting the same marker as it moves through the
#' maps. Additionally, the reorder parameter tau (see `reorder_tau`) is also
#' printed between map comparisons.
#'
#' @param maplist A list of maps containing at least "position" and "marker" columns.
#' The map from a single chromosome is expected. Names of this list are used as
#' axis names.
#' @param ... Additional parameters to be passed to `plot()` such as `main`
#' for the title.
#'
#' @return A comparison plot between maps.
#' @export
#'
#' @examples
#' maplist <- list(pre = read.table("test/test_map.txt"),
#' flat = readRDS("test/test_flat_map.RDS")$locimap,
#' sphere = readRDS("test/test_sphere_map.RDS")$locimap)
#' maplist$flat$marker <- maplist$flat$locus
#' maplist$sphere$marker <- maplist$sphere$locus
#'
#' iterplot(maplist,main = "Example")
#'
iterplot <- function(maplist,...){
  have_pos <- sapply(maplist,hasName,"position")
  if(!all(have_pos)) stop("Not all maps have the position column")
  have_mark <- sapply(maplist,hasName,"marker")
  if(!all(have_mark)) stop("Not all maps have the marker column")

  nmaps <- length(maplist)
  #Calculate coordinates
  coo <- seq(0,(nmaps*2-1))
  coo <- data.frame(start = coo[-length(coo)],end = coo[-1])
  coo$mid <- (coo$start + coo$end)/2
  coo$left <- coo$start + 0.05
  coo$right <- coo$end - 0.05

  #Plot the maps
  yrange <- rev(range(extract(maplist,"position")))
  plot(0,type = "n",ylim = yrange, xlim = range(coo),
       xlab = "",ylab = "position",axes = F,...)

  for(i in coo$start){
    if(i%%2 == 0){
      #It's even we must plot a map
      map <- maplist[[i/2+1]]
      co <- coo[i+1,]
      rect(xleft = co$left,xright = co$right,
           ytop = min(map$position),ybot = max(map$position),
           col = "lightgrey",border = NA)
      segments(x0 = co$left, x1 = co$right,y0 = map$position)

    }else{
      #We must plot segments between maps
      map1 <- maplist[[(i-1)/2+1]]
      map2 <- maplist[[(i+1)/2+1]]
      co <- coo[i+1,]

      tau <- round(reorder_tau(map1,map2),2)
      if(tau < 0){
        warning("Map ",(i+1)/2+1," is plotted in reversed order")
        map2$position <- scales::rescale(map2$position,
                                         rev(range(map2$position)))
        maplist[[(i+1)/2+1]] <- map2
        tau <- round(reorder_tau(map1,map2),2)
      }
      map_match <- match(map2$marker,map1$marker)
      segments(x0 = co$left, x1 = co$right,
               y0 = map1$position, y1 = map2$position[map_match])
      text(x = co$mid,y = 0,pos = 3,labels = bquote(tau == .(tau)),
           xpd = T,cex = 0.75)
    }
  }

  mapnames <- names(maplist)
  if(is.null(mapnames)) mapnames <- paste0("map",1:length(maplist))
  axis(1, at = coo$mid[coo$start%%2 == 0],labels = mapnames)
  axis(2)
}

#IBD plots -------

#' Graphical genotype plotter
#'
#' Given an IBD matrix or a list of matrices, a graphical genotype
#' plot is returned. Graphical genotypes allow to assess which homologue
#' an individual has obtained, they highlight putative errors and show
#' non-informative genotypes.
#'
#' @param IBD IBD matrix or list of IBD matrices. If IBD is a list, the
#' names are used as subplot titles.
#' @param only_informative logical, should only informative markers be plotted?
#' False by default.
#' @param ... Additional parameters to be passed to `plot`. For instance
#' `cex.axis` or `main`.
#'
#' @return A plot containing graphical representations of IBD. Dark
#' is a probability close to 1 and light is a probability close to 0.
#' @export
#'
#' @examples
graphical_genotype <- function(IBD,only_informative = F,...){
  UseMethod("graphical_genotype",IBD)
}

#' Title
#'
#' @param IBD
#' @param cex.axis
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
graphical_genotype.matrix <- function(IBD,cex.axis = 0.8,...){
  image(t(IBD),axes = F,...)
  at <- seq(0,1,length.out = ncol(IBD))
  axis(1,labels = colnames(IBD),las = 2,at = at,cex.axis = cex.axis)
}

#' Title
#'
#' @param IBD
#' @param ex.axis
#' @param main
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
graphical_genotype.list <- function(IBD,ex.axis = 0.8,main = NULL,...){
  par(mfrow = c(1,length(IBD)),
      mar = c(4,1,3,1),
      oma = c(0,0,3,0))
  for(i in seq_along(IBD)){
    ib <- IBD[[i]]
    graphical_genotype.matrix(ib,cex.axis = cex.axis,main = names(IBD)[i])
  }
  mtext(3,outer = T,text = main,cex = 1.3)
}


#Reorder tau
tau <- reorder_tau(maplist$pre,maplist$flat)
if(round(tau,4) != 0.89) stop("Obtained reorder tau and stored reorder tau do not coincide")

#Recdist plot
map <- read.table("test/test_map.txt")
linkdf <- readRDS("test/test_linkdf.RDS")
recdist <- recdist_calc(map,linkdf)
recdist_plot(recdist)


#Iterplot
maplist$flat$marker <- maplist$flat$locus
maplist$sphere$marker <- maplist$sphere$locus
iterplot(maplist,main = "Testing")


