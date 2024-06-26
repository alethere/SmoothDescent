#Mapping functions ---------
#These functions are mostly made to simplify the process of using polyploid mapping
#procedures defined in MDSmap


#' MDSmap caller
#'
#' Shortcut function to launch MDSmap
#'
#' @param linkdf Linkage data.frame as defined in the `MDSMap::estimate.map()` function.
#' @param ndim 2 or 3, number of dimensions to use in MDS-mapping.
#'
#' @return Estimated linkage map, including position and other parameters.
#' @keywords internal
mdsmap <- function(linkdf,ndim = 2){
  linkdf <- polymapR:::prepare_pwd(linkdf)
  file <- paste0("tmp.",paste0(sample(c(letters,LETTERS),10),collapse = ""),".map")
  count <- length(unique(unlist(linkdf[,1:2])))
  write(count,file = file)
  write.table(linkdf, file = file, append = T,
              col.names = F, row.names = F, quote = F)

  newmap <- MDSMap::calc.maps.pc(file,weightfn = "lod2",ndim = ndim)
  rem <- file.remove(file)
  newmap <- newmap$locimap
  colnames(newmap)[2] <- "marker"

  return(newmap)
}

#' Linkdf shortcut
#'
#' Shortcut to calculate the linkage data.frame using polymapR's approach
#' to autopolyploid linkage calculation. Random pairing is assumed, it only
#' considers one linkage group at a time and returns linkage both from
#' parent1 to parent2 and viceversa. Importantly, markers for which
#' parental genotype does not match offspring segregation are eliminated.
#'
#' @param geno numeric data.frame or matrix, with marker names as rownames and
#' individual names as column names.
#' @param ploidy numeric, ploidy of both parents.
#' @param p1name character, name of the first parent.
#' @param p2name characther, name of the second parent.
#' @param ncores numeric, number of threads to use in the linakge calculation.
#'
#' @return a linkage dataframe. Some estimates might be duplicated,
#' since they are calculated for both parents.
#'
linkdf_shortcut <- function(geno,ploidy,p1name,p2name,ncores = 1){
  ma <- data.frame(markers = rownames(geno),Assigned_LG = 1)
  rownames(ma) <- ma$markers
  linkdf_p1 <- polymapR::finish_linkage_analysis(dosage_matrix = as.matrix(geno),
                                                 parent1 = p1name,
                                                 parent2 = p2name,
                                                 ploidy = ploidy,
                                                 G2_test = T,
                                                 ncores = ncores,
                                                 marker_assignment = ma,
                                                 LG_number = 1,verbose = F,
                                                 convert_palindrome_markers = T)$LG1
  linkdf_p2 <- polymapR::finish_linkage_analysis(dosage_matrix = as.matrix(geno),
                                                 parent1 = p2name,
                                                 parent2 = p1name,
                                                 ploidy = ploidy,
                                                 G2_test = T,
                                                 ncores = ncores,
                                                 marker_assignment = ma,
                                                 LG_number = 1,verbose = F,
                                                 convert_palindrome_markers = T)$LG1

  linkdf <- rbind(linkdf_p1,linkdf_p2)
  return(linkdf)
}

