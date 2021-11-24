#Dataset documentation

#' Homologue matrix
#'
#' A homologue matrix is a markers x homologues matrix, with marker names
#' as row names and homologue names (as many as ploidy for each parent) for
#' columns. This example is from one chromosome of a diploid organism (strabwerry).
#'
#' @name hom
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("homologue")
#' @keywords data
NULL

#' Map data.frame
#'
#' A map data.frame should contain at least two columns: marker (with marker names)
#' and position, with physical or genetic positions. Each row corresponds to a marker.
#' This example is from one chromosome of a diploid organism (strabwerry).
#'
#' @name map
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("map")
#' @keywords data
NULL

#' Linkage data.frame
#'
#' This linkage data.frame contains pairwise linkage information as obtained by
#' `linkdf_shortcut()`, which calculates linkage using polymapR. It contains coluns
#' marker_a, marker_b, and LOD_independence, which is used (rather than LOD score).
#'
#' @name linkdf
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("linkage")
#' @keywords data
NULL

#' Genotype matrix
#'
#' A genotype matrix is a markers x individuals matrix, with marker names as row names
#' and individual names as column names. Dosages (as integers) are expected.
#' This example is from one chromosome of a diploid organism (strabwerry).
#'
#' @name geno
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("genotype")
#' @keywords data
NULL

#' Smooth Descent output
#'
#' Example output from a single run of the Smooth Descent algorithm (using `smooth_descent()`).
#' The example data from this package were used (`geno`,`hom`,`map`). The rest of
#' parameters were kept on default.
#'
#' @name sdescent
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("smooth_descent")
#' @keywords data
NULL

#' Smooth Descent iterations output
#'
#' Example output from 5 iterations of Smooth Descent (using `smooth_map_iter()`).
#' The example data from this package were used (`geno`,`hom`,`map`). Mapping was
#' performed with `ndim = 3`. The rest of parameters were kept on default.
#'
#' @name sd_iter
#' @docType data
#' @author Alejandro Therese Navarro \email{alejandro.theresenavarro@@wur.nl}
#' @usage data("iter_smooth_descent")
#' @keywords data
NULL
