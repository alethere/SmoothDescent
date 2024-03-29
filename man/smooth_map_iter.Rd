% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smooth_wrapper.R
\name{smooth_map_iter}
\alias{smooth_map_iter}
\title{Iterative Smooth Mapping}
\usage{
smooth_map_iter(
  geno,
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
  non_inf = c(0.3, 0.7),
  verbose = T,
  obs.method = "naive",
  pred.method = "prediction",
  hmm.error = 0.01
)
}
\arguments{
\item{geno}{matrix with markers on the rows, individuals on the columns.
Rownames and columnames expected.}

\item{homologue}{matrix with markers on the rows, homologue names on the columns.
Rownames and columnames expected.}

\item{iters}{numeric indicating the number of iterations to perform.}

\item{map}{optionally, data.frame with at least columns "marker" and "position". If it is
not specified, a map will be estimated from the uncorrected genotype data with polymapR.}

\item{ploidy}{numeric indicating the ploidy. Both parents must be of the same ploidy,}

\item{p1name}{character, name of the first parent. Must be present in the geno columnames.
If it's not specified it will be taken as the name of the first column.}

\item{p2name}{character, name of the second parent. Must be present in the geno columnames.
If it's not specified it will be taken as the name of the second column.
and it is assumed that "homologue" has 2*ploidy columns.}

\item{prediction_interval}{numeric, interval to be used during the IBD prediction
step. It should be specified in the same units as the "position" in the map.}

\item{prediction_threshold}{float, probability threshold for imputing new genotypes.
All new genotypes with a probability under this threshold will be considered uncertain.
Defaults to 0.8.}

\item{prediction_points}{numeric, number of points to use for IBD prediction. If NULL, all
points in map$position are used, otherwise n equally spaced points are used. Greatly
improves efficiency if the number of markers is very large.}

\item{error_threshold}{numeric, threshold over which a marker is considered erroneous.
Usually 0.8 should be good enough to be sensitive but stringent (not have false positives)}

\item{ncores}{number of cores to use for linkage estimation.}

\item{mapping_ndim}{2 or 3, number of dimensions to use for multi-dimensional mapping}

\item{estimate_premap}{logical, whether to use polymapR to estimate a preliminary map.}

\item{max_distance}{numeric, markers that have near neighbours will be eliminated. This
parameter defines the maximum neighbour distance allowed. A warning will be issued if some markers
are eliminated.}

\item{non_inf}{numeric, lower and upper probability boundaries to consider an
IBD probability non-informative (if they fall within the threshold they will be
ignored during prediction). Defaults to 0.3 - 0.7. Symmetrical boundaries are
recommended but not necessary.}

\item{verbose}{logical, should smooth descent report the steps it takes?}

\item{obs.method}{character, either "naive" or "heuristic" (or substrings). This parameter allows to
switch between using the IBD calculation (for observed IBDs) described in the Smooth Descent paper, or the
heuristic method from `polyqtlR`. However, our research has shown better results with the
naive method.}

\item{pred.method}{character, either "prediction" or "hmm" (or substrings). This parameter
allows to switch between using the IBD calculation (for predicted IBDs) between the weighted
average method or the Hidden Markov model implemented in `polyqtlR`. Our research shows better
results in polyploids with the HMM, although for high marker densities the
weighted average method is faster (specially if `prediction_points` is used)}
}
\value{
list of lists, where each element contains the following items:
* obsIBD: list of observed IBD matrices (marker x individual) for each parental homologue
* predIBD: list of predicted IBD matrices (marker x individual) for each parental homologue
* oldmap: data.frame containing the original map
* error: list of error matrices (marker x individual) for each parental homologue
* newIBD: list of observed IBD matrices of the corrected genotypes.
* newmap: data.frame containing the new updated map.
* newgeno: matrix containing the new genotypes.
* rec: list containing the recombination counts using the observed IBD (`obs`),
the predicted IBD (`pred`) and the updated IBD (`new`). `obs` and `pred` use the old
map order while `new` uses the re-estimated map. For more information see `rec_count()`.
* recdist: data.frame containing pair-wise recombination and distance between markers,
useful to plot using `recdist_plot()`
* r2: R-squared parameter of pair-wise recombination and final map distance. A higher
value indicates a better newmap.
* tau: reordering parameter tau (Kendall's rank correlation). Obtained with `reorder_tau()`
* eliminated: markers eliminated due to a large average neighbour distance, they tend
to be problematic when re-mapping.
}
\description{
This function allows to apply multiple interations of the Smooth Descent
algorithm, including a re-mapping step using new genotypes. Optionally, the
preliminary map can also be estimated instead of being provided. Both mapping
procedures are performed with packages `polymapR` and `MDSMap`, which perform
multi-dimensional scaling mapping. For usage of the Smooth Descent algorithm
with a different map algorithm see `smooth_descent()`, which performs the
genotype correction without re-mapping.
}
\examples{
\dontrun{

data("genotype")
data("homologue")
data("map")

res <- smooth_descent(geno,hom,iters = 5, ploidy = 2, p1name = "P1", p2name = "P2",
   estimate_premap = F, mapping_ndim = 3, ncores = 1)

}
}
