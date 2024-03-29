% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{extract}
\alias{extract}
\title{Extraction from list}
\usage{
extract(l, index = 1, simplify = T)
}
\arguments{
\item{l}{list of elements with identical structure.}

\item{index}{numeric or character to be applied to each element of list.
Logicals can also be used, but will not work well with lists of lists.}

\item{simplify}{should the output be simplified?}
}
\value{
list, vector or matrix of extracted objects
}
\description{
Given a list of identical objects, extract the same index of each object.
This can be used to extract the same column from a list of matrices,
the same index from a list of vectors or the same element from a list of lists.
Multiple indeces can be provided.
}
\examples{
example <- list(data.frame(A = "first dataframe", B = "second argument"),
    data.frame(A = "second dataframe", B = "second argument"))
#From each data.frame extract column A
extract(example,"A")
}
