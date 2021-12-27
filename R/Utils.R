#Utilities --------------
#Useful little function to extract same-name elements from a list
#' Extraction from list
#'
#' Given a list of identical objects, extract the same index of each object.
#' This can be used to extract the same column from a list of matrices,
#' the same index from a list of vectors or the same element from a list of lists.
#' Multiple indeces can be provided.
#'
#' @param l list of elements with identical structure.
#' @param index numeric or character to be applied to each element of list.
#' Logicals can also be used, but will not work well with lists of lists.
#' @param simplify should the output be simplified?
#'
#' @return list, vector or matrix of extracted objects
#' @export
#'
#' @examples
#' example <- list(data.frame(A = "first dataframe", B = "second argument"),
#'     data.frame(A = "second dataframe", B = "second argument"))
#' #From each data.frame extract column A
#' extract(example,"A")
extract <- function(l,index = 1,simplify = T){
  if(length(index) > 1){
    res <- lapply(index,function(i){
      extract(l,i,simplify = T)
    })
    names(res) <- index
    return(res)
  }else{
    sapply(l,function(element){
      if(is.list(element)){
        element[[index]]
      }else if(is.matrix(element)){
        element[,index]
      }else if(is.vector(element)){
        element[index]
      }
    },simplify = simplify)
  }
}

#Function to cbind two matrices based on columnames. Similar to merge.

#' Matrix merging of rows with equal names
#'
#' Similar to merge, given two matrices, row names are matched
#' and columns are bound. All columns are preserved and all row.names
#' present in both matrices are present in the result. Empty cells,
#' belonging to a row that is present in one matrix and not the other
#' are filled in with the `fill` value-
#'
#' @param mat1 matrix 1, with dimnames.
#' @param mat2 matrix 2, with dimnames.
#' @param fill value to put inside the output matrix for those rows
#' that were present in only one of the two matrices.
#'
#' @return fused matrix
#' @export
rowfuse <- function(mat1,mat2,fill = NA){
  rows <- unique(c(rownames(mat1),rownames(mat2)))
  resmat <- matrix(fill,ncol = ncol(mat1) + ncol(mat2),
                   nrow = length(rows))
  rownames(resmat) <- rows
  colnames(resmat) <- c(colnames(mat1),colnames(mat2))
  resmat[rownames(mat1),1:ncol(mat1)] <- mat1
  resmat[rownames(mat2),(ncol(mat1)+1):ncol(resmat)] <- mat2
  return(resmat)
}
