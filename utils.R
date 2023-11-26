# utils.R
# Functions for supporting utilities in the DynaMICE model (Dynamic Measles Immunisation Calculation Engine)

# ------------------------------------------------------------------------------
#' Expand the matrix to a different dimension
#
#' This function returns an expanded contact matrix with the specified age
#' structure for model inputs.
# ------------------------------------------------------------------------------
#' @param A A matrix to be expanded.
#' @param expand_rows Number of times to repeat each row.
#' @param expand_cols Number of times to repeat each column.
#' @param rescale_rows A logical variable to control whether to re-scale the
#' expanded rows.
#' @param rescale_cols A logical variable to control whether to re-scale the
#' expanded columns.
#' @return An expanded matrix with re-scaling if applicable.
#' @examples
#' expandMatrix (matrix(1:9,3,3), 2, 1, FALSE, FALSE)
expandMatrix <- function (A,
                          expand_rows  = 1,
                          expand_cols  = 1,
                          rescale_rows = F,
                          rescale_cols = F) {
  
  if(!is.matrix(A)){
    stop("A is not a matrix")
  }
  
  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)
  
  for(c in 1:cols) {
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }
  
  B <- matrix (matvals,
               nrow = rows * expand_rows,
               ncol = cols * expand_cols)
  
  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols)
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }
  
  return (B)
  
} # end of function -- expandMatrix
# ------------------------------------------------------------------------------

