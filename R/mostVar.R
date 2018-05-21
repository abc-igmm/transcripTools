#' mostVar
#' Calculate the top n most variable genes in a matrix of gene expression data
#'
#' @param data An expression matrix (samples in columns, genes in rows)
#' @param n The number of variable genes to return
#' @param i_want_most_var Logical value to return most variable (TRUE) or least
#' variable (FALSE) genes
#' @return A subset of n most variable rows from the original data
#' @examples
#' #synthesise example data matrix
#' mtx <- matrix(rnorm(mean = 7, n = 90), ncol = 10) 
#' row.names(mtx) <- LETTERS[1:9]
#'
#' mv_3 <- mostVar(data = mtx, n = 3)
#'
#' @export
mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}
