#' quickBuildeSet
#'
#' Build an minimal expressionSet from matched expression and phenotypic 
#' matrices. Caution should be employed as this eSet will lack feature data and
#' proper phenotypic data variable descriptions.
#' @param xpr Expression matrix, samples in columns, genes in rows. Column order
#' must be identical to pheno row order
#' @param pheno Phenotypic data matrix, samples in rows, variables in columns. 
#' Row order must be identical to xpr column order
#' @return An eSet of the combine matrices
#' @examples
#' # Here we'll first dismantle the NKI dataset before re-combining as a minimal
#' # eSet
#' library(breastCancerNKI)
#' library(Biobase)
#' data(nki)
#' xpr <- exprs(nki)
#' pheno <- pData(nki)
#' 
#' minimal_nki <- quickBuildeSet(xpr, pheno)
#' # note that feature data and variable descriptions are now missing
#' @export
quickBuildeSet <- function(xpr, pheno){
    library(Biobase)
    library(testthat)
    test_that("dimnames are correct between xpr and pheno", {
                  expect_identical(row.names(pheno), colnames(xpr))
                        })
    metadata <- data.frame(labelDescription = colnames(pheno), 
                           row.names = colnames(pheno))
    phenoData <- new("AnnotatedDataFrame", data = pheno, varMetadata = metadata)
    ExpressionSet(
                  assayData = as.matrix(xpr),
                  phenoData = phenoData
                  )
}
