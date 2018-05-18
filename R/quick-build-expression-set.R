#' quickBuildeSet
#' Build an expressionSet from matched expression and phenotypic matrices
#' @param xpr Expression matrix, samples in columns, genes in rows. Column order
#' must be identical to pheno row order
#' @param pheno Phenotypic data matrix, samples in rows, variables in columns. 
#' Row order must be identical to xpr column order
#' @return An eSet of the combine matrices
#' @export
quickbuildeSet <- function(xpr, pheno){
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
