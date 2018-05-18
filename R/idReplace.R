#' A function to convert a matrix/dataframes's genomic identifiers (found as row
#' names) 
#' 
#' Gene ids are converted before being substituted for the original ids in the 
#' input data.
#'
#' Provided with a matrix/dataframe of gene expression with genomic identifiers; 
#' and a biomaRt hsapiens_gene_ensembl specific set of filters and attributes 
#' (i.e. the id format of the input and output respectively); this function 
#' returns the original data object with corresponding mappings of the given 
#' replacement identifiers.
#' 
#' Multiple input ids mapping to a single unique output id are aggregated by 
#' mean expression.
#'
#' The 'format_in' argument defaults to ensembl gene id, and the mapped
#' 'format_out' argument defaults to hgnc symbol.
#'
#' Common filters/attributes are "ensembl_gene_id", "hgnc_symbol", 
#' "entrezgene", "illumina_humanht_12_v4" and "affy_hg_u133a".
#' 
#' Requires biomaRt and internet access.
#' @inheritParams id2GeneSymbolBM
#' @param input A matrix/dataframe. Samples as columns, genes as rows
#' @param format_in Format/type of genomic identifiers of the input object
#' @return A id-converted matrix otherwise identical to input
#' @examples
#' 
#' mtx <- matrix(rnorm(mean = 7, n = 9), ncol = 3) 
#' row.names(mtx) <- c("ESR1", "ERBB2", "AURKA")
#' 
#' idReplace(mtx,
#'           format_in = "hgnc_symbol",
#'           format_in = "entrezgene")
#' 
#' @export
idReplace <- function(input,
                      format_in = "ensembl_gene_id",
                      format_out = "hgnc_symbol"){

  #Function to retrieve corresponding new ids for original ids
  id2GeneSymbolBM <- function(ids,
                              filters = format_in,
                              attributes = format_out){
    library(biomaRt)
    ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "www.ensembl.org")
    getBM(attributes = c(filters, attributes),
          filters = filters,
          values = ids,
          mart = ensembl)
  }

  #Merge and aggregate (mean) multiple IDs per format_out
  input <- as.data.frame(input)
  input[[format_in]] <- row.names(input)

  input.mrg <- merge(x = id2GeneSymbolBM(row.names(input),
                                      filters = format_in,
                                      attributes = format_out),
                  y = input,
                  by = format_in)

  input.mrg[input.mrg==""] <- NA
  input.mrg <- na.omit(input.mrg)
  ## added data.table conversion to speed up aggregation step [defunct]
  #input.mrg <- data.table(input.mrg)

  input.agg <- aggregate(input.mrg,
                        list(input.mrg[[format_out]]),
                        FUN = mean)

  #input.agg <- input.mrg[, Group.1 := mean(x), by = get(format_out)]


  #Set row names and tidy now-redundant columns
  row.names(input.agg) <- input.agg$Group.1

  #data.table aggregation makes this unecessary
  input.agg$Group.1 <- NULL
  input.agg[[format_in]] <- NULL
  input.agg[[format_out]] <- NULL

  input.agg
}
