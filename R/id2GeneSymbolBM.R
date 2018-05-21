#' id2GeneSymbolBM
#'
#' Provided with a character vector of genomic identifiers; and a biomaRt
#' hsapiens_gene_ensembl specific set of filters and attributes (i.e. the
#' id format of the input and output respectively); this function returns
#' a data frame with the corresponding mappings of the given identifiers
#' to each of specified attributes.
#'
#' The 'format_in' argument defaults to ensembl gene id, and the mapped
#' 'format_out' argument defaults to hgnc symbol.
#'
#' Common filters/attributes are "ensembl_gene_id", "hgnc_symbol", 
#' "entrezgene", "illumina_humanht_12_v4" and "affy_hg_u133a".
#' 
#' Credit to Gil Tomas for the original code. Requires biomaRt and 
#' internet access.
#' 
#' @param ids A vector of genomic identifiers
#' @param format_in Format/type of genomic identifiers of the ids input
#' @param format_out Format/type of genomic identifiers to be returned
#' @param format_out Format/type of genomic identifiers to be returned
#' @return A 2-column dataframe of original and converted ids
#' @examples
#' id_vec <- c("ESR1", "ERBB2", "AURKA")
#' 
#' id2GeneSymbolBM(ids = id_vec,
#'                 format_in = "hgnc_symbol",
#'                 format_out = "entrezgene"
#'                 )
#' 
#' @export
id2GeneSymbolBM <- function(ids,
                            format_in = "ensembl_gene_id",
                            format_out = "hgnc_symbol"){
    #define mart to query
    ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", 
                                dataset = "hsapiens_gene_ensembl", 
                                host = "www.ensembl.org")
    #convert ids
    biomaRt::getBM(attributes = c(format_in, format_out),
                   filters = format_in,
                   values = ids,
                   mart = ensembl)
}
