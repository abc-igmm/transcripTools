CreateAnnotationDataframe_hgnc <- function(expressionDF){
  df <- data.frame(probe=row.names(expressionDF), Gene.Symbol=row.names(expressionDF))
  row.names(df) <- df$Gene.Symbol
  df
}