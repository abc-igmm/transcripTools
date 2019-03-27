CreateAnnotationDataframe_entrez <- function(expressionDF){
  df <- data.frame(probe=row.names(expressionDF), EntrezGene.ID=row.names(expressionDF))
  row.names(df) <- df$EntrezGene.ID
  df
}
