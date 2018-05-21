#' heatmapArrange 
#' 
#' A function to take a matrix/dataframe of gene expression data and return a 
#' dataframe suitable for plotting as a heatmap with ggplot2
#' 
#' @param data_in is a dataframe of gene expression data
#' @param cluster_row A logical value for whether row clustering is performed
#' @param cluster_column A logical value for whether column clustering is 
#' performed
#' @param scale A logical value for whether scaling is performed 
#' @param by_row A logical value for whether scaling is performed by row (TRUE)
#' of column (FALSE)
#' @param by_row determines whether scaling is performed on rows or columns
#' @return A dataframe compatible with ggplot2
#' @examples
#synthesise example data matrix for two treatment groups X and Y
#' mtx <- matrix(rnorm(mean = 7, n = 90), ncol = 10) 
#' row.names(mtx) <- LETTERS[1:9]
#' colnames(mtx) <- c(paste0("X", 1:5), paste0("Y", 1:5))
#' #add in some differential expression 
#' mtx[,6:10] <- mtx[,6:10] + 3
#' 
#' #Arrange and cluster data for heatmap plotting via ggplot2
#' hm_dfr <- heatmapArrange(mtx, 
#'                          cluster_row = TRUE,
#'                          cluster_column = FALSE,
#'                          scale = TRUE,
#'                          by_row = TRUE
#' ) 
#' 
#' library(ggplot2)
#' ggplot(hm_dfr, aes(x = col_var, y  = row_var, fill = value)) +
#'     geom_tile() +
#'     scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')
#' 
#' @export
heatmapArrange <- function(data_in, 
                           class_by = NULL, 
                           cluster_row = FALSE, 
                           cluster_column = FALSE, 
                           scale = TRUE, 
                           by_row = TRUE){
    #clustering : if all samples have the same value (i.e. 0), clustering fails 
    #so remove those which do
    homogeneous_rows <- sapply(1:nrow(data_in), function(x) sd(data_in[x,]))
    data_in <- data_in[homogeneous_rows != 0,]

    #SCALING
    #scale data by row or column (or neither)
    if(scale == TRUE){
        if(by_row == TRUE){
            dfr <- data.frame(t(scale(t(data_in))))
            colnames(dfr) <- colnames(data_in)
            } else {
            dfr <- data.frame(scale(data_in))
            colnames(dfr) <- colnames(data_in)
            }
        } else {
            dfr <- data.frame(data_in)
        }

    #ARRANGING
    #melt to ggplot2-able format
    dfr$row_var <- row.names(dfr)
    dfr_mlt <- reshape2::melt(dfr)

    #CLUSTERING
    #cluster by row
    if(cluster_row == TRUE){
        r_clst <- hclust(as.dist(1-cor(t(data_in), method = "pearson")), 
                         method = "complete", 
                         members = NULL)
        dfr_mlt$row_var <- factor(dfr_mlt$row_var, 
                                    levels = r_clst$labels[r_clst$order])
    }
    #cluster by column
    if(cluster_column == TRUE){
        c_clst <- hclust(as.dist(1-cor(data_in, method = "pearson")), 
                         method = "complete", 
                         members = NULL)
        dfr_mlt$variable <- factor(dfr_mlt$variable, 
                                   levels = c_clst$labels[c_clst$order])
    }
    colnames(dfr_mlt) <- c("row_var", "col_var", "value")

    print("scale_fill_gradient2 green-red colour reminder: scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850') ")
    dfr_mlt
}
