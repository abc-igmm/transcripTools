#' heatmapArrange 
#' 
#' A function to take a matrix/dataframe of gene expression data and return a 
#' dataframe suitable for plotting as a heatmap with ggplot2
#' 
#' @param data_in is a dataframe of gene expression data
#' @param class_by is a vector of strings/regular expression that are used to
#' parse column names on the fly to produce a class column in the dataframe 
#' output. These must reflect some common element within the input data's sample
#' names that can be used to define classes. See examples...
#' @param cluster_row A logical value for whether row clustering is performed
#' @param cluster_column A logical value for whether column clustering is 
#' performed
#' @param scale A logical value for whether scaling is performed 
#' @param by_row A logical value for whether scaling is performed by row (TRUE)
#' of column (FALSE)
#' @param by_row determines whether scaling is performed on rows or columns
#' @return A dataframe compatible with ggplot2
#' @examples
#' #synthesise example data matrix
#' mtx <- matrix(rnorm(mean = 7, n = 90), ncol = 10) 
#' row.names(mtx) <- LETTERS[1:9]
#' colnames(mtx) <- c(paste0("X", 1:5), paste0("Y", 6:10))
#'
#' hm_dfr <- heatmapArrange(mtx, 
#'                          class_by = c("X", "Y"),
#'                          cluster_row = TRUE,
#'                          cluster_column = TRUE,
#'                          scale = TRUE,
#'                          by_row = TRUE
#' ) 
#'
#' #library(ggplot2)
#' #ggplot(hm_dfr, aes(x = x, y  = y, fill = )) 
#'
#' @export
heatmapArrange <- function(data_in, 
                           class_by = "", 
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
    dfr$row_value <- row.names(dfr)
    dfr_mlt <- reshape2::melt(dfr)

    #CLUSTERING
    #cluster by row
    if(cluster_row == TRUE){
        r_clst <- hclust(as.dist(1-cor(t(data_in), method = "pearson")), 
                         method = "complete", 
                         members = NULL)
        dfr_mlt$row_value <- factor(dfr_mlt$row_value, 
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

    #SAMPLE CLASSIFYING
    #assign a sample type if wanted
    if(!identical(class_by, "")){
        for(x in 1:length(class_by)){
            class_vec <- c()
            for(i in class_by[[x]]){
                class_vec[grep(i, droplevels(dfr_mlt$variable))] <- i
            }
            dfr_mlt[[names(class_by)[x]]] <- class_vec
        }
    }
    print("scale_fill_gradient2 green-red colour reminder: 
          scale_fill_gradient2(high = '#d73027', mid = 'black', low = '#1a9850')
          ")
    dfr_mlt
}
