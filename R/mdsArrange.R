#'  mdsArrange
#' 
#' A function to take a matrix/dataframe of gene expression data and return a 
#' dataframe suitable for plotting as an MDS plot with ggplot2
#' 
#' @param d is a dataframe of gene expression data
#' @param isAlreadyScaled Logical value for whether d is already scaled 
#' @return A dataframe compatible with ggplot2
#' 
#' @examples
#' #Using the NKI dataset
#' library(breastCancerNKI)
#' library(Biobase)
#' data(nki)
#' mtx <- exprs(nki)
#' 
#' #Calculate PCs 1 & 2
#' mds_dfr <- mdsArrange(d = mtx, isAlreadyScaled = TRUE)
#' #Merge with phenotypic data for plotting
#' mds_mrg <- merge(mds_dfr, pData(nki), by = 0)
#' 
#' #Plot MDS and colour by estrogen receptor positivity
#' library(ggplot2)
#' ggplot(mds_mrg, aes(x = x, y  = y, colour = as.factor(er))) + 
#'     geom_point()
#' 
#' @export 
mdsArrange <- function(d, isAlreadyScaled = FALSE){
    #calculate PCs 1&2
    if(isAlreadyScaled == FALSE){
        dst <- dist(t(scale(d)))
        dN <- dimnames(d)[[2]]
        dst.m <- as.matrix(dst)
        dimnames(dst.m) <- list(dN, dN)
        dst.cmd <- cmdscale(dst.m, k=2)
    } else {
        dst <- dist(t(d))
        dN <- dimnames(d)[[2]]
        dst.m <- as.matrix(dst)
        dimnames(dst.m) <- list(dN, dN)
        dst.cmd <- cmdscale(dst.m, k=2)
    }
    #arrange as a single dataframe
    data.frame(ids = row.names(dst.cmd), x = dst.cmd[,1], y = dst.cmd[,2])
}
