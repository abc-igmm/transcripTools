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
#' #synthesise example data matrix
#' mtx <- matrix(rnorm(mean = 7, n = 90), ncol = 10) 
#' row.names(mtx) <- LETTERS[1:9]
#' colnames(mtx) <- c(paste0("X", 1:5), paste0("Y", 6:10))
#'
#' hm_dfr <- mdsArrange(d = mtx, isAlreadyScaled = FALSE)
#'
#' #library(ggplot2)
#' #ggplot(hm_dfr, aes(x = x, y  = y, fill = )) 
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
