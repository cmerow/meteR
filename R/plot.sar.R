## function to plot objects of class sar

plot.sar <- function(x, add=FALSE, ...) {
    if(add) {
    	points(x[, 'A'], x[, 'S'], type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    } else {
        plot(x[, 'A'], x[, 'S'], xlab='Area', 
             ylab=sprintf('Number of %s', ifelse(attr(x, 'type') == 'ear', 'endemics', 'species')), 
             type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    }
}
