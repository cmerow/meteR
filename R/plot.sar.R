## function to plot objects of class sar

plot.sar <- function(x, ...) {
    plot(x[, 'A'], x[, 'S'], xlab='Area', 
         ylab=sprintf('Number of %s', ifelse(attr(x, 'type') == 'ear', 'endemics', 'species')), 
         type=ifelse(attr(x, 'source')=='emperical', 'p', 'l'), ...)	
}