plot.meteRelat <- function(x, legend=TRUE, ...) {
    plot(x$obs, ...)
    lines(x$pred, col='red')
    
    if(legend) legend('topleft', c('METE prediction', 'Data'), col=c('red', 'black'), 
                      pch=c(NA, 1), lty=c(1, NA))
}