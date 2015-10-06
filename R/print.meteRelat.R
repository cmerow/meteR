## print method for objects of class `meteRelat'

print.meteRelat <- function(x) {
    cat(switch(attr(x$pred, 'type'), 
               'sar' = 'Species area relationship',
               'ear' = 'Endemcis area relationship'),
        sprintf('predicted using %s', ifelse(is.null(x$obs),
                                             'state variables only',
                                             'raw data')),
        '\n')

    print(x$pred)
    if(!is.null(x$obs)) print(x$obs)
    
    invisible(x)
}