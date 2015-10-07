## generic
sipd <- function(x, ...) {
	UseMethod('sipd')
}

## function to extract species level individual power distribution 
## from ESF object and return object inheriting from meteDist
sipd.meteESF <- function(esf, sppID, n) {
    if(is.na(esf$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
    
    if(!missing(n)) { # n provided, only return theoretical (cause how deal with multiple spp with same n?)
        if(!missing(sppID)) warning('`sppID` ignored if `n` provided; to explore individual species use `sppID` instead')
        X <- NULL # data null if n provided
    } else if(!is.null(esf$data)) { # data provided
        if(!missing(sppID)) {
            spp <- as.factor(esf$data$s)
            if(is.numeric(sppID)) sppID <- levels(spp)[sppID]
            
            n <- sum(esf$data$n[esf$data$s == sppID])
            x <- esf$data$e
            if(is.null(x)) { # no energy data
                X <- NULL
            } else {
                X <- sort(x, decreasing=TRUE)
            }
        } else {
            
            stop('must provide either `sppID` or `n`')
        }
    } else {
        stop('must provide `n` if `esf` has no empirical data')
    }
    
    this.eq <- function(epsilon, log=FALSE) {
        out <- meteTheta(epsilon, n=n, la2=esf$La[2])
        if(log) out <- log(out)
        return(out)
    }
    
    this.p.eq <- function(epsilon, log=FALSE, lower.tail=TRUE) {
        out <- .meteThetaCum(epsilon, n, esf$La[2], log)
        if(!lower.tail) out <- 1 - out
        
        return(out)
    }
    
    FUN <- distr::AbscontDistribution(d=this.eq, p=this.p.eq,
                                      low1=1, low=1, up=esf$state.var[3], up1=esf$state.var[3])
    
    out <- list(type='sipd', data=X,
                d=this.eq, p=this.p.eq, q=FUN@q, r=FUN@r,
                state.var=esf$state.var, La=esf$La)
    out$state.var['n'] <- n
    class(out) <- c('ipd', 'meteDist')
    
    return(out)
}

#' @title Title of function
#'
#' @description
#' \code{function.name} what it does
#'
#' @details
#' how it works
#' etc.
#' 
#' @param arg description of arg
#' @param arg description of arg
#' @keywords manip
#' @export
#' 
#' @examples
#' #code to run
#' 
#' @return - the type of object that the function returns
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
                                        #  @note other junk to mention
                                        #  @seealso - to provide pointers to other related topics
                                        #  @references - references to scientific literature on this topic
                                        #  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
                                        #  @family - a family name. All functions that have the same family tag will be linked in the documentation.

##	Intra-specific metabolic rate distribution

meteTheta <- function(e, n, la2) {
    la2 * n * exp(-la2 * n * (e-1))
}

## integral of meteTheta
.meteThetaCum <- function(e, n, la2, log=FALSE) {
    out <- 1 - exp(-la2 * n * (e - 1))
    if(log) out <- log(out)
    
    return(out)
}
