
#' @title Generic method to obtain the species-level individual power distribution (SIPD)
#'
#' @description Extract species level individual power distribution 
#' from ESF object and return object inheriting from meteDist. This distribution (\eqn{\Theta}) describes the distribution of metabolic rates across the individuals of a species with n individuls
#'
#' @details 
#' \deqn{
#'    \Theta( \epsilon \mid n, S_{0}, N_{0}, E_{0} ) \approx \lambda_{2} n e^{- \lambda_[2] n (\epsilon -1)}
#' }
#' 
#' 
#' @param esf an object of class meteESF (i.e. the fitted distribution \eqn{R(n,e)})
#' @param sppID
#' @param n integer. number of individuals in the species
#'
#' @rdname sipd
#' @export
#' 
#' @examples
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' sipd1 <- sipd(arth.esf, sppID=5)
# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

sipd <- function(x, ...) {
	UseMethod('sipd')
}

# @return \code{NULL}
#' @rdname sipd
#' @method sipd meteESF
#' @S3method sipd meteESF
#' 
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

#================================================================
#' @title Intra-specific metabolic rate distribution
#'
#' @description Distributi􏴜on of metabolic rates over individuals within a species of abundance n0 (Theta)
 
# $\Theta (\epsilon | n0,S0,N0,E0) = \dfrac{R}{\phi}$#'
# @details  
#' 
#' @param e 
#' @param n
#' @param la2 
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))

# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

meteTheta <- function(e, n, la2) {
    la2 * n * exp(-la2 * n * (e-1))
}

## integral of meteTheta
.meteThetaCum <- function(e, n, la2, log=FALSE) {
    out <- 1 - exp(-la2 * n * (e - 1))
    if(log) out <- log(out)
    
    return(out)
}
