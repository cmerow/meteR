#' @title METE species abundance distribution
#'
#' @description
#' \code{sad.mete} returns the species abundance distribution
#' predicted by METE (Phi(n))
#'
#' @details
#' See Examples.
#' 
#' @param esf an object of class mete. 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure 
#' function
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' sad=sad.meteESF(esf1)
#' 
#' @return An object of class \code{meteDist} which inherits from objects of class \code{sad}. The object contains a list with the following elements.
#' \describe{
#'    \item{\code{type}}{'sad'}
#'    \item{\code{data}}{X}
#'    \item{\code{d}}{this.eq}
#'    \item{\code{p}}{FUNp}
#'    \item{\code{q}}{FUNq}
#'    \item{\code{r}}{FUNr}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @note
#' @seealso metePhi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


sad.meteESF <- function(esf) {
    x <- esf$data$n
    otu <- esf$data$s
    
    if(is.null(x)) {
      X <- NULL
    } else {
      X <- sort(tapply(x,otu,sum),decreasing=TRUE)
    }
    
    
    this.eq <- function(n, log=FALSE) {
      out <- metePhi(n=n,la1=esf$La[1],la2=esf$La[2],Z=esf$Z,
                     S0=esf$state.var[1],N0=esf$state.var[2],
                     E0=esf$state.var[3])
      if(log) out <- log(out)
      return(out)
    }
    
    FUN <- distr::DiscreteDistribution(supp=1:esf$state.var["N0"], 
                                prob=this.eq(1:esf$state.var["N0"]))
    
    # rankFun <- qphi2rank(FUN@q,esf$state.var["S0"])
    
    out <- list(type='sad', data=X, 
                d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
                state.var=esf$state.var, La=esf$La)
    class(out) <- c('sad', 'meteDist')
    
    return(out)
}



#==============================================================================
#' @title Equation of the METE species abundance distribution
#'
#' @description
#' \code{metePhi} returns the species abundance distribution 
#' (Phi(n)) predicted by METE; vectorized in n
#'
#' @details
#' See Examples
#' 
#' @param n the value (number of individuals) at which to calculate \deqn{\Phi}
#' @param la1,la2 Lagrange multipliers
#' @param Z partition function
#' @param S0 Total number of species
#' @param N0 Total number of individuals
#' @param E0 Total metabolic rate
#' 
#' @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' metePhi(min(arth$mass^(.75)),
#'        esf1$La[1],esf1$La[2],
#'        esf1$Z,esf1$state.var['S0'],
#'        esf1$state.var['N0'],
#'        esf1$state.var['E0'])
#' 
#' @return numeric
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
#' @seealso \code{sad.mete}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.


metePhi <- function(n, la1, la2, Z, S0, N0, E0) {
    if(missing(Z)) Z <- .meteZ(la1, la2, S0, N0, E0)
    
    beta <- la1 + la2
    sigma <- la1 + E0*la2
    
    return((exp(-beta*n) - exp(-sigma*n))/(la2*Z*n))
}
