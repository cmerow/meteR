
#' @title Relationship between mean metabolic rate (\eqn{\bar{\epsilon}}) and abundance
#'  
#' @description \code{ebar} calculates the relationship between average metabolic rate of a species and that species' abundance. Also known as the Damuth relationship
#' @details
#' See examples.
#' 
#' @param esf an object of class meteESF. 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' damuth <- ebar(esf1)
#' 
#' @return An object of class \code{meteRelaT}. The object contains a list with the following elements.
#' \describe{
#'   \item{\code{pred}}{predicted relationship}
#'   \item{\code{obs}}{observed relationship}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteDist, sad.meteESF, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.


ebar <- function(esf) {
  if(is.na(esf$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
  
  x <- esf$data$e
  
  if(is.null(x)) {
    X <- NULL
  } else {
    X <- aggregate(list(n=esf$data$n, e=esf$data$e), esf$data$s, sum)
    X$e <- X$e/X$n
  }
  
  thr <- data.frame(n=1:min(esf$state.var['N0'], X$n), e=1 + 1/(1:min(esf$state.var['N0'], X$n) * la2))
  
  out <- list(obs=X, pred=thr)
  class(out) <- 'meteRelat'
  
  return(out)
}