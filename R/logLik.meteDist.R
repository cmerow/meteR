#' @title Compute log-likelihood of a meteDist object
#'
#' @description
#' \code{logLik.meteDist} computes log-likelihood of a meteDist object
#'
#' @details
#' See Examples.
#' 
#' @param object a \code{meteDist} object
#' @param ... arguments to be passed
# @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' ipd1=ipd.meteESF(esf1)
#' logLik(ipd1)
#' 
#' @return object of class \code{logLik}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

logLik.meteDist <- function(object,...) {
	lik <- sum(object$d(object$data,log=TRUE))
	attr(lik, 'df') <- 1
	class(lik) <- 'logLik'
	return(lik)
}

#==============================================================================
#' @title Compute log-likelihood z-score
#'
#' @description
#' \code{logLikZ.meteDist} computes a log-likelihood z-score
#'
#' @details
#' \code{logLikZ.meteDist} simulates from a fitted METE distribution (e.g. a species abundance distribution or individual power distribution) and calculates the likelihood of these simulated data sets. The distribution of these values is compared against the likelihood of the data to obtain a z-score. 
#' 
#' @param x a \code{meteDist} object
#' @param nrep number of simulations from the fitted METE distribution 
#' @param return.sim logical; return the simulated liklihood values
# @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' ipd1=ipd.meteESF(esf1)
#' llz=logLikZ.meteDist(ipd1, nrep=100, return.sim=TRUE)
#' plot(density(llz$sim),xlim=range(c(llz$sim,llz$obs)),
#'      xlab='log(likelihood)',col='red')
#' abline(v=llz$obs,lty=2)
#' legend('top',legend=c('data','simulated'),col=c('black','red'),
#'       lty=c(1,1),bty='n') 
#' 
#' @return list with elements
#' \describe{
#'    \item{z}{The z-score}
#'    \item{obs}{log-likelihood}
#'    \item{sim}{\code{nrep} Simulated values}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso mseZ.meteDist
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.
logLikZ.meteDist <- function(x, nrep, return.sim=FALSE) {
	lik.obs <- logLik(x)
	lik.sim <- replicate(nrep, {
		                  new.dat <- x$r(length(x$data))
		                  sum(x$d(new.dat, log=TRUE))
	})
	
	if(return.sim) {
		return(list(z=(lik.obs-mean(lik.sim))/sd(lik.sim), 
		            obs=lik.obs,
		            sim=lik.sim))
	} else {
		return(list(z=(lik.obs-mean(lik.sim))/sd(lik.sim)))
	}
}
