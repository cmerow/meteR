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

##	function to compute log-likelihood of a meteDist object
logLik.meteDist <- function(x) {
	return(sum(x$fun@d(x$data,log=TRUE)))
}

logLikZ.meteDist <- function(x,nrep) {
	lik.obs <- logLik(x)
	lik.sim <- replicate(nrep, {
		new.dat <- x$fun@r(length(x$data))
		sum(x$fun@d(new.dat,log=TRUE))
	})
	
	(lik.obs-mean(lik.sim))/sd(lik.sim)
}


##	function to compute deviance of a meteDist object
deviance.meteDist <- function(x) {
	mod.logLik <- logLik.meteDist(x)
	
	x.prob <- table(x$data)
	x.prob <- x.prob/sum(x.prob)
	perf.logLik <- sum(log(as.numeric(x.prob[as.character(x$data)])))
	
	return(-2 * (mod.logLik - perf.logLik))
}