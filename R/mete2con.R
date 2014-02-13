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

###########  functions for calculating METE SAD under 2nd resource constraint  ###########

mete2con.lik <- function(zeta,data) {
	this.fun <- Vectorize(function(x) -sum(log(mete.phi2con(data,sum(data),x))))
	
	out <- this.fun(zeta)
	out[!is.finite(out)] <- NA
	
	return(out)
}

mete2con.mle <- function(data) {
	nlm(mete2con.lik,p=1,data=data)
}

##	normalized SAD
mete.phi2con <- function(n,N0,zeta) {
	kern <- function(x) exp(-zeta*x)/(x^2)
	Z <- sum(kern(1:N0))
	
	return(kern(n)/Z)
}

##	normalized CDF
mete.phi2con.cdf <- function(N0,zeta) {
	return(list(x = 1:N0, y = cumsum(mete.phi2con(1:N0,N0,zeta))))
}
