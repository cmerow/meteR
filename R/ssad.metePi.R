#' @title Species Spatial Abundance Distribution
#'
# @description 
#'
# @details
#' 
#' 
#' @param ssf 
#' 
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

ssad.metePi <- function(ssf) {
	x <- ssf$data$n
	
	if(is.null(x)) {
		X <- NULL
	} else {
		X <- sort(x, decreasing=TRUE)
	}
	
	this.eq <- function(n, log=FALSE) {
		out <- .mete.Pi(n, ssf$La, ssf$state.var['n0'])
		if(log) out <- log(out)
		
		return(out)
	}
	
	FUN <- distr::DiscreteDistribution(supp=0:ssf$state.var['n0'],
	                                   prob=this.eq(0:ssf$state.var['n0']))
	
	out <- list(type='ssad', data=X,
	            d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
	            state.var=ssf$state.var, La=ssf$La)
	
	class(out) <- c('ssad', 'meteDist')
	
	return(out)
}
