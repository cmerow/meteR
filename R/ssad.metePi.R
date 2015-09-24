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
