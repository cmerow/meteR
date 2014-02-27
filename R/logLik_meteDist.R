##	function to compute log-likelihood of a meteDist object
logLik.meteDist <- function(x) {
	return(sum(x$fun@d(x$data,log=TRUE)))
}

##  function to compute log-likelihood z-score
logLikZ.meteDist <- function(x, nrep, return.sim=FALSE) {
	lik.obs <- logLik(x)
	lik.sim <- replicate(nrep, {
		new.dat <- x$fun@r(length(x$data))
		sum(x$fun@d(new.dat,log=TRUE))
	})
	
	if(return.sim) {
		return(lik.sim)
	} else {
		return((lik.obs-mean(lik.sim))/sd(lik.sim))
	}
}


##	function to compute deviance of a meteDist object
deviance.meteDist <- function(x) {
	mod.logLik <- logLik.meteDist(x)
	
	x.prob <- table(x$data)
	x.prob <- x.prob/sum(x.prob)
	perf.logLik <- sum(log(as.numeric(x.prob[as.character(x$data)])))
	
	return(-2 * (mod.logLik - perf.logLik))
}