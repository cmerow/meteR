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

#library(nleqslv)

mete.lambda <- function(S0,N0,E0) {
	##	reasonable starting values
	init.la2 <- S0/(E0-N0)
	beta.guess <- 0.01
	init.beta <- nlm(function(b) {(b*log(1/b) - 1/2^8)^2},p=0.001)
	
	if(init.beta$code < 4) {	# was there some level of convergence?
		init.beta <- init.beta$estimate
	} else {
		init.beta <- beta.guess
	}
	
	init.la1 <- init.beta - init.la2
	
	##	the solution
	la.sol <- nleqslv(x=c(la1 = init.la1, la2 = init.la2), fn = la.syst2, S0=S0,N0=N0,E0=E0)
	
	return(list(lambda=c(la1=la.sol$x[1], la2=la.sol$x[2]), syst.vals=la.sol$fvec, converg=la.sol$termcd,
				mesage=la.sol$message, nFn.calc=la.sol$nfcnt, nJac.calc=la.sol$njcnt))
}


la.syst2 <- function(La,S0,N0,E0) {
	##	params
	b <- La[1] + La[2]
	s <- La[1] + E0*La[2]
	
	n <- 1:N0
	
	##	expressions
	g.bn <- exp(-b*n)
	g.sn <- exp(-s*n)
	
	univ.denom <- sum((g.bn - g.sn)/n)
	rhs.7.19.num <- sum(g.bn - g.sn)
	rhs.7.20.num <- sum(g.bn - E0*g.sn)
	
	##	the two functions to solve
	f <- rep(NA,2)
	
	f[1] <- rhs.7.19.num/univ.denom - N0/S0
	f[2] <- (1/La[2]) + rhs.7.20.num/univ.denom - E0/S0
	
	return(f)
}

#
#tryS0 <- 64
#tryN0 <- 2^4*tryS0
#tryE0 <- 2^10*tryN0
#mete.lambda(tryS0,tryN0,tryE0)
