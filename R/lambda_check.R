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

##	right hand side of eq 7.19 from John's book
rhs.7.19 <- function(la1,la2,N0,E0) {
	b <- la1 + la2
	s <- la1 + E0*la2
	
	n <- 1:N0
	
	sum(exp(-b*n) - exp(-s*n))/sum((exp(-b*n) - exp(-s*n))/n)
}

##	right hand side of eq 7.20 from John's book
rhs.7.20 <- function(la1,la2,N0,E0) {
	b <- la1 + la2
	s <- la1 + E0*la2
	
	n <- 1:N0
	
	(1/la2) + sum(exp(-b*n) - E0*exp(-s*n))/sum((exp(-b*n) - exp(-s*n))/n)
}

##  NOTE: not the most efficient to keep these sepparate, but i do it for readability

##	system to solve, corresponds to setting lhs of eqs. 7.19 and 7.20 equal to zero
la.syst <- function(La,S0,N0,E0) {
	f <- rep(NA,2)
	
	f[1] <- rhs.7.19(La[1],La[2],N0,E0) - N0/S0
	f[2] <- rhs.7.20(La[1],La[2],N0,E0) - E0/S0
	
	return(f)
}


##	some state variable values from table 7.2 of John's book
tryS0 <- 64
tryN0 <- 2^4*tryS0
tryE0 <- 2^10*tryN0

##	the solution
la.sol <- nleqslv(x=c(la1=-0.0016,la2=3.82e-06),fn=la.syst,S0=tryS0,N0=tryN0,E0=tryE0)


##########  sanity check  ##########
tryla1 <- la.sol$x[1]
tryla2 <- la.sol$x[2]

tryN0/tryS0;rhs.7.19(tryla1,tryla2,tryN0,tryE0)
tryE0/tryS0;rhs.7.20(tryla1,tryla2,tryN0,tryE0)


