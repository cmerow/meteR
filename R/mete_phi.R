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

##	function to return the species abundance distribution (phi)
##	predicted by METE; vectorized in n

mete.phi <- function(n,la1,la2,Z,S0,N0,E0) {
	if(missing(Z)) Z <- meteZ(la1,la2,S0,N0,E0)
	
	beta <- la1 + la2
	sigma <- la1 + E0*la2
	
	return((exp(-beta*n) - exp(-sigma*n))/(la2*Z*n))
}


##	function to return Z, the normalizing constant of R as well as
##	the simplifying parameters beta and sigma

meteZ <- function(la1,la2,S0,N0,E0) {
	beta <- la1 + la2
	sigma <- la1 + E0*la2
	
	t1 <- S0/(la2*N0)
	t2 <- (exp(-beta) - exp(-beta*(N0+1)))/(1-exp(-beta))
	t3 <- (exp(-sigma) - exp(-sigma*(N0+1)))/(1-exp(-sigma))
	
	Z <- t1*(t2 - t3)
	
	return(Z)
}