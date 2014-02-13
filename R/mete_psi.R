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

mete.psi <- function(e,la1,la2,Z,S0,N0,E0) {
	gamma <- la1 + e*la2
	
	t1 <- S0/(N0*Z)
	t2 <- exp(-gamma)/((1-exp(-gamma))^2)
	t3 <- exp(-gamma*N0)/(1-exp(-gamma))
	t4 <- N0 + exp(-gamma)/(1-exp(-gamma))
	
	return(t1*(t2 - t3*t4))
}