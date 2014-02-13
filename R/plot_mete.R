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

# source("~/R_functions/my_ecdf.R")

plot.mete <- function(x,type=c("abundance","power"),pers=c("rank","cumulative"),...) {
	pers <- match.arg(pers, choices=c("rank","cumulative"))
	
	X <- extract.mete(x,type)
	
	if(pers=="rank") {
		plot(X$data,...)
		points(X$rankFun,type="l")
	} else if(pers=="cumulative") {
		this.cdf <- my.ecdf(X$data)
		
		plot(this.cdf,...)
		curve(X$fun@p(x),add=TRUE)
	}
}
