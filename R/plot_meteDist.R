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

##	plot an object of meteDist (i.e. make rank distribution plot and 
##	show theoretical prediction).

plot.meteDist <- function(x,ptype=c("cdf","rad"),th.col="red",lower.tail=TRUE,...) {
	ptype <- match.arg(ptype,c("cdf","rad"))
	
	if(ptype=="cdf") {
		plot(my.ecdf(x$data,!lower.tail),...)
		this.curve <- x$fun@p
		
		if(x$type %in% c("richness","abundance")) {
			this.supp <- 1:sum(x$data)
			points(this.supp,this.curve(this.supp,lower.tail=lower.tail),type="l",col=th.col)
		} else {
			curve(this.curve(x,lower.tail=lower.tail),add=TRUE,col=th.col)
		}
	} else {
		if('ylim' %in% names(list(...))) {
			plot(sort(x$data,decreasing=TRUE),...)
		} else {
			plot(sort(x$data,decreasing=TRUE), ylim=range(x$data, x$rankFun), ...)
		}
		points(x$rankFun,type="l",col=th.col)
	}
}
