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

qphi2rank <- function(qfun,S0,...) {
	rel.abund <- qfun(seq(1-1/S0,0,by=-1/S0),...)
	rel.abund <- rel.abund[rel.abund > 0]
	return(rel.abund)
}

#ephi2rank <- function(comm,S0,...) {
#	xtab <- table(comm)
#	xval <- as.numeric(names(xtab))
#	xcum <- cumsum(as.numeric(xtab))/sum(xtab)
#
#xspp <- cut(xcum,breaks=seq(0,1,by=1/100),labels=1:100)
#head(xspp)
#
#xabund <- split(xval,xspp)
#}
#
#
#
#x <- qphi2rank(qlnorm,S0=40,meanlog=3,sdlog=2)
#x <- qphi2rank(qgeom,S0=400,prob=0.001)
#plot(x,log="y")
#
#
#
#curve(dgeom(x,0.9),from=0,to=40)