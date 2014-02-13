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

residuals.meteDist <- function(X, type=c("rank","cumulative"),relative=TRUE,log=FALSE) {
	type <- match.arg(type,choices=c("rank","cumulative"))
	
	if(type=="rank") {
		obs <- X$data
		pred <- X$rankFun
	} else if(type=="cumulative") {
		obs <- my.ecdf(X$data)
		pred <- X$fun@p(obs[,1],log.p=log)
		if(log) obs <- log(obs[,2])
		else obs <- obs[,2]
	}
	
	out <- obs - pred
	if(relative) out <- out/abs(pred)
	
	return(out)
} 