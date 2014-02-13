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

##	print and summary methods for METE object

print.METE <- function(x) {
	cat("METE object with state variables:\n")
	print(x[[1]]$state.var)
	
	cat("\n");cat("with Lagrange multipliers:\n")
	print(x[[1]]$La)
	
	cat("\n");cat("predicting the distributions:\n")
	cat(names(x)[-1]);cat("\n")
}

# summary.METE <- function(x) {
	
# }