#================================================================
#' @title Plot the species abundance distribution (SAR), i.e. objects of class sar
#'
# @description
#'
# @details
#' 
#' 
#' @param x an object of class SAR made with 
#' 
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))

#' @return NULL
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

plot.sar <- function(x, add=FALSE, ...) {
    if(add) {
    	points(x[, 'A'], x[, 'S'], type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    } else {
        plot(x[, 'A'], x[, 'S'], xlab='Area', 
             ylab=sprintf('Number of %s', ifelse(attr(x, 'type') == 'ear', 'endemics', 'species')), 
             type=ifelse(attr(x, 'source')=='empirical', 'p', 'l'), ...)
    }
}
