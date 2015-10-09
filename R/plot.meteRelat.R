#================================================================
#' @title plot the ???
#'
# @description 
#'
# @details
#' 
#' 
#' @param x 
#' @param legend Logical 
#' @param ... arguments to be passed to plot
#' 
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))

# @return 
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

plot.meteRelat <- function(x, legend=TRUE, ...) {
    plot(x$obs, ...)
    lines(x$pred, col='red')
    
    if(legend) legend('topleft', c('METE prediction', 'Data'), col=c('red', 'black'), 
                      pch=c(NA, 1), lty=c(1, NA))
}