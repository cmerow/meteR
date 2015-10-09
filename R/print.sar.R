#' @title print.meteESF
#'  
#' @description \code{print.meteESF} prints an object of class \code{meteESF}
#'
#' @details
#' See Examples
#' 
#' @param x an object of class \code{meteESF}
#' @param ... arguments to be passed
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(anbo)
#' anbo.sar <- meteSAR(anbo.new$spp, anbo.new$count, anbo.new$row, anbo.new$col, Amin=1, A0=16)
#' print(anbo.sar)
#' anbo.sar # alternatively
#' 
#' @return Returns the object silently
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow

print.sar <- function(x) {
	cat(sprintf('%s %s area relationship ranging from \n', attr(x, 'source'), 
	            ifelse(attr(x, 'type')=='sar', 'species', 'endemics')))
	
	cat(sprintf('A: [%s, %s] \nS: [%s, %s] \n', min(x[, 'A']), max(x[, 'A']), min(x[, 'S']), max(x[, 'S'])))
	
	invisible(x)	
}