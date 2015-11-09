#' @title print.sar
#'  
#' @description S3 method for class \code{sar}
#'
#' @details
#' See Examples
#' 
#' @param x an object of class \code{sar}
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

print.damuth <- function(x) {
  cat('Abundance metabolic rate relationship ranging from \n')
  
  cat(sprintf('n: [%s, %s] \ne: [%s, %s] \n', 
              min(x[['n']], na.rm=TRUE), 
              max(x[['n']], na.rm=TRUE), 
              min(x[['e']], na.rm=TRUE), 
              max(x[['e']], na.rm=TRUE)))
  
  invisible(x)	
}