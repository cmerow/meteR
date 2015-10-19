#' @title Species Spatial Abundance Distribution
#'
# @description 
#'
# @details
#' 
#' 
#' @param ssf An objects of class meteSSF; i.e. the spatial structure function \eqn{\Pi(n)}
#' 
#' @export
#' 
#' @examples
#' data(anbo)
#' samp2mat <- function(site,spp,abund) {
#'     y <- tapply(abund, list(site, spp), sum)
#'     y[is.na(y)] <- 0
#'     return(y)
#' }
#' anbo.mat <- samp2mat(paste(anbo[, 1], anbo[,2]), anbo$spp, anbo$count)
#' anbo.new <- data.frame(t(sapply(strsplit(rownames(anbo.mat), ' ', fixed=TRUE), as.numeric)), spp = rep(colnames(anbo.mat), each=nrow(anbo.mat)), count = as.vector(anbo.mat))
#' colnames(anbo.new)[1:2] <- colnames(anbo)[1:2]
#' ## anbo.new now has 0 abundance where needed
#' pi1 <- meteSSF(anbo.new$count[anbo.new$spp=='crcr'], A=1, A0=16)
#' plot(ssad(pi1))

# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

ssad <- function(x) {
  UseMethod('ssad')
}

#' @rdname ssad
# @method ssad meteSSF
# @S3method ssad meteSSF
#' @export 

ssad.meteSSF <- function(ssf) {
	x <- ssf$data$n
	
	if(is.null(x)) {
		X <- NULL
	} else {
		X <- sort(x, decreasing=TRUE)
	}
	
	this.eq <- function(n, log=FALSE) {
		out <- metePi(n, ssf$La, ssf$state.var['n0'])
		if(log) out <- log(out)
		
		return(out)
	}
	
	FUN <- distr::DiscreteDistribution(supp=0:ssf$state.var['n0'],
	                                   prob=this.eq(0:ssf$state.var['n0']))
	
	out <- list(type='ssad', data=X,
	            d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
	            state.var=ssf$state.var, La=ssf$La)
	
	class(out) <- c('ssad', 'meteDist')
	
	return(out)
}
