#' @title Plot to compare METE predictions to data 
#'
#' @description
#' \code{plot.meteDist} plots both the theoretical prediction and data for a \code{meteDist} object using either a rank or cumulative distribution plot
#'
#' @details
#' \code{plot.meteDist} automatically extracts the prediction and data from the \code{meteDist} objects. Additional plotting arguments can be passed from \code{par}.
#' 
#' @param x a \code{meteDist} object
#' @param ptype either "cdf" or "rad"
#' @param th.col line color of theoretical prediction 
#' @param lower.tail logical; choose TRUE to highlight differences between data and theory at low abundance; choose FALSE to highlight differences at high abundance.
#' @param add.legend logical; add a legend
#' @param ... arguments to be passed to \code{plot}
#' @param add.line add the curve for a fitted model to the existing plot
# @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' ipd1=ipd.meteESF(esf1)
#' plot(ipd1)
#' plot(ipd1, ptype='rad')
#' 
# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

#   if(length(x$data)==0 | add.line){ # diffent routine if no data...
#     if(ptype=="cdf") {
#       this.curve <- x$p
#       # start with an empty plot
#       if(!add.line){do.call(plot, c(list(this.supp,y=this.curve(this.supp,lower.tail=lower.tail), plot.par,col='white')))}
#       if(x$type %in% c("gsd", "sad")) {
#         this.supp <- 1:sum(x$data)
#         points(this.supp, this.curve(this.supp,lower.tail=lower.tail),
#                type="l", col=th.col)
#       } else {
#         curve(this.curve(x,lower.tail=lower.tail), add=TRUE, col=th.col)
#       }
#     }
#     else {
#       if(!('ylim' %in% names(plot.par))) {
#         plot.par$ylim <- x$q(c(1-1/x$state.var['S0'],1/x$state.var['S0']))
#       }
#       if(!add.line){do.call(plot, c(list(x=sort(x$data,decreasing=TRUE)), plot.par))}
#       points(meteDist2Rank(x), type="l", col=th.col)
#     }
#       
# 
#   } else { # for plotting with data


plot.meteDist <- function(x, ptype=c("cdf","rad"), th.col="red", 
                          lower.tail=TRUE, add.legend=TRUE, add.line=FALSE, ...) {
	ptype <- match.arg(ptype,c("cdf", "rad"))
	
  plot.par <- list(...)
  
  if(!('ylab' %in% names(plot.par))) {
    ylab <- ifelse(ptype=='cdf', 'Cumulative probability', '%s')
    ylab <- sprintf(ylab, switch(x$type,
                                 'sad' = 'Abundance',
                                 'ipd' = 'Metabolic rate'))
    plot.par$ylab <- ylab
  }
  
	if(!('xlab' %in% names(plot.par))) {
	  xlab <- ifelse(ptype=='cdf', '%s', 'Rank')
	  xlab <- sprintf(xlab, switch(x$type,
	                               'sad' = 'Abundance',
	                               'ipd' = 'Metabolic rate'))
    plot.par$xlab <- xlab
	}
  
	if(ptype=="cdf") {
	  this.curve <- x$p
    ## if no data, don't plot it, just plot the curve
	  if(is.null(x$data)) {
	    xmax <- max(x$data)
	    X <- .ecdf(x$data, !lower.tail)
	    plot.par$type <- 'n'
	  } else {
	    xmax <- ifelse(is.finite(max(plot.par$xlim)), 
	                   max(plot.par$xlim), 
	                   x$state.var['N0']/x$state.var['S0'])
	    X <- cbind(c(1, floor(xmax)), this.curve(c(1, floor(xmax))))
	  }
	  
	  do.call(plot, c(list(x=X, plot.par)))
	  
	  if(x$type %in% c("gsd", "sad")) {
	    this.supp <- 1:xmax
	    points(this.supp, this.curve(this.supp,lower.tail=lower.tail),
	           type="l", col=th.col)
	  } else {
	    curve(this.curve(x,lower.tail=lower.tail), add=TRUE, col=th.col)
	  }
	} else {
	  ## if no data, don't plot it, just plot the rank fun
	  if(is.null(x$data)) {
	    X <- x$rankFun
	    plot.par$type <- 'n'
	  } else {
	    X <- x$data
	  }
	  
	  ## if ylim not already specified make sure both data and theory fit
	  if(!('ylim' %in% names(plot.par))) {
	    plot.par$ylim <- range(X, x$rankFun)
	  }
	  
	  ## do plotting
	  do.call(plot, list(x=X, plot.par))
	  points(meteDist2Rank(x), type="l", col=th.col)
	}

  if(add.legend) legend('right', legend=c('data', 'METE'), col=c('black', 'red'),
                        lty=c(NA, 1), pch=c(21, NA), bty='n') 
}

