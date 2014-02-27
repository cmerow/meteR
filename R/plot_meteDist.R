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
