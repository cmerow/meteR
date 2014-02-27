# source("~/R_functions/my_ecdf.R")

plot.mete <- function(x,type=c("abundance","power"),pers=c("rank","cumulative"),...) {
	pers <- match.arg(pers, choices=c("rank","cumulative"))
	
	X <- extract.mete(x,type)
	
	if(pers=="rank") {
		plot(X$data,...)
		points(X$rankFun,type="l")
	} else if(pers=="cumulative") {
		this.cdf <- my.ecdf(X$data)
		
		plot(this.cdf,...)
		curve(X$fun@p(x),add=TRUE)
	}
}
