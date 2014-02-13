#' @title Title of function
#'
#' @description
#' \code{function.name} what it does
#'
#' @details
#' how it works
#' etc.
#' 
#' @param arg description of arg blah blah blah blah
#' @param arg description of arg blah blah
#' @keywords manip
#' @export
#' 
#' @examples
#' #code to run
#' 
#' @return An object of class "mete"
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#  @note other junk to mention
#  @seealso - to provide pointers to other related topics
#  @references - references to scientific literature on this topic
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.

# library(distr)

##	function to take raw data (x) and create all needed data, as
##	well as all parameters and density functions
##
##	DATA:
##		species abundance
##		species energy
##		individual energy
makeMete <- function(otu,abund,power,min.e) {
	otu <- as.character(otu)
	
	if(!(class(abund) %in% c("numeric","integer","double"))) stop("argument `abund' must be numeric")
	if(!(class(power) %in% c("numeric","integer","double"))) stop("argument `power' must be numeric")
	
	if(missing(min.e)) min.e <- min(power)
	power <- power/min.e
	
	thisS0 <- length(unique(otu))
	thisN0 <- sum(abund)
	thisE0 <- sum(power)
	
	thisESF <- makeESF(s0=thisS0, n0=thisN0, e0=thisE0)
	thisESF$emin <- min.e
	SAD <- makeMeteDist(x=abund, otu=otu, esf=thisESF, type="abund")
	IPD <- makeMeteDist(x=power, otu=otu, esf=thisESF, type="power")
	ISP <- makeMeteDist(x=power, otu=otu, esf=thisESF, type="theta")
	
	return(list(ESF=thisESF, SAD=SAD, IPD=IPD, ISP=ISP))
}


##	function to compute lambdas, normalizing constant and eventually
##	R(n,epsilon), the `Ecosystem Structure Function' (ECF).
makeESF <- function(s0,n0,e0) {
	esf.par <- mete.lambda(s0,n0,e0)
	esf.par.info <- esf.par[-1]
	esf.par <- esf.par[[1]]
	
	names(esf.par) <- c("la1","la2")
	
	thisZ <- meteZ(esf.par["la1"],esf.par["la2"],s0,n0,e0)
	
	return(list(La=esf.par,La.info=esf.par.info,Z=thisZ,state.var=c(S0=s0,N0=n0,E0=e0)))
}


##	function to compute the various distributions of METE,
##	including tabulating the raw data for comparison
makeMeteDist <- function(x,otu,esf,type=c("abundance","power","theta")) {
	type <- match.arg(type,c("abundance","power","theta"))
	
	if(type=="abundance") {
		X <- sort(tapply(x,otu,sum),decreasing=TRUE)
		this.eq <- function(n) mete.phi(n=n,la1=esf$La[1],la2=esf$La[2],Z=esf$Z,
										S0=esf$state.var[1],N0=esf$state.var[2],
										E0=esf$state.var[3])
		
		FUN <- DiscreteDistribution(supp=1:esf$state.var["N0"],prob=this.eq(1:esf$state.var["N0"]))
		rankFun <- qphi2rank(FUN@q,esf$state.var["S0"])
	} else if(type=="power") {
		X <- sort(x,decreasing=TRUE)
		
		this.eq <- function(epsilon) mete.psi(epsilon, la1=esf$La[1], la2=esf$La[2], Z=esf$Z,
											  S0=esf$state.var[1], N0=esf$state.var[2], E0=esf$state.var[3])
		
		FUN <- AbscontDistribution(d=this.eq,low1=1,low=1,up=esf$state.var[3],up1=esf$state.var[3])
		# FUN <- AbscontDistribution(d=this.eq,p=thisP,low1=1,low=1,up=esf$state.var[3],up1=esf$state.var[3])
		######### somehow lost line where I get the p function through numeric integration `thisP'
		######### so above where p=NULL in AbsconDistribution is **stub**
		
		##	distr can't seem to handle this one...
		FUN@p <- Vectorize(function(epsilon,lower.tail=TRUE,log.p=FALSE) {
			out <- integrate(this.eq,lower=1,upper=epsilon)$value
			
			if(!lower.tail) out <- 1 - out
			
			if(log.p) out <- log(out)
			
			return(out)
		},vectorize.args="epsilon")
		
		rankFun <- qphi2rank(FUN@q,esf$state.var["N0"])  # won't be too accurate
	} else {
		X <- tapply(x,otu,c)
		
		this.eq <- function() "stub"	# function(epsilon,n)  mete.theta(epsilon = epsilon, n = n, la2 = esf$La[2])
		
		FUN <- "stub"	# AbscontDistribution(d=this.eq,low1=1,up=esf$state.var[3],up1=esf$state.var[3])
		
		rankFun <- NA
	}
	
	##	returning:	`type'		= what distribution is of
	##				`data'		= the raw data to be used
	##				`fun.eq'	= the density function
	##				`fun'		= Distribution class (functions for r, d, p, q)
	##				`rankFun'	= expected predict rank abund dist
	return(list(type=type, data=X, fun.eq=this.eq, fun=FUN, rankFun=rankFun))
}




