##	function to make Pi distribution
##	based on function for R(n,epsilon) in `makeMete'
metePi <- function(n, n0=sum(n), A, A0) {
	# n0 <- sum(n)
	
	thisSSF <- .makeSSF(n0, A, A0)
	# SSAD <- makeSSAD(n, thisSSF)
	
	out <- thisSSF
	out$data$n <- n
	class(out) <- 'metePi'
	return(out)
}

##	PMF for Pi
.mete.Pi <- function(n,la,n0) {
	1/.mete.Pi.Z(la, n0) * exp(-la*n)
}

##	normalization constant for Pi
.mete.Pi.Z <- function(la,n0) {
	if(la != 0) {
		return((1-exp(-la*(n0+1)))/(1-exp(-la)))
	} else {
		return(n0+1)
	}
}

##	constraint function for Pi lagrange multiplier
##	function: x = e^(-la)
##  vectorized over `x'
.pi.cons <- function(x,n0,A,A0) {
	lhs <- rep(NA,length(x))
	
	case1 <- x > 1 & (n0+1)*log(x) <= log(2e+64)
	case2 <- x > 1 & (n0+1)*log(x) > log(2e+64)
	case3 <- x < 1
	case4 <- x == 1
	case0 <- case1 | case3
	
	lhs[case0] <- x[case0]/(1-x[case0]) - ((n0 + 1)*x[case0]^(n0+1))/(1-x[case0]^(n0+1))
	lhs[case2] <- x[case2]/(1-x[case2]) + n0 + 1
	lhs[case4] <- n0/2
	
	return(lhs - n0*A/A0)
}


##	make *S*patial *S*tructure *F*unction, like the ESF slot in mete class
.makeSSF <- function(n0, A, A0, eq52=.useEq52(n0,A,A0)) {
	# browser()
	if(A/A0 == 0.5) {
		return(list(La=0, La.info='analytic solution',
                state.var=c(n0=n0,A=A,A0=A0)))
	} else if(eq52) {
		La <- -log((n0*A/A0)/(1+n0*A/A0))
		return(list(La=La, La.info='analytic solution',
                state.var=c(n0=n0,A=A,A0=A0)))
	} else {
		if(.pi.cons(1-.Machine$double.eps^0.45,n0,A,A0) > 0) {
			upper <- 1-.Machine$double.eps^0.45
		} else {
			upper <- 1.01*(n0*(A/A0-1)-1)/(n0*(A/A0-1))
		}
		
		sol <- uniroot(.pi.cons,c(0,upper),
                   n0, A, A0,
                   tol=.Machine$double.eps^0.75)
		La <- -log(sol$root)
		
		return(list(La=La,La.info=sol[-1],
                state.var=c(n0=n0,A=A,A0=A0)))
	}
}


## anything with n0 > 8000 and A/A0 < 0.5 use approx,
## otherwise follow this eq (validation in `check_eq52.R')
## vectorized over `n0'
.useEq52 <- function(n0,A,A0) {
	test <- exp(4.9 -1.36*log(n0) + 
                0.239*log(n0)^2 -0.0154*log(n0)^3)
	res <- A0/A >= test
	
	res[n0 > 2^16 & A/A0 < 0.5] <- TRUE
	res[A/A0 == 0.5] <- TRUE
	
	return(res)
}

