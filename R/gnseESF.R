#' @title gsneESF
#'  
#' @description \code{gsneESF} Calculates the ``ecosystem structure
#' function'' using a higher taxon (e.g. genera) as an additional constraint
#'
#' @details
#' Use follows that of \code{meteESF}
#' 
#' @param gen A vector of genera names
#' @param spp A vector of species names 
#' @param abund A vector of abundances 
#' @param power A vector of metabolic rates 
#' @param G0 Total number of genera
#' @param S0 Total number of species
#' @param N0 Total number of individuals
#' @param E0 Total metabolic rate; defaults to N0*1e6 if not specified or 
#'                  calculated from \code{power} to allow one to fit models that 
#'                  do not depend on metabolic rates
#' @param minE Minimum possible metabolic rate
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
# @examples
#' 
#' @return An object of class \code{meteESF} with elements
#' \describe{
#'   \item{\code{data}}{The data used to construct the ESF}
#'   \item{\code{emin}}{The minimum metabolic rate used to rescale metabolic rates}
#'   \item{\code{La}}{Vector of Lagrange multipliers}
#'   \item{\code{La.info}}{Termination information from optimization procedure}
#'   \item{\code{state.var}}{State variables used to constrain entropy maximization}
#'   \item{\code{Z}}{Normalization constant for ESF}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteESF
# @references 

#' @useDynLib meteR
#' @importFrom Rcpp sourceCpp

gsneESF <- function(gen, spp, abund, power, G0 = NULL, S0 = NULL, N0 = NULL, E0 = NULL, minE) {
  ## case where real data given
  if(!missing(gen) & !missing(spp) & !missing(abund)) {
    gen <- as.character(gen)
    spp <- as.character(spp)
    
    ## set minE to min of power unless power is missing, then set minE to 1
    if(missing(minE) & !missing(power)) {
      minE <- min(power)
    } else {
      minE <- 1
    }
    
    ## simple stat variable calculations
    G0 <- length(unique(gen))
    S0 <- length(unique(spp))
    N0 <- sum(abund)
    
    ## calculating E0 depends on if power is missing
    if(missing(power)) {
      ## E0 could be specified with missing power, in which case use it, otherwise don't
      if(missing(E0)) {
        E0 <- N0*10^6
        e.given <- FALSE
      } else {
        E0 <- E0 / minE
        e.given <- TRUE
      }
    } else {
      power <- power/minE
      E0 <- sum(power*abund)
      e.given <- TRUE
    }
    
  ## case where only state variables given (E0 may or may not be given, set useful default if missing)
  } else {
    if(missing(G0) | missing(S0) | missing(N0)) 
      stop('must provide either data or state variables')
    if(missing(E0)) {
      E0 <- N0 * 10^6
      e.given <- FALSE
    } else {
      E0 <- E0 / ifelse(missing(minE), 1, minE)
      e.given <- TRUE
    }
  }
  
  thisESF <- .makegsneESF(G0, S0, N0, E0)
  thisESF$data <- NULL
  thisESF$emin <- ifelse(missing(minE), 1, minE)
  if(!e.given) thisESF$state.var['E0'] <- NA
  
  return(thisESF)
}


## makes ecosystem structure funciton for GSNE
.makegsneESF <- function(g0, s0, n0, e0) {
  esf.par <- .gsne.la(g0,s0,n0,e0)
  
  esf.par.info <- esf.par[-(1:2)]
  thisZ <- esf.par$Z
  esf.par <- esf.par$La
  
  return(list(La=esf.par, La.info=esf.par.info, Z=thisZ, 
              state.var=c(G0=g0,S0=s0,N0=n0,E0=e0)))
}


## function to find Lagrange multipliers numerically (uses optimization on functions wrapping C++ code)
.gsne.la <- function(G0, S0, N0, E0, ftol=1e-08) {
  this.mn.grid <- expand.grid(1:S0,1:N0)
  
  # best guesses for la1 and la2
  init.la1 <- .xlogx(G0/S0)
  init.beta <- .xlogx(S0/N0)
  init.la2 <- init.beta - G0/(E0-N0)
  
  sol <- nleqslv::nleqslv(c(0.1, 0.01), fn=.gsne.la.eq,
                 G0=G0, S0=S0, N0=N0, E0=E0,
                 method="Newton", xscalm="auto",
                 control=list(ftol=ftol))
  
  sol$x <- c(sol$x, G0/(E0-N0))
  thisZ <- .gsneZ(sol$x, S0, N0)
  
  return(list(La=sol$x, Z=thisZ,
              term.code=sol$termcd, term.message=sol$message,
              fn.val=sol$fvec))
}

## function to compute GSNE constraint equations
.gsne.la.eq <- function(La, G0, S0, N0, E0) {
  la1 <- La[1]
  la2 <- La[2]
  la3 <- G0/(E0-N0)
  
  beta <- la2 + la3
  
  .Call('meteR_gsne_la', PACKAGE = 'meteR', G0, S0, N0, la1, beta)
}

## function to compute normalizing constant
.gsneZ <- function(La,S0,N0) {
  la1 <- La[1]
  la2 <- La[2]
  la3 <- La[3]
  
  .Call('meteR_gsne_Z', PACKAGE = 'meteR', la1, la2, la3, S0, N0)
}

## function to solve general case of y = -x*log(x) when y < 1
## used to find starting values for la1 and beta in METE
.xlogx <- function(y) {
  if(y > 0.367) {
    return(0.4)
  } else {
    out <- uniroot(function(x) {-y -x*log(x)},interval=c(1e-16,0.367))
    return(out$root)
  }
}


## clean up C++ stuff on unloading the package
.onUnload <- function (libpath) {
  library.dynam.unload("mypackage", libpath)
}
