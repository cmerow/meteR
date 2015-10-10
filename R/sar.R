#' @title Compute METE species area relationship (SAR)
#'
# @description
#'
# @details
#' 
#' 
#' @param spp
#' @param abun 
#' @param row 
#' @param col
#' @param x 
#' @param y 
#' @param S0
#' @param N0
#' @param Amin
#' @param A0
#' @param upscale
#' @param EAR
#' 
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))

# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

meteSAR <- function(spp, abund, row, col, x, y, S0 = NULL, N0 = NULL,
                    Amin, A0, upscale=FALSE, EAR=FALSE) {    
  ## figure out vector of sizes in units of cells; right now only doublings supported
  ## not needed if upscale is TRUE
  if(!upscale) {
    areaInfo <- .findAreas(
      spp=if(missing(spp)) NULL else spp,
      abund=if(missing(abund)) NULL else abund, 
      row=if(missing(row)) NULL else row, 
      col=if(missing(col)) NULL else col, 
      x=if(missing(x)) NULL else x, 
      y=if(missing(y)) NULL else y, 
      Amin=if(missing(Amin)) NULL else Amin, 
      A0=if(missing(A0)) NULL else A0)
    areas <- areaInfo$areas
    row <- areaInfo$row
    col <- areaInfo$col
    nrow <- areaInfo$nrow
    ncol <- areaInfo$ncol
    Amin <- areaInfo$Amin
    A0 <- areaInfo$A0
  }
  
  if(upscale & EAR) stop('upscaling EAR not currently supported')
  
  ## the ESF
  if(!missing(spp) & !missing(abund)) {
    S0 <- length(unique(spp))
    N0 <- sum(abund)
  }
  if(is.null(S0) | is.null(N0)) stop('must provide spp and abund data or state variables S0 and N0')
  thisESF <- meteESF(S0=S0, N0=N0)
  
  ## calculate empirical SAR
  if(!missing(spp) & !missing(abund)) {
    eSAR <- empiricalSAR(spp, abund, row=row, col=col, Amin=Amin, A0=A0, EAR=EAR)
  } else {
    eSAR <- NULL
  }
  
  ## calculate theoretical SAR
  if(upscale) {
    thrSAR <- upscaleSAR(thisESF, Amin, A0, EAR)
  } else {
    thrSAR <- downscaleSAR(thisESF, areas*Amin, A0, EAR)
  }
  
  out <- list(obs=eSAR, pred=thrSAR)
  class(out) <- 'meteRelat'
  
  return(out)
}



#================================================================
#' @title Empirical SAR
#'
# @description
#'
# @details
#' 
#' 
#' @param spp
#' @param abun 
#' @param row 
#' @param col
#' @param x 
#' @param y 
#' @param Amin
#' @param A0
#' @param EAR logical. TRUE computes the endemics-area relatinship
#' @param upscale logical. Do upscaling?
#' 
#' @export
#' 
#' @examples
#' esf=meteESF(spp=anbo$spp,
#'              abund=anbo$count)

# @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.


empiricalSAR <- function(spp, abund, row, col, x, y, Amin, A0, EAR=FALSE) {
  ## figure out vector of sizes in units of cells; right now only doublings supported
  areaInfo <- .findAreas(
    spp=if(missing(spp)) NULL else spp,
    abund=if(missing(abund)) NULL else abund, 
    row=if(missing(row)) NULL else row, 
    col=if(missing(col)) NULL else col, 
    x=if(missing(x)) NULL else x, 
    y=if(missing(y)) NULL else y, 
    Amin=if(missing(Amin)) NULL else Amin, 
    A0=if(missing(A0)) NULL else A0)
  areas <- areaInfo$areas
  row <- areaInfo$row
  col <- areaInfo$col
  nrow <- areaInfo$nrow
  ncol <- areaInfo$ncol
  Amin <- areaInfo$Amin
  A0 <- areaInfo$A0
  
  ## loop over areas
  out <- lapply(areas, function(a) {
    nspp <- .getSppInGroups(spp, abund, row, col, .getNeighbors(a, nrow, ncol), EAR)
    cbind(A=a*Amin, S=nspp)
  })
  out <- do.call(rbind, out)
  
  ## make output of class `sar' and tell it about empirical v. theoretical and ear v. sar
  attr(out, 'source') <- 'empirical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}



#================================================================
#' @title Downscale the species area relationship (SAR)
#'
# @description
#'
# @details
#' 
#' 
#' @param x an object of class meteESF
#' @param A
#' @param A0
#' @param EAR logical. TRUE computes the endemics-area relatinship
#'  
#' @export
#' 
#' @examples
#' esf=meteESF(spp=anbo$spp,
#'              abund=anbo$count)
#'              
#' @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

downscaleSAR <- function(x, A, A0, EAR=FALSE) {
  n0 <- 1:x$state.var['N0']
  
  ## difference between EAR and SAR is for EAR we get Pi(n0) [fun .getPin0]
  ## and for SAR we get 1 - Pi(0) [1 - .getPi0]
  if(EAR) {
    piFun <- function(a) .getPin0(n0, a, A0)
  } else {
    piFun <- function(a) 1 - .getPi0(n0, a, A0)
  }
  
  ## function to get species number at scale `a'
  getspp <- function(a) {
    probs <- piFun(a) * 
      with(x, 
           metePhi(n0, La[1], La[2], Z, 
                   state.var['S0'], state.var['N0'], 
                   ifelse(is.na(state.var['E0']), 1e+06, state.var['E0'])))
    
    return(x$state.var['S0'] * sum(probs))
  }
  
  ## loop over A
  nspp <- sapply(A, getspp)
  
  ## should return matrix with column for area and column for spp
  out <- cbind(A=A, S=nspp)
  attr(out, 'source') <- 'theoretical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}




#================================================================
#' @title upscale SAR
#'
# @description
#'
# @details
#' 
#' 
#' @param x an object of class meteESF
#' @param A0
#' @param Aup
#' @param EAR logical. TRUE computes the endemics-area relatinship
#' 
#' @export
#' 
#' @examples
#' esf=meteESF(spp=anbo$spp,
#'              abund=anbo$count)
#'              
#' @return list
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# @seealso sad.mete, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

upscaleSAR <- function(x, A0, Aup, EAR=FALSE) {
  ## vector of areas starting with anchor area A0
  Aups <- A0 * 2^(0:ceiling(log(Aup/A0)/log(2)))
  
  ## vector of abundances at each area
  N0s <- x$state.var['N0'] * 2^(0:ceiling(log(Aup/A0)/log(2)))
  
  ## vector of number of species at each area
  S0s <- numeric(length(Aups))
  S0s[1] <- x$state.var['S0']
  
  ## vector to hold termination codes from nleqslv about whether optimization succeeded
  termcodes <- numeric(length(Aups))
  
  ## need to recursively solve constraint fun (solution in `.solveUpscale') up to Aup
  for(i in 2:length(Aups)) {
    S0s[i] <- .solveUpscale(S0s[i-1], N0s[i-1])
  }
  
  ## should return matrix with column for area and column for spp
  out <- cbind(A=Aups, S=S0s)
  attr(out, 'source') <- 'theoretical'
  attr(out, 'type') <- ifelse(EAR, 'ear', 'sar')
  class(out) <- 'sar'
  
  return(out)
}
