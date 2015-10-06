## figure out areas
.findAreas <- function(spp, abund, row, col, x, y, Amin, A0) {
    if(missing(spp) | missing(abund)) { # no data
        if(missing(Amin)) { # use row and col as nrow and ncol
            if(length(row) == 1 & length(col) == 1) {
                ## for simplicity make sure nrow <= ncol
                nrow <- min(row, col)
                ncol <- max(row, col)
                Amin <- A0/(nrow*ncol)
            } else {
                stop('row and col should be scalers and interpreted as number of rows and cols when A0 provided')
            }
        } else { # use Amin to figure out nrow and ncol
            nrow <- ncol <- floor(sqrt(A0/Amin))
        }
        row <- col <- NULL
    } else { # data provided
        if(!missing(x) & !missing(y)) { # case where x,y data provided turn point data into row col data
            if(!missing(row) & !missing(col)) { # if row and col given use them to get nrow and ncol
                if(length(row) == 1 & length(col) == 1) {
                    ## for simplicity make sure nrow <= ncol
                    nrow <- min(row, col)
                    ncol <- max(row, col)
                } else {
                    ## if row and col are vectors doesn't make sense with point data
                    stop('if using x,y location data, row and col should be scalers indicating number of desired rows and columns')
                }
            } else { # if row and col not given use max extent of x,y data and Amin to make grid
                xrng <- diff(range(x))
                yrng <- diff(range(y))
                ## make sure row (=y) <= col (=x)
                if(xrng < yrng) {
                    temp <- y
                    y <- x
                    x <- temp
                    rng <- diff(range(x))
                    yrng <- diff(range(y))
                }

                ## try to get nrow and ncol such that min cell is as close to square as possible
                A0 <- xrng*yrng
                ncell <- floor(A0/Amin)
                nrow <- floor(yrng/sqrt(Amin))
                ncol <- round(ncell/nrow)
                ncell <- nrow * ncol
                Amin <- A0/ncell
            }
            
            ## now we have nrow and ncol, use those to make grid ##############  NEED TO DO  ################
            
        } else if(length(row) != length(spp) | length(col) != length(spp)) { # case where `row' and `col' provided but not for all records
            stop('either row and column must be given for each spp or individual, or x,y coordinates given')
        } else { # case where row and col provided for each record
            ## make sure row and col IDs are 1:nrow and 1:ncol
            if(!is.numeric(row) | max(row) != length(unique(row))) row <- as.numeric(as.factor(row))
            if(!is.numeric(col) | max(col) != length(unique(col))) col <- as.numeric(as.factor(col))
            
            ## number of cells in each direction, transpose data if needed so always more columns
            if(max(row) > max(col)) {
                temp <- col
                col <- row
                row <- temp
            }
            
            ## get nrow and ncol
            nrow <- max(row)
            ncol <- max(col)
        }
    }
    
    ## now regardless of input we have nrow and ncol and if there are real
    ## data we have row and col identities for each entry in the data.
    
    ## figure out vector of sizes in units of cells; right now only doublings supported
    maxDoubling <- .calcMaxDoubling(floor(log(nrow*ncol) / log(2)), nrow)
    
    return(list(areas = 2^(0:maxDoubling), row=row, col=col, nrow=nrow, ncol=ncol, Amin=Amin, A0=A0))
}

## funciton uses recursion to calculate maximum possible doubling for an area with 
## given minimum dimension
.calcMaxDoubling <- function(doub, mindim) {
	if(doub %% mindim != 0) doud <- doub - 1
	
	if(2^doub > 2*mindim^2) {
		doub <- .calcMaxDoubling(floor(log(2*mindim^2) / log(2)), mindim)
	}
	
	return(doub)
}

## function to find neighbors and return groups of neighbors of given size `a'
.getNeighbors <- function(a, nr, nc) {
    foo <- function() {
        addToCol <- 1:(nc - c2 + 1)
        addToRow <- 1:(nr - c1 + 1)
        
        groups <- vector('list', length(addToRow) * length(addToCol))
        for(i in 1:length(addToRow)) {
            for(j in 1:length(addToCol)) {
                temp <- expand.grid(col=addToCol[j] + 0:(c2-1),
                                    row=addToRow[i] + 0:(c1 - 1))
                groups[[i + length(addToRow)*(j-1)]] <- paste(temp[, 2], temp[, 1], sep=',')
            }
        }
        
        return(data.frame(group=rep(1:(length(addToRow) * length(addToCol)), 
                              each=a),
                          cells=unlist(groups), stringsAsFactors=FALSE))
    }
    
    c1 <- max((1:nr)[a %% (1:nr) == 0])
    c2 <- a/c1
    groups1 <- foo()
    
    if(c2 < nr & c1 != c2) {
        c1 <- c2
        c2 <- a/c1
        groups2 <- foo()
        groups2[, 1] <- groups2[, 1] + max(groups1[, 1])
        return(rbind(groups1, groups2))
    } else {
        return(groups1)
    }
}

## function to get number of species in groups of cells
.getSppInGroups <- function(spp, abund, row, col, groups, endemics=FALSE) {
    cellID <- paste(row, col, sep=',')
    cellGroup <- groups$group[match(cellID, groups$cells)]

    sppByGroup <- tapply(abund, list(cellGroup, spp), sum)
    sppByGroup[is.na(sppByGroup)] <- 0
    sppByGroup[sppByGroup > 0] <- 1
    
    if(endemics) {
    	theseEndem <- colSums(sppByGroup) == 1
    	return(rowSums(sppByGroup[, theseEndem, drop=FALSE]))
    } else {
    	return(rowSums(sppByGroup))
    }
}

## makeSSF is not vectorized and too slow so make special
## function to extract \lambda_\Pi and make it vectorized
## over n0
.getPi0 <- function(n0,A,A0) {
	if(A/A0 == 0.5) {
		pi0 <- 1/(1+n0)
	} else {
		eq52 <- .useEq52(n0,A,A0)
		pi0 <- numeric(length(eq52))
		
		if(any(eq52)) {
			# cat(sprintf('using eq52 approx %s times \n', sum(eq52)))
			pi0[eq52] <- 1 - n0[eq52]/(n0[eq52] + A0/A)
		}
		
		if(any(!eq52)) {
			# cat(sprintf('using exact sol %s times \n', sum(!eq52)))
			pi0[!eq52] <- sapply(n0[!eq52], function(n) {
				.mete.Pi(0, metePi(0, n0=n, A, A0)$La, n0=n)
			})
		}
	}
	
	return(pi0)
}

## two helper functions giving upscaling constraints
.upCon1 <- function(beta, Sup, N0, S0) {
	## eq 8 in Harte et al. 2009 Ecol Lett
	2*Sup*exp(beta) - 2*N0 * 
	    (1 - exp(-beta))/(exp(-beta) - exp(-beta*(2*N0+1))) * 
	    (1 - exp(-beta*2*N0) / (2*N0+1)) -
	    S0 # when solved should return 0
}

.upCon2 <- function(beta, N0, Sup) {
	## simplification of eq 9 in Harte et al. 2009 Ecol Lett
	Sup/(2*N0) *
	    ((1-exp(-beta)^(2*N0-1))/(1-exp(-beta)) - 1) -
	    log(1/beta)
}
