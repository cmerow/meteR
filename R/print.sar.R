print.sar <- function(x) {
	cat(sprintf('%s %s area relationship ranging from \n', attr(x, 'source'), 
	            ifelse(attr(x, 'type')=='sar', 'species', 'endemics')))
	
	cat(sprintf('A: [%s, %s] \nS: [%s, %s] \n', 
	            min(x[, 'A'], na.rm=TRUE), 
	            max(x[, 'A'], na.rm=TRUE), 
	            min(x[, 'S'], na.rm=TRUE), 
	            max(x[, 'S'], na.rm=TRUE)))
	
	invisible(x)	
}