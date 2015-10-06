print.sar <- function(x) {
	cat(sprintf('%s %s area relationship ranging from \n', attr(x, 'source'), 
	            ifelse(attr(x, 'type')=='sar', 'species', 'endemics')))
	
	cat(sprintf('A: [%s, %s] \nS: [%s, %s] \n', min(x[, 'A']), max(x[, 'A']), min(x[, 'S']), max(x[, 'S'])))
	
	invisible(x)	
}