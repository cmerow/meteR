##	print and summary methods for METE object

print.METE <- function(x) {
	cat("METE object with state variables:\n")
	print(x[[1]]$state.var)
	
	cat("\n");cat("with Lagrange multipliers:\n")
	print(x[[1]]$La)
	
	cat("\n");cat("predicting the distributions:\n")
	cat(names(x)[-1]);cat("\n")
}

# summary.METE <- function(x) {
	
# }