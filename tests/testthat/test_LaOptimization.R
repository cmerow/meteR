# mete.lambda <- function(S0, N0, E0) {
#   ## reasonable starting values
#   init.la2 <- S0/(E0-N0)
#   beta.guess <- 0.001
#   
#   ## uses eq 7.29 from Harte 2011, catches possible error and assigns a default guess
#   init.beta <- try(uniroot(function(b) {
#     (1 - exp(-b)) / (exp(-b) - exp(-b*(N0+1))) * log(1/b) - S0/N0
#   }, interval=c(1/N0, S0/N0)), silent=TRUE)
#   
#   if(class(init.beta) == 'try-error') {
#     init.beta <- beta.guess
#   } else {
#     init.beta <- init.beta$root
#   }
#   
#   init.la1 <- init.beta - init.la2
#   
#   ## the solution
#   la.sol <- nleqslv::nleqslv(x=c(la1 = init.la1, la2 = init.la2), 
#                              fn = .la.syst2, S0=S0,N0=N0,E0=E0,
#                              control=list(ftol=.Machine$double.eps))
#   
#   return(list(lambda=c(la1=la.sol$x[1], la2=la.sol$x[2]), 
#               syst.vals=la.sol$fvec, converg=la.sol$termcd,
#               mesage=la.sol$message, nFn.calc=la.sol$nfcnt, nJac.calc=la.sol$njcnt))
# }
# 
# 
# meteTable <- read.csv('~/Dropbox/mete/mete_table7.2.csv')
# 
# comp <- t(mapply(function(S0, N0, E0) {
#   out <- mete.lambda(S0, N0, E0)
#   c(out$lambda, out$converg)
# }, S0=s0, N0=n0, E0=e0))
# 
# comp <- as.data.frame(comp)
# colnames(comp) <- c('la1', 'la2', 'conv')
# plot(comp$conv, comp$la1 - meteTable$la1);abline(h=0)
