# # @title Class "meteESF"
# #' 
# #' @description
# #' Class to hold the Ecosystem Structure Function $R(n, epsilon)$
# #'
# #' @return Contains the elements
# #' \describe{
# #'   \item{data}{The data used to construct the ESF}
# #'   \item{emin}{The minimum metabolic rate used to rescale metabolic rates}
# #'   \item{La}{Vector of Lagrange multipliers}
# #'   \item{La.info}{Termination information from optimization procedure}
# #'   \item{state.var}{State variables used to constrain entropy maximization}
# #'   \item{Z}{Normalization constant for ESF}
# #' }
# #'
# #' @method Methods defined for meteESF
# #' \describe{
# #'   \item{print}{print summary of ESF}
# #'   \item{sad}{extract species abundance distribution from ESF}
# #'   \item{ipd}{extract species abundance distribution from ESF}
# #'   \item{sipd}{extract species abundance distribution from ESF}
# #' }
# #'
# #' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# #' @seealso meteDist, meteRelat, meteSSF
# #' @references Harte, J. 2011. Maximum entropy and ecology: 
# #' a theory of abundance, distribution, and energetics. Oxford University Press.
# "meteESF-class"
# 
# #==========================================================================
# 
# 
# #' @title Class "meteSSF"
# #' 
# #' @description
# #' Class to hold the Spatial Structure Function $Pi(n)$
# #'
# #' @return Contains the elements
# #' \describe{
# #'   \item{data}{The data used to construct the SSF}
# #'   \item{emin}{The minimum metabolic rate used to rescale metabolic rates}
# #'   \item{La}{Vector of Lagrange multipliers}
# #'   \item{La.info}{Termination information from optimization procedure}
# #'   \item{state.var}{State variables used to constrain entropy maximization}
# #'   \item{Z}{Normalization constant for SSF}
# #' }
# #'
# #' @method Methods defined for meteSSF
# #' \describe{
# #'   \item{print}{print summary of SSF}
# #'   \item{ssad}{extract species abundance distribution from SSF}
# #' }
# #'
# #' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# #' @seealso meteDist, meteRelat, meteESF
# #' @references Harte, J. 2011. Maximum entropy and ecology: 
# #' a theory of abundance, distribution, and energetics. Oxford University Press.
# "meteSSF-class"
# 
# #==========================================================================
# 
# 
# #' @title Class "meteDist"
# #' 
# #' @description
# #' Class to hold probability distributions predicted by METE
# #'
# #' @return Contains the elements
# #' \describe{
# #'   \item{data}{The data used to construct the prediction}
# #'   \item{d}{density funciton}
# #'   \item{p}{cumulative density function}
# #'   \item{q}{quantile funtion}
# #'   \item{r}{random number generator}
# #'   \item{La}{Vector of Lagrange multipliers}
# #'   \item{state.var}{State variables used to constrain entropy maximization}
# #'   \item{type}{Type of distribution returned}
# #' }
# #'
# #' @method Methods defined for meteDist
# #' \describe{
# #'   \item{print}{print summary of distribution}
# #'   \item{plot}{plot theoretical prediction and data}
# #'   \item{logLik}{calculate log likelihood}
# #'   \item{residuals}{calculate residuals between data and prediction, either the CDF or rank distribution}
# #' }
# #'
# #' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# #' @seealso meteESF, meteRelat, meteSSF
# #' @references Harte, J. 2011. Maximum entropy and ecology: 
# #' a theory of abundance, distribution, and energetics. Oxford University Press.
# "meteDist-class"
# 
# 
# #==========================================================================
# 
# 
# #' @title Class "meteRelat"
# #' 
# #' @description
# #' Class to hold relationships predicted by METE; relationships are defined as predictions 
# #' that relate one variable to another, e.g. area to species
# #'
# #' @return Contains the elements
# #' \describe{
# #'   \item{obs}{The data used to construct the prediction}
# #'   \item{pred}{The minimum metabolic rate used to rescale metabolic rates}
# #' }
# #'
# #' @method Methods defined for meteRelat
# #' \describe{
# #'   \item{print}{print summary of distribution}
# #'   \item{plot}{plot theoretical prediction and data}
# #'   \item{residuals}{calculate residuals between observed and predicted values}
# #' }
# #'
# #' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
# #' @seealso meteESF, meteDist, meteSSF
# #' @references Harte, J. 2011. Maximum entropy and ecology: 
# #' a theory of abundance, distribution, and energetics. Oxford University Press.
# "meteRelat-class"
