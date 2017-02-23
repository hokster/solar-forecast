
# randomly permute the states
# beta and kappa is a matrix with the regression coefficients for each state 
# passed as a row.
# cov is a list of matrices
random_permute <- function(num_states,kappa,beta,cov) {
  state.permute <- sample(1:num_states)
  kappa.permute <- kappa[state.permute,]
  beta.permute <- beta[state.permute,,]
  cov.permute <- cov[state.permute]
  return(list(kappa.permute,beta.permute,cov.permute))
}