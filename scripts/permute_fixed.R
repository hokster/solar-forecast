# randomly permute the states
# beta is a list of matrices with the regression coefficients for each state 
# passed as a row. beta[[1]] is the first output variable.
# cov is a list of matrices
# trans_prob is a matrix of transition probabilities
random_permute <- function(num_states,trans_prob,beta,cov) {
  state.permute <- sample(1:num_states)
  trans_prob.permute <- trans_prob[state.permute,state.permute]
  beta.permute <- beta[[1]][state.permute,,]
  cov.permute <- cov[state.permute]
  return(list(trans_prob.permute,beta.permute,cov.permute))
}