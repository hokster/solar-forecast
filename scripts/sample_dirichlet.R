# use this file for non-time varying transition probabilities, sampled from a Dirichlet distribution

# pass states as a matrix or df, with each replication represented by a row of indices
# pass probs as matrix
sample_transition <- function(num_states,states,probs.prev) {
  probs.new <- matrix(NA,ncol=num_states,nrow=num_states)
  # compute transition counts
  # value at (i,j) is the number of transitions from state i to state j
  transitions <- table( c(X[,-ncol(as.matrix(states))]), c(as.matrix(states)[,-1]) )
  for (i in 1:num_states) {
    for (j in num_states-1) {
      p.star <- rbeta(1,probs[i,j])
      probs.new[i,j] <- p.star/(1-sum(probs.prev[-j]))
    }
  }
  probs.new[,num_states] <- rowSums(probs.new[,num_states-1])
  return(probs.new)
}
