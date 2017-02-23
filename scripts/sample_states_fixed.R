require(mvtnorm)

# function to forward-filter the states
filter_forward <- function(data,X,beta,sigma,num_states,trans_prob) {
  T <- dim(data)[1]  
  # allocate memory for probabilities
  filter_prob <- mat.or.vec(T,num_states)
  
  # loop forward in time to compute filtered probabilities
  for (t in 1:T) {
    # predict initial state probabilities, assuming uniform distribution of S0
    if (t == 1) {
      pred_prob <- colMeans(trans_prob)
    }
    else {
      pred_prob <- filter_prob[t-1,] %*% trans_prob
    }
    
    # filter state probabilities
    means <- rowSums(sweep(beta,1,X[t,],FUN='*')) # calculate mean of data normal distribution
    lik <- sapply(1:num_states, function(k) {dmvnormd(data[t,],means[,k],sigma[k,,])})
    filter_prob[t,] <- lik*pred_prob
    filter_prob[t,] <- filter_prob[t,]/sum(filter_prob[t,])
  }
  
  return(filter_prob)
}

# function to multi-move sample the state indicators
sample_backward <- function(T,num_states,filter_prob,trans_prob) {
  states <- mat.or.vec(T,1)
  # sample the terminal state S_T
  states[T] <- sample(1:num_states,1,prob=filter_prob[T,])
  
  # backwards sample from the relevant conditional distribution
  for (t in T-1:1) {
    # calculate conditional probabilities
    cond_prob <- sweep(trans_prob[,states[t+1]],1,filter_prob[t,],FUN='*')
    cond_prob <- sweep(cond_prob,1,rowSums(cond_prob),FUN='*')
    # sample state S_t
    states[t] <- sample(1:num_states,1,prob=cond_prob)
  }
  
  return(states)
}