# use this file for time-varying transition probabilities

require(MASS)

# randomly sample differenced utilities, per Fruhwith & Fruhwith-Schnattner (2010)
# Inputs:
#   T: number of points in the time series
#   num_states: number of states
#   Z: matrix of transition covariates. Each row corresponds to the covariates for a time t.
#   kappa: matrix of transition regression cofficients. Each column corresponds to the coefficients
#         for a state k.
#   states: vector with sampled state sequence.
sample_utilities <- function(T,num_states,Z,kappa,states) {
  # create state matrix
  states <- as.factor(states)
  D <- model.matrix(~states-1)
  # sample uniform values
  W <- matrix(runif(T*num_states),nrow=T,ncol=num_states)
  # calculate regression probabilities
  reg <- Z%*%kappa
  # calculate log(lambda)
  log.lambda <- sapply(1:num_states,function(i) {log(rowSums(exp(reg)[,-i]))})
  # calculate the probability of being in each state
  state.prob <- 1-(exp(log_lambda-reg)/(1+exp(log_lambda-reg)))
  # sample differenced utilities
  logit.quant <- D+W*(1-D-state.prob)
  utils <- reg-log.lambda+log(logit.quant)-log(1-logit.quant)
  return (list(utils,reg,log.lambda))
}

sample_component <- function(T,num_states,utils,reg,log.lambda) {
  w <- c(1.8446,17.268,37.393,31.697,10.89,0.90745)/100
  s.sq <- c(0.68159,1.2419,2.2388,4.0724,7.4371,13.772)
  # calculate probabilities of being in each component
  mult.prob <- sapply(1:6,function(r) {(w[r]/sqrt(s.sq[r])*exp(0.5*(utils+log.lambda-reg)^2/s.sq[r]))})
  mult.prob <- sweep(mult.prob,3,apply(mult.prob,3,sum),FUN='/')
  components <- outer(1:T,1:num_states-1,FUN=Vectorize(function(t,k) {
    sample(1:6,1,prob=mult.prob[t,k,])
  }))
  return(components)
}

sample_kappa <- function(k,utils,Z,comps) {
  s.sq <- c(0.68159,1.2419,2.2388,4.0724,7.4371,13.772)
  cov <- solve(diag(s.sq[comps[,k]]))
  Sk <- utils[,k]
  kappa.var <- solve(t(Z)%*%cov%*%Z)
  kappa.mean <- Sigma%*%t(Z)%*%cov%*%Sk
  kappa <- mvrnorm(n=1,mu=kappa.mean,Sigma=kappa.var)
  return(kappa)
}