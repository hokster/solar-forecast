# file to generate synthetic data and try to fit model
library(MASS)

# load MCMC files
source("sample_states_fixed.R")
source("sample_dirichlet.R")
source("sample_regression.R")
source("permute.R")

set.seed(1)

# generate synthetic data

# set parameters
T <- 100
trans_prob.actual <- matrix(c(0.2,0.5,0.3,
                              0.05,0.6,0.35,
                              0.5,0.4,0.1),
                            nrow=3,byrow=TRUE)
beta.actual <- array(NA,dim=c(3,3,2))
beta.actual[,,1] <- matrix(c(0.75,0.25,0.5,
                        1,-0.5,-0.05,
                        0.25,0.5,0),
                      nrow=3,byrow=TRUE)
beta.actual[,,2] <- matrix(c(0.25,0.5,0.25,
                             0.75,-0.5,-0.05,
                             0,0.45,0),
                           nrow=3,byrow=TRUE)

cov.actual <- list()
cov.actual[[1]] <- matrix(c(0.8,0.75,
                            0.75,0.8),
                          nrow=2,byrow=TRUE)
cov.actual[[2]] <- matrix(c(0.3,0.25,
                            0.25,0.3),
                          nrow=2,byrow=TRUE)
cov.actual[[3]] <- matrix(c(0.5,0.2,
                            0.2,0.5),
                          nrow=2,byrow=TRUE)

# sample states and generate data
states.actual <- mat.or.vec(T,1)
dat <- mat.or.vec(T,2)
dat[1,] <- c(1,0.8)
dat[2,] <- c(0.8,0.5)
state0 <- sample(1:3,size=1)
for (t in 1:T) {
  if (t == 1) {
    states.actual[t] <- sample(1:3,size=1,prob=trans_prob.actual[state0,])
  }
  else if (t == 2) {
    states.actual[t] <- sample(1:3,size=1,prob=trans_prob.actual[states.actual[t-1],])
  }
  else {
    states.actual[t] <- sample(1:3,size=1,prob=trans_prob.actual[states.actual[t-1],])
    dat[t,1] <- sum(c(1,dat[c(t-1,t-2),1])*beta.actual[states.actual[t],,1])
    dat[t,2] <- sum(c(1,dat[c(t-1,t-2),2])*beta.actual[states.actual[t],,2])
    dat[t,] <- dat[t,] + mvrnorm(n=1,mu=mat.or.vec(1,2),Sigma=cov.actual[[states.actual[t]]])
  }
}

save.image('synthdata.RData')