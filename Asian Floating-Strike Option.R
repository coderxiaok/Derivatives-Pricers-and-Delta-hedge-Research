## Asian Floating-Strike Call Option
## Author: Ming Tian
## Date: 2015-05-06
## Copyright @2015 MingTian

#assign the parameters firstly
mu = 0  #mean
sig = 0.2  #sigma/annulaized volitility
S0 = 100  #start stock price
r = 0  #interest rate
q = 0  #dividends of stock
dt = 1/252  #time step
T_m = 1  #maturity time
type = "c"  #option type

# 20150913 update version
Random.Path.New <- function(S0_ = S0, sig_ = sig, mu_ = mu, dt_ = dt, T_m_ = T_m, n_path) {
  # generate the random vector
  rand_vec <- rlnorm(n=(T_m_/dt_)*n_path, meanlog = mu_ - 0.5*((sig_*sqrt(1/252))^2), sdlog = sig_*sqrt(1/252))
  
  # build the result
  result <- matrix(data=rand_vec, nrow=n_path, ncol=(T_m_/dt_))
  
  # apply the function to transform to the price ts
  result[, 1] <- S0_
  result_new <- t(apply(result, 1, function(x) {
    cumprod(x)
  }))
  
  # return the result
  result_new
}

#define a payoff function for Asian Floating Strike Option
AFS.Payoff <- function(Price_ts, type = "c")
{
  if(type == "c")  #Call
  {
    payoff_tmp <- Price_ts[length(Price_ts)] - mean(Price_ts)
    payoff <- max(0, payoff_tmp)
  }
  else  #Put
  {
    payoff_tmp <- mean(Price_ts) - Price_ts[length(Price_ts)]
    payoff <- max(0, payoff_tmp)
  }
  
  # return the result
  payoff
}

# 20150917 updated
AFS.Payoff.Vectorization <- function(PayoffFunc, S0_ = S0, sig_ = sig, sig_DH = sig, 
                                     mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000, type = "c")
{
  #Call the Random.Path function firstly
  rand_path <- Random.Path.New(n_path = n_path_)
  
  #build a result container
  result <- matrix(1:(n_path_*2), ncol = 2)
  result[, 1] <- rand_path[, ncol(rand_path)]
  #decide the option type
  if(type == "c")  #Call
  {
    payofffunc <- function(x) {
      # calculate option Payoff
      Opt_payoff <- PayoffFunc(Price_ts = x, type = type)
      # return
      Opt_payoff
    }
    
    result[, 2] <- apply(rand_path, 1, payofffunc)
  }
  else  #Put
  {
    payofffunc <- function(x) {
      # calculate option Payoff
      Opt_payoff <- PayoffFunc(Price_ts = x, type = type)
      # return
      Opt_payoff
    }
    
    result[, 2] <- apply(rand_path, 1, payofffunc)
  }
  #return
  result
}

# call test
test <- AFS.Payoff.Vectorization(PayoffFunc = AFS.Payoff)
plot(test[, 1], test[, 2], col="blue", cex=0.5)

# put test
test <- AFS.Payoff.Vectorization(PayoffFunc = AFS.Payoff, type = "p")
plot(test[, 1], test[, 2], col="blue", cex=0.5)
