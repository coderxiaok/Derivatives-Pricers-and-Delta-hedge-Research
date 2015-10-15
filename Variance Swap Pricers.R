## Variance Swap Pricers
## Author: Ming Tian
## Copyright @2015 Ming Tian

# Black-Scholes Pricer
BSPricer <- function(S, K, r, Ttm, vol, type) {
  d1 <- (1 / (vol * sqrt(Ttm))) * (log(S / K) + (r + vol^2/2) * Ttm)
  d2 <- d1 - vol * sqrt(Ttm)
  if (type == "c") {
    # call price
    CallPrice <- S * pnorm(d1) - K * exp(-r * Ttm) * pnorm(d2)
    # return
    return(CallPrice)
  }
  else if (type == "p") {
    # put price
    PutPrice <- pnorm(-d2) * K * exp(-r * Ttm) - pnorm(-d1) * S
    # return
    return(PutPrice)
  }
  else {
    stop("Wrong input type!")
  }
}

# Black-Scholes Gamma
BSGamma <- function(S, K, r, Ttm, vol) {
  d1 <- (1 / (vol * sqrt(Ttm))) * (log(S / K) + (r + vol^2/2) * Ttm)
  return(exp(0)*dnorm(d1)/(S*vol*sqrt(Ttm)))
}

# JPMorgan Model Pricer
VarSwapPricer.SimpleVersion <- function(S, put_strikes, call_strikes, put_IVs, call_IVs, Ttm=1, r, Kvs) {
  # check the inputs 
  if (length(put_strikes) != length(put_IVs)) {
    stop("Put strikes are not as long as put IVs, please check!")
  }
  if (length(call_strikes) != length(call_IVs)) {
    stop("Call strikes are not as long as Call IVs, please check!")
  }
  
  # compute forward price firstly
  fprice <- S * exp(r*Ttm)
  
  # calculate the price of puts
  put_strikes_abs <- fprice * put_strikes
  put_series <- apply(cbind(put_strikes_abs, put_IVs), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "p")
  })
  put_series <- put_series / fprice
  
  # calculate the price of calls
  call_strikes_abs <- fprice * call_strikes
  call_series <- apply(cbind(call_strikes_abs, call_IVs), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "c")
  })
  call_series <- call_series / fprice
  
  # sum the puts and calls
  len1 <- length(put_strikes)
  put_k0 <- put_strikes[1] - (put_strikes[2] - put_strikes[1])
  put_weight <- (put_strikes[1:len1] - c(put_k0, put_strikes[1:(len1-1)])) / (put_strikes^2)
  put_weight[put_strikes == 1] <- put_weight[put_strikes == 1] / 2
  sum1 <- sum(put_series * put_weight)
  len2 <- length(call_strikes)
  call_k0 <- call_strikes[1] - (call_strikes[2] - call_strikes[1])
  call_weight <- (call_strikes[1:len2] - c(call_k0, call_strikes[1:(len2-1)])) / (call_strikes^2)
  call_weight[call_strikes == 1] <- call_weight[call_strikes == 1] / 2
  sum2 <- sum(call_series * call_weight)
  
  # var swap price
  varswap <- (2/Ttm) * (sum1 + sum2) - exp(-r*Ttm) * (Kvs^2)
  
  # return 
  return(varswap)
}

# # JPMorgan Model Pricer test
# put_strikes_rate <- seq(from = 0.5, to = 1, by = 0.05)
# call_strikes_rate <- seq(from = 1, to = 1.5, by = 0.05)
# put_IVs_tmp <- c(0.276, 0.264, 0.252, 0.240, 0.227, 0.214, 0.200, 0.187, 0.173, 0.160, 0.148)
# call_IVs_tmp <- c(0.148, 0.137, 0.129, 0.122, 0.119, 0.118, 0.119, 0.121, 0.125, 0.129, 0.134)
# 
# # two input parameters
# discount_factor <- 0.977368853
# forward_price <- 2935.02
# Kvs_tmp <- 0.16625
# 
# rf=log(1/discount_factor)
# S0 <- forward_price / exp(rf)
# VarSwapPricer.SimpleVersion(S = S0, put_strikes = put_strikes_rate, call_strikes = call_strikes_rate, 
#                             put_IVs = put_IVs_tmp, call_IVs = call_IVs_tmp, r = rf, Kvs = Kvs_tmp) 

# Figure on page 13

#assign the parameters firstly
K = 110  #strike price
mu = 0  #mean
sig = 0.3  #sigma/annulaized volitility
S0 = 100  #start stock price
r = 0  #interest rate
q = 0  #dividends of stock
dt = 1/252  #time step
T_m = 1  #maturity time
type = "c"  #option type

Path <- function(S0_ = S0, mu_ = mu, dt_ = dt, T_m_ = T_m, seed_) {
  # set random seed
  set.seed(seed_)
  # lognormal random number
  sig1 <- 0.2
  rand_vec1 <- rlnorm(n=200, meanlog = mu_ - 0.5*((sig1*sqrt(1/252))^2), sdlog = sig1*sqrt(1/252))
  sig2 <- 0.3
  rand_vec2 <- rlnorm(n=10, meanlog = mu_ - 0.5*((sig2*sqrt(1/252))^2), sdlog = sig2*sqrt(1/252))
  sig3 <- 0.4
  rand_vec3 <- rlnorm(n=41, meanlog = mu_ - 0.5*((sig3*sqrt(1/252))^2), sdlog = sig3*sqrt(1/252))
  # return and price
  rt_vec <- c(1, rand_vec1, rand_vec2, rand_vec3)
  price_vec <- S0_ * cumprod(rt_vec)
  # return
  return(price_vec)
}
# path test
path <- Path(seed_ = 5)
plot(path, type="l")
abline(h=110, col="red")

#define a payoff function for vanilla options
Vanilla.Payoff <- function(K_ = K, S_end, type = "c")
{
  if(type == "c")  #Call
  {
    payoff <- max(0, (S_end - K_))
  }
  else  #Put
  {
    payoff <- max(0, (K_ - S_end))
  }
  payoff
}

Delta.Hedge.PnL <- function(PayoffFunc, path, K_ = K, S0_ = S0, sig_ = sig, sig_DH = sig, rf=r,
                            mu_ = mu, dt_ = dt, T_m_ = T_m, type = "c") {
  #decide the option type
  if(type == "c")  #Call
  {
    # call option price
    CallPrice <- BSPricer(S=S0_, K=K_, r=rf, Ttm=T_m_, vol=sig_, type="c")
    pnlfunc <- function(x) {
      # calculate option P&L
      S_end <- x[length(x)]
      Opt_payoff <- PayoffFunc(K_ = K_, S_end, type = type)
      
      # calculate delta hedge P&L
      S_path <- ts(x)
      S_diff <- lag(S_path) - S_path
      t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
      delta_vec <- pnorm( (1 / (sig_DH * sqrt(t_ts))) * (log(S_path / K_) + (0 + sig_DH^2/2) * t_ts) )
      delta_vec <- delta_vec[-length(delta_vec)]
      deltapayoff_vec <- delta_vec * S_diff
      #deltapayoff <- sum(deltapayoff_vec)
      
      # return
      cumsum(deltapayoff_vec) + CallPrice - Opt_payoff
    }
    pnl <- pnlfunc(path)
  }
  else  #Put
  {
    #do nothing for now
  }
  pnl
}
# test
pnl <- Delta.Hedge.PnL(PayoffFunc=Vanilla.Payoff, path=path)*100000
plot(pnl, type="l")
abline(h=0, col="blue")
# realized vol
rl_vol <- sd(diff(log(path), 1))*sqrt(252)
rl_50dvol <- sapply(c(51:252), function(x) {
  sd(diff(log(path[(x-50):x]), 1) * sqrt(252))
})

# gamma
gamma <- path^2 * BSGamma(S=path, K=K, r=r, Ttm=rev(seq(from=1/252, to=1, by=1/252)), vol=sig)
plot(gamma, type="l")

# graph (a)
figa_data <- cbind(path, rep(K, 252), c(1,pnl))



# Goldman Sachs Model Pricer
FairVariance <- function(S, put_strikes, call_strikes, put_IVs, call_IVs, Ttm=1, r, ref_S) {
  # check the inputs 
  if (length(put_strikes) != length(put_IVs)) {
    stop("Put strikes are not as long as put IVs, please check!")
  }
  if (length(call_strikes) != length(call_IVs)) {
    stop("Call strikes are not as long as Call IVs, please check!")
  }
  
  # compute weights firstly
  f <- function(ST) (2/Ttm)*((ST - ref_S)/ref_S - log(ST/ref_S))
  # for call weights
  call_weights <- c()
  call_weights[1] <- (f(call_strikes[2]) - f(call_strikes[1])) / (call_strikes[2] - call_strikes[1])
  n_tmp <- length(call_strikes)
  for (i in 2:(n_tmp - 1)) {
    call_weights[i] <- (f(call_strikes[i+1]) - f(call_strikes[i])) / (call_strikes[i+1] - call_strikes[i]) - sum(call_weights[1:(i-1)])
  }
  k_tmp <- call_strikes[n_tmp] + (call_strikes[n_tmp] - call_strikes[n_tmp - 1])
  call_weights[n_tmp] <- (f(k_tmp) - f(call_strikes[n_tmp])) / (k_tmp - call_strikes[n_tmp]) - sum(call_weights[1:(n_tmp-1)])
  # for put weights
  put_weights <- c()
  put_weights[1] <- (f(put_strikes[2]) - f(put_strikes[1])) / (put_strikes[1] - put_strikes[2])
  n_tmp <- length(put_strikes)
  for (i in 2:(n_tmp - 1)) {
    put_weights[i] <- (f(put_strikes[i+1]) - f(put_strikes[i])) / (put_strikes[i] - put_strikes[i+1]) - sum(put_weights[1:(i-1)])
  }
  k_tmp <- put_strikes[n_tmp] + (put_strikes[n_tmp] - put_strikes[n_tmp - 1])
  put_weights[n_tmp] <- (f(k_tmp) - f(put_strikes[n_tmp])) / (put_strikes[n_tmp] - k_tmp) - sum(put_weights[1:(n_tmp-1)])
  
  # compute call and put price series
  put_series <- apply(cbind(put_strikes, put_IVs), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "p")
  })
  call_series <- apply(cbind(call_strikes, call_IVs), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "c")
  })
  
  # sum results
  sum_cp <- sum(put_weights * put_series) + sum(call_weights * call_series)
  
  # price formula
  kvar <- (2/Ttm) * (r*Ttm - (S*exp(r*Ttm)/ref_S - 1) - log(ref_S/S)) + exp(r*Ttm)*sum_cp
  
  # return 
  return(kvar)
  #return(sum_cp)
  #return(list(pweight=put_weights, pprice=put_series, cweight=call_weights, cprice=call_series))
}

# # test GS paper P21 example
# put_strikes_tmp <- rev(seq(from=50, to=100, by=5))
# call_strikes_tmp <- seq(from=100, to=135, by=5)
# put_vol <- rev(seq(from=30, to=20) * 0.01)
# call_vol <- seq(from=20, to=13) * 0.01
# rf <- 0.05
# S0 <- 100
# T_m <- 90/365
# sqrt(10000*FairVariance(S=S0, put_strikes=put_strikes_tmp, call_strikes=call_strikes_tmp, 
#                              put_IVs=put_vol, call_IVs=call_vol, Ttm=T_m, r=rf, ref_S=100))

# 
# # test JPMorgan Paper P17 example
# put_strikes_new <- rev(c(1467.51, 1614.26, 1761.01, 1907.76, 2054.51, 2201.26, 2348.01, 2494.76, 2641.51, 2788.26, 2935.02))
# call_strikes_new <- c(2935.02, 3081.77, 3228.52, 3375.27, 3522.02, 3668.77, 3815.52, 3962.27, 4109.02, 4255.77, 4402.52)
# put_IVs_new <- rev(c(0.276, 0.264, 0.252, 0.240, 0.227, 0.214, 0.200, 0.187, 0.173, 0.160, 0.148))
# call_IVs_new <- c(0.148, 0.137, 0.129, 0.122, 0.119, 0.118, 0.119, 0.121, 0.125, 0.129, 0.134)
# discount_factor <- 0.977368853
# rf <- log(1/discount_factor)
# T_m <- 1
# S0 <- 2935.02 * exp(-rf*T_m)
# sqrt(FairVariance(S=S0, put_strikes = put_strikes_new, call_strikes = call_strikes_new,
#                              put_IVs=put_IVs_new, call_IVs = call_IVs_new, Ttm = T_m, r=rf, ref_S = 2935.02))

# Theoritical Fair Variance
TheoryFairVar.IVlinearK <- function(S, sigma0=0.3, b_slope=0.2, Ttm=1, r, ref_S, K_lower, K_upper, step=0.001) {
  # compute the Sf price
  Sf <- S * exp(r*Ttm)
  # obtain the upper bound of the call strikes
  #   if(b_slope == 0) {
  #     upper <- 200
  #   }
  #   else {
  #     upper <- floor(sigma0*Sf/b_slope + Sf)
  #   }
  put_bound <- c(K_lower, ref_S)
  call_bound <- c(ref_S, K_upper)
  # cut the put and call range
  put_range <- seq(from=put_bound[1], to=put_bound[2], by=step)
  call_range <- seq(from=call_bound[1], to=call_bound[2], by=step)
  # integral
  dk_put <- step
  dk_call <- step
  iv_put <- rep(sigma0, length(put_range)) - b_slope*(put_range - Sf)/Sf
  iv_call <- rep(sigma0, length(call_range)) - b_slope*(call_range - Sf)/Sf
  price_put <- apply(cbind(put_range, iv_put), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "p")
  })
  price_call <- apply(cbind(call_range, iv_call), 1, function(x) {
    BSPricer(S = S, K = x[1], r = r, Ttm = Ttm, vol = x[2], type = "c")
  })
  put_integral <- sum((1/put_range^2) * price_put * dk_put)
  call_integral <- sum((1/call_range^2) * price_call * dk_call)
  # fair variance 
  kvar <- (2/Ttm) * (r*Ttm - (S*exp(r*Ttm)/ref_S - 1) - log(ref_S/S) + exp(r*Ttm)*(put_integral+call_integral))
  # return result
  return(kvar)
}
# test
sqrt(10000*TheoryFairVar.IVlinearK(S=100, sigma0=0.3, b_slope=0.2, Ttm=90/365, r=0.05, ref_S=100, K_lower=10, K_upper=200, step=0.0001))
# # test GS paper P24 example
# b_slope <- 0.2
# sigma0 <- 0.3
# r <- 0.05
# S0 <- 100
# T_m <- 90/365
# Sf <- S0 * exp(r*Ttm)
# put_strikes_tmp <- rev(seq(from=10, to=100, by=1))
# call_strikes_tmp <- seq(from=100, to=200, by=1)
# put_vol <- rep(sigma0, length(put_strikes_tmp)) - b_slope*(put_strikes_tmp - Sf)/Sf
# call_vol <- rep(sigma0, length(call_strikes_tmp)) - b_slope*(call_strikes_tmp - Sf)/Sf
# sqrt(10000*FairVariance(S=S0, put_strikes=put_strikes_tmp, call_strikes=call_strikes_tmp, 
#                              put_IVs=put_vol, call_IVs=call_vol, Ttm=T_m, r=rf, ref_S=100))

# the first approx pricing model
Approx.IVlinearK <- function(S, sigma0=0.3, b_slope=0.2, Ttm=1, r, ref_S) {
  kvar <- sigma0^2 * (1 + 3 * Ttm * b_slope^2)
  return(kvar)
}
# test
sqrt(10000*Approx.IVlinearK(S=100, sigma0=0.3, b_slope=0.2, Ttm=1))

# the second approx pricing model
Approx.IVlinearDelta <- function(S, sigma0=0.3, b_slope=0.2, Ttm=1, r, ref_S) {
  kvar <- sigma0^2 * (1 + (b_slope*sqrt(Ttm)/sqrt(pi)) + (1/12)*(b_slope^2/sigma0^2))
  return(kvar)
}
# test
sqrt(10000*Approx.IVlinearDelta(S=100, sigma0=0.3, b_slope=0.2, Ttm=0.25))