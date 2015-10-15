## Derivatives Pricers
## Author: Ming Tian
## Date: 2015-05-06
## Copyright @2015 MingTian

#assign the parameters firstly
K = 100  #strike price
mu = 0  #mean
sig = 0.2  #sigma/annulaized volitility
S0 = 100  #start stock price
r = 0  #interest rate
q = 0  #dividends of stock
dt = 1/252  #time step
T_m = 1  #maturity time
type = "c"  #option type

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
# # test
# BSPricer(100, 100, 0.05, 0.25, 0.2, "c")
# library(VarSwapPrice)
# black_scholes(100, 100, 0.05, 0.25, 0.2)
# library(fOptions)
# GBSOption(TypeFlag="c", S=100, X=100, Time=0.25, r=0.05, b=0.05, sigma=0.2)


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
# # this is the new test
# test <- Random.Path.New(n_path = 10000)
# hist(test[,ncol(test)], breaks=50)


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

# 20150917 updated
Delta.Hedge.Return.Vectorization <- function(PayoffFunc, K_ = K, S0_ = S0, sig_ = sig, sig_DH = sig, 
                                             mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000, type = "c")
{
  #Call the Random.Path function firstly
  rand_path <- Random.Path.New(n_path = n_path_)
  #build a return result container
  result <- c(1:n_path_)
  #decide the option type
  if(type == "c")  #Call
  {
    # call option price
    d1 <- (1 / (sig_ * sqrt(T_m_))) * (log(S0_ / K_) + (0 + sig_^2/2) * T_m_)
    d2 <- d1 - sig_ * sqrt(T_m_)
    CallPrice <- S0_ * pnorm(d1) - K_ * exp(0) * pnorm(d2)
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
      deltapayoff_vec <- -delta_vec * S_diff
      deltapayoff <- sum(deltapayoff_vec)
      
      # return
      (Opt_payoff + deltapayoff - CallPrice)
    }
    
    result <- apply(rand_path, 1, pnlfunc)
  }
  else  #Put
  {
    #do nothing for now
  }
  #return
  result
}
# # this is new Vectorization test
# test <- Delta.Hedge.Return.Vectorization(PayoffFunc = Vanilla.Payoff)
# hist(test, breaks=50)

# Hedge to Market or Hedge to model Comparison
#modify the function above
Market.vs.Model <- function(PayoffFunc, K_ = K, S0_ = S0, sig_ = sig, sig_DH = sig, 
                            mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000, type = "c")
{
  #Call the Random.Path function firstly
  rand_path <- Random.Path.New(n_path = n_path_, sig_ = sig_)
  #build a return result container
  result <- c(1:n_path_)
  #decide the option type
  if(type == "c")  #Call
  {
    # call option price
    d1 <- (1 / (sig_ * sqrt(T_m_))) * (log(S0_ / K_) + (0 + sig_^2/2) * T_m_)
    d2 <- d1 - sig_ * sqrt(T_m_)
    CallPrice <- S0_ * pnorm(d1) - K_ * exp(0) * pnorm(d2)
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
      deltapayoff_vec <- -delta_vec * S_diff
      deltapayoff <- sum(deltapayoff_vec)
      
      # return
      (Opt_payoff + deltapayoff - CallPrice)
    }
    
    result <- apply(rand_path, 1, pnlfunc)
  }
  else  #Put
  {
    #do nothing for now
  }
  #result
  PnL <- result
  lastp <- rand_path[,ncol(rand_path)]
  func_vol <- function(vec)
  {
    vec <- ts(vec)
    vec_l <- lag(vec)
    vol <- sd(log(vec_l/vec))
    return(vol*sqrt(252))
  }
  vol <- apply(rand_path, 1, func_vol)
  result <- cbind(PnL, vol, lastp)
  return(result)
}

#build the four models
Model1 <- Market.vs.Model(PayoffFunc = Vanilla.Payoff, sig_ = 0.1, sig_DH = 0.1)
Model2 <- Market.vs.Model(PayoffFunc = Vanilla.Payoff, sig_ = 0.1, sig_DH = 0.2)
Model3 <- Market.vs.Model(PayoffFunc = Vanilla.Payoff, sig_ = 0.2, sig_DH = 0.1)
Model4 <- Market.vs.Model(PayoffFunc = Vanilla.Payoff, sig_ = 0.2, sig_DH = 0.2)
#12 figures arranged in 4 rows and 3 columns
par(mar=c(1,1,1,1))
par(mfrow=c(4,4))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Simulated at 10%,\n
     hedge at 10%", cex=1.2)
hist(Model1[,1], breaks=60 ,main = "")
plot(x=Model1[,3], y=Model1[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model1[,2], y=Model1[,1], col = "blue", abline(lm(Model1[,1]~Model1[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Simulated at 10%,\n
     hedge at 20%", cex=1.2)
hist(Model2[,1], breaks=60 , main = "")
plot(x=Model2[,3], y=Model2[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model2[,2], y=Model2[,1], col = "blue", abline(lm(Model2[,1]~Model2[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Simulated at 20%,\n
     hedge at 10%", cex=1.2)
hist(Model3[,1], breaks=60 , main = "")
plot(x=Model3[,3], y=Model3[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model3[,2], y=Model3[,1], col = "blue", abline(lm(Model3[,1]~Model3[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Simulated at 20%,\n
     hedge at 20%", cex=1.2)
hist(Model4[,1], breaks=60 , main = "")
plot(x=Model4[,3], y=Model4[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model4[,2], y=Model4[,1], col = "blue", abline(lm(Model4[,1]~Model4[,2])), cex=0.5)


# A new question about delta hedge as bellow:
"
Should we hedge to model if OOTM and market if ITM? And average if ATM/NTM?

Or maybe market if it is far OOTM or far ITM and model if NTM (near the money)?

Make reasonable assumptions and answer this question.

(At a minimum make the three graphs as before, for whatever you scenario is.)
"

Market.vs.Model1 <- function(PayoffFunc, K_ = K, S0_ = S0, sig_ = sig, sig_DH = sig, 
                             mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000, type = "c", level = 20)
{
  #Call the Random.Path function firstly
  rand_path <- Random.Path.New(n_path = n_path_, sig_ = sig_)
  #build a return result container
  result <- c(1:n_path_)
  #decide the option type
  if(type == "c")  #Call
  {
    # call option price
    d1 <- (1 / (sig_ * sqrt(T_m_))) * (log(S0_ / K_) + (0 + sig_^2/2) * T_m_)
    d2 <- d1 - sig_ * sqrt(T_m_)
    CallPrice <- S0_ * pnorm(d1) - K_ * exp(0) * pnorm(d2)
    
    # modify the function for dynamic hedgeing
    pnlfunc <- function(x) {
      # calculate option payoff
      S_end <- x[length(x)]
      Opt_payoff <- PayoffFunc(K_ = K_, S_end, type = type)
      
      # generate a new sigma vector
      sig_vec <- rep(0, n=length(x))
      sig_vec[x > (K_ + level) | x < (K_ - level)] <- sig_DH
      sig_vec[x >= (K_ - level) & x <= (K_ + level)] <- sig_
      #sig_vec[x == K_] <- mean(c(sig_, sig_DH))
      
      # calculate delta hedge P&L
      S_path <- ts(x)
      S_diff <- lag(S_path) - S_path
      t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
      delta_vec <- pnorm( (1 / (sig_vec * sqrt(t_ts))) * (log(S_path / K_) + (0 + sig_vec^2/2) * t_ts) )
      delta_vec <- delta_vec[-length(delta_vec)]
      deltapayoff_vec <- -delta_vec * S_diff
      deltapayoff <- sum(deltapayoff_vec)
      
      # return
      (Opt_payoff + deltapayoff - CallPrice)
    }
    
    result <- apply(rand_path, 1, pnlfunc)
  }
  else  #Put
  {
    #do nothing for now
  }
  #result
  PnL <- result
  lastp <- rand_path[,ncol(rand_path)]
  func_vol <- function(vec)
  {
    vec <- ts(vec)
    vec_l <- lag(vec)
    vol <- sd(log(vec_l/vec))
    return(vol*sqrt(252))
  }
  vol <- apply(rand_path, 1, func_vol)
  result <- cbind(PnL, vol, lastp)
  return(result)
}

# # test function
# Modeltest <- Market.vs.Model1(PayoffFunc = Vanilla.Payoff, sig_ = 0.1, sig_DH = 0.1)
# hist(Modeltest[,1], breaks=60 ,main = "")
# plot(x=Modeltest[,3], y=Modeltest[,1], col = "blue", abline(h=0), cex=0.5)
# plot(x=Modeltest[,2], y=Modeltest[,1], col = "blue", abline(lm(Modeltest[,1]~Modeltest[,2])), cex=0.5)

#build the four models
Model1 <- Market.vs.Model1(PayoffFunc = Vanilla.Payoff, sig_ = 0.1, sig_DH = 0.1, level = 10)
Model2 <- Market.vs.Model1(PayoffFunc = Vanilla.Payoff, sig_ = 0.1, sig_DH = 0.2, level = 10)
Model3 <- Market.vs.Model1(PayoffFunc = Vanilla.Payoff, sig_ = 0.2, sig_DH = 0.1, level = 10)
Model4 <- Market.vs.Model1(PayoffFunc = Vanilla.Payoff, sig_ = 0.2, sig_DH = 0.2, level = 10)
#12 figures arranged in 4 rows and 3 columns
par(mar=c(1,1,1,1))
par(mfrow=c(4,4))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Market at 10%,\n
     Model at 10%", cex=1.2)
hist(Model1[,1], breaks=60 ,main = "")
plot(x=Model1[,3], y=Model1[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model1[,2], y=Model1[,1], col = "blue", abline(lm(Model1[,1]~Model1[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Market at 10%,\n
     Model at 20%", cex=1.2)
hist(Model2[,1], breaks=60 , main = "")
plot(x=Model2[,3], y=Model2[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model2[,2], y=Model2[,1], col = "blue", abline(lm(Model2[,1]~Model2[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Market at 20%,\n
     Model at 10%", cex=1.2)
hist(Model3[,1], breaks=60 , main = "")
plot(x=Model3[,3], y=Model3[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model3[,2], y=Model3[,1], col = "blue", abline(lm(Model3[,1]~Model3[,2])), cex=0.5)
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x=0.5,y=0.5, labels = "Market at 20%,\n
     Model at 20%", cex=1.2)
hist(Model4[,1], breaks=60 , main = "")
plot(x=Model4[,3], y=Model4[,1], col = "blue", abline(h=0), cex=0.5)
plot(x=Model4[,2], y=Model4[,1], col = "blue", abline(lm(Model4[,1]~Model4[,2])), cex=0.5)

# Why it is impossible to arbitrage from volatility skew
"figure out why selling 70% puts at 30 vol, buying 130% calls at 10 vol, 
and hedging both to 20 vol, won't necessarily make 20 vol points of profit"

# "best trade" function - the false one
LongCallSellPutStart <- function(K_call = 130, K_put = 70, S0_ = 100, sig_sim = 0.2, sig_call = 0.1, sig_put = 0.3, sig_DH = 0.2,
                                 mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000) {
  #Call the Random.Path function firstly
  rand_path <- Random.Path.New(n_path = n_path_, S0_ = S0_, sig_ = sig_sim, mu_ = mu_, dt_ = dt_, T_m_ = T_m_)
  
  # long call part
  # call option price
  d1 <- (1 / (sig_call * sqrt(T_m_))) * (log(S0_ / K_call) + (0 + sig_call^2/2) * T_m_)
  d2 <- d1 - sig_call * sqrt(T_m_)
  CallPrice <- S0_ * pnorm(d1) - K_call * exp(0) * pnorm(d2)
  
  # delta hedge function
  pnlfunc <- function(x) {
    # calculate option P&L
    S_end <- x[length(x)]
    Opt_payoff <- Vanilla.Payoff(K_ = K_call, S_end, type = "c")
    
    # calculate delta hedge P&L
    S_path <- ts(x)
    S_diff <- lag(S_path) - S_path
    t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
    delta_vec <- pnorm( (1 / (sig_DH * sqrt(t_ts))) * (log(S_path / K_call) + (0 + sig_DH^2/2) * t_ts) )
    delta_vec <- delta_vec[-length(delta_vec)]
    deltapayoff_vec <- -delta_vec * S_diff
    deltapayoff <- sum(deltapayoff_vec)
    
    # return
    (Opt_payoff + deltapayoff - CallPrice)
  }
  
  cresult <- apply(rand_path, 1, pnlfunc)
  # long call part end
  
  # sell put part
  # put option price
  d1 <- (1 / (sig_put * sqrt(T_m_))) * (log(S0_ / K_put) + (0 + sig_put^2/2) * T_m_)
  d2 <- d1 - sig_put * sqrt(T_m_)
  PutPrice <- K_put * exp(0) * pnorm(-d2) - S0_ * pnorm(-d1)
  
  # delta hedge function
  pnlfunc <- function(x) {
    # calculate option P&L
    S_end <- x[length(x)]
    Opt_payoff <- Vanilla.Payoff(K_ = K_put, S_end, type = "p")
    
    # calculate delta hedge P&L
    S_path <- ts(x)
    S_diff <- lag(S_path) - S_path
    t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
    delta_vec <- pnorm( (1 / (sig_DH * sqrt(t_ts))) * (log(S_path / K_put) + (0 + sig_DH^2/2) * t_ts) ) - 1
    delta_vec <- delta_vec[-length(delta_vec)]
    deltapayoff_vec <- delta_vec * S_diff
    deltapayoff <- sum(deltapayoff_vec)
    
    # return
    (-Opt_payoff + deltapayoff + PutPrice)
  }
  
  presult <- apply(rand_path, 1, pnlfunc)
  # sell put part end
  
  # result
  result <- cresult + presult
  
  # return 
  return(result)
}
# showing result
lcsp <- LongCallSellPutStart(mu_ = 0, dt_ = 1/252, T_m_ = 1, n_path_ = 10000)
hist(lcsp, breaks=50)

# 20151005 update version
# local volitility function - assume linear relationship with price S
GetLocalVol <- function(S, sigma80 = 0.4, sigma120 = 0.05) {
  (S - 80)*(sigma120 - sigma80) / 40 + sigma80
}
# generate the random path from local volitility
Random.Path.LV <- function(S0_ = S0, mu_ = mu, dt_ = dt, T_m_ = T_m, n_path) {
  # generate the random vector
  # rand_vec <- rlnorm(n=(T_m_/dt_)*n_path, meanlog = mu_ - 0.5*((sig_*sqrt(1/252))^2), sdlog = sig_*sqrt(1/252))
  
  # build the result
  result <- matrix(data=0, nrow=n_path, ncol=(T_m_/dt_))
  
  # apply the function to transform to the price ts
  result[, 1] <- S0_
  result_new <- t(apply(result, 1, function(x) {
    #cumprod(x)
    # start from the day1 and enter for loop
    for (i in 2:length(x)) {
      # get the local vol
      lv <- GetLocalVol(S=x[i-1])
      x[i] <- x[i-1] * rlnorm(n=1, meanlog = mu_ - 0.5*(lv * sqrt(1/252))^2, sdlog = lv * sqrt(1/252))
    }
    x
  }))
  # return the result
  result_new
}

# correct "best trade" function - arbitrage is impossible
LongCallSellPutStart.Correct <- function(K_call = 130, K_put = 70, S0_ = 100, sig_sim = 0.2, sig_call = 0.1, sig_put = 0.3, sig_DH = 0.2,
                                         mu_ = mu, dt_ = dt, T_m_ = T_m, n_path_ = 10000) {
  #Call the Random.Path function firstly
  rand_path <- Random.Path.LV(n_path = n_path_, S0_ = S0_, mu_ = mu_, dt_ = dt_, T_m_ = T_m_)
  
  # long call part
  # call option price
  d1 <- (1 / (sig_call * sqrt(T_m_))) * (log(S0_ / K_call) + (0 + sig_call^2/2) * T_m_)
  d2 <- d1 - sig_call * sqrt(T_m_)
  CallPrice <- S0_ * pnorm(d1) - K_call * exp(0) * pnorm(d2)
  
  # delta hedge function
  pnlfunc <- function(x) {
    # calculate option P&L
    S_end <- x[length(x)]
    Opt_payoff <- Vanilla.Payoff(K_ = K_call, S_end, type = "c")
    
    # calculate delta hedge P&L
    S_path <- ts(x)
    S_diff <- lag(S_path) - S_path
    t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
    delta_vec <- pnorm( (1 / (sig_DH * sqrt(t_ts))) * (log(S_path / K_call) + (0 + sig_DH^2/2) * t_ts) )
    delta_vec <- delta_vec[-length(delta_vec)]
    deltapayoff_vec <- -delta_vec * S_diff
    deltapayoff <- sum(deltapayoff_vec)
    
    # return
    (Opt_payoff + deltapayoff - CallPrice)
  }
  
  cresult <- apply(rand_path, 1, pnlfunc)
  # long call part end
  
  # sell put part
  # put option price
  d1 <- (1 / (sig_put * sqrt(T_m_))) * (log(S0_ / K_put) + (0 + sig_put^2/2) * T_m_)
  d2 <- d1 - sig_put * sqrt(T_m_)
  PutPrice <- K_put * exp(0) * pnorm(-d2) - S0_ * pnorm(-d1)
  
  # delta hedge function
  pnlfunc <- function(x) {
    # calculate option P&L
    S_end <- x[length(x)]
    Opt_payoff <- Vanilla.Payoff(K_ = K_put, S_end, type = "p")
    
    # calculate delta hedge P&L
    S_path <- ts(x)
    S_diff <- lag(S_path) - S_path
    t_ts <- rev(seq(from=dt_, to=T_m_, by=dt_))
    delta_vec <- pnorm( (1 / (sig_DH * sqrt(t_ts))) * (log(S_path / K_put) + (0 + sig_DH^2/2) * t_ts) ) - 1
    delta_vec <- delta_vec[-length(delta_vec)]
    deltapayoff_vec <- delta_vec * S_diff
    deltapayoff <- sum(deltapayoff_vec)
    
    # return
    (-Opt_payoff + deltapayoff + PutPrice)
  }
  
  presult <- apply(rand_path, 1, pnlfunc)
  # sell put part end
  
  # result
  result <- cresult + presult
  
  # return 
  return(result)
}
# showing result
lcsp <- LongCallSellPutStart.Correct(mu_ = 0, dt_ = 1/252, T_m_ = 1, n_path_ = 10000)
hist(lcsp, breaks=50)
summary(lcsp)

