library(ggplot2)
library(xts)
library(portfolioBacktest)
library(riskParityPortfolio)
library(reshape2)
library(quadprog)
#load data online
faang_data <- stockDataDownload(c('AAPL','MSFT','AMZN','FB','GOOGL','GOOG','BRK-B',
                                  'JNJ','JPM','V'),from = "2014-01-01", to = "2020-12-13")

#vanilla RP vs markowitz tengency portfolio

#vanilla risk parity portfolio
risk_parity <- function(dataset) {
  prices <- dataset$adjusted
  log_returns <- diff(log(prices))[-1]
  return(riskParityPortfolio(cov(log_returns))$w)
}

#rp_return
rp_return <- function(dataset, lmd_mu = 0.01) {
  prices <- dataset$adjusted
  log_returns <- diff(log(prices))[-1]
  mu <- colMeans(log_returns)
  return(riskParityPortfolio(cov(log_returns),mu = mu, lmd_mu = lmd_mu, formulation = "rc-over-sd vs b-times-sd")$w)
}

#rp_var
rp_var <- function(dataset, lmd_var = 500) {
  prices <- dataset$adjusted
  log_returns <- diff(log(prices))[-1]
  return(riskParityPortfolio(cov(log_returns), lmd_var = lmd_var, formulation = "rc-over-sd vs b-times-sd")$w)
}

#rp_BC
rp_BC <- function(dataset) {
  prices <- dataset$adjusted
  log_returns <- diff(log(prices))[-1]
  Dmat <- matrix(0, 2, 10)
  Dmat[1, ] <- c(-1,0,-1,-1,0,0,0,0,0,0)
  Dmat[2, ] <- c(0,0,0,0,1,1,0,0,0,0)
  dvec <- c(-0.5, 0.1)
  return(riskParityPortfolio(cov(log_returns), Dmat = Dmat, dvec = dvec)$w)
}

# tangency portfolio (maximum sharpe ratio)
max_sharpe_ratio <- function(dataset) {
  prices <- dataset$adjusted
  log_returns <- diff(log(prices))[-1]
  N <- ncol(prices)
  Sigma <- cov(log_returns)
  mu <- colMeans(log_returns)
  if (all(mu <= 1e-8))
    return(rep(0, N))
  Dmat <- 2 * Sigma
  Amat <- diag(N)
  Amat <- cbind(mu, Amat)
  bvec <- c(1, rep(0, N))
  dvec <- rep(0, N)
  res <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  w <- res$solution
  return(w/sum(w))
}

#call portfolioBacktest and benchmark against the uniform (1/N) portfolio
#rollng: 240 #rebalance: every 40 trading day
bt <- portfolioBacktest(list("risk parity portfolio" = risk_parity,
                             "markowitz portfolio"    = max_sharpe_ratio),
                        list(faang_data),
                        T_rolling_window = 12*20, 
                        optimize_every = 2*20, rebalance_every = 2*20)

#check performance summary
backtestSummary(bt)$performance

# plot cumulative returns chart
backtestChartCumReturns(bt)

#plot max drawdown chart
backtestChartDrawdown(bt)

#plot assets exposures over time
backtestChartStackedBar(bt, portfolio = "risk parity portfolio", legend = TRUE)
backtestChartStackedBar(bt, portfolio = "markowitz portfolio" , legend = TRUE)

#back test of diff RP methods
bt1 <- portfolioBacktest(list("RP" = risk_parity,
                              "RP_return" = rp_return,
                              "RP_var" = rp_var,
                              "RP_BC" = rp_BC),
                         list(faang_data),
                         T_rolling_window = 12*20, 
                         optimize_every = 2*20, rebalance_every = 2*20)

#check performance summary
backtestSummary(bt1)$performance

#plot cumulative returns chart
backtestChartCumReturns(bt1)

#plot max drawdown chart
backtestChartDrawdown(bt1)

#risk parity VS S&P500
rp_sp <- function(dataset) {
  prices <- dataset$adjusted[, -11]
  log_returns <- diff(log(prices))[-1]
  w <- riskParityPortfolio(cov(log_returns))$w
  w <- c(w, 0)
  return(w)
}

sp_500 <- function(dataset){
  return(c(0,0,0,0,0,0,0,0,0,0,1))
}

#load sp500 data
bt_data <- stockDataDownload(c('AAPL','MSFT','AMZN','FB','GOOGL','GOOG','BRK-B',
                               'JNJ','JPM','V','^GSPC'),from = "2014-01-01", to = "2020-12-13")

#back test for RP and SP500
bt2 <- portfolioBacktest(list("RP" = rp_sp,
                              "S&P500" = sp_500),
                         list(bt_data),
                         T_rolling_window = 12*20, 
                         optimize_every = 2*20, rebalance_every = 2*20)

#check performance summary
backtestSummary(bt2)$performance

#plot cumulative returns chart
backtestChartCumReturns(bt2)

#plot max drawdown chart
backtestChartDrawdown(bt2)
