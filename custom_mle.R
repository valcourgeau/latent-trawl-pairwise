CustomMarginalMLE <- function(data){
  init_guess <- eva::gpdFit(data, 0.0)$par.sum$Estimate
  
  fn_mle <- function(par){
    data_for_mle <- data
    p <- length(which(data_for_mle>0)) / length(data_for_mle)
    kap <- (1-p^{par[1]})*par[2]/abs(par[1])
    p_non_zero <- 1-(1+kap/(par[2]/abs(par[1])-kap))^{-par[1]}
    like <- eva::dgpd(data_for_mle[data_for_mle>0.0], scale = par[2], shape=par[1],
                      loc=0, log.d = F) * p_non_zero
    
    log.like <- sum(log(like)) + length(data_for_mle == 0.0) * log(1-p_non_zero)
    return(-log.like)
  }
  
  return(stats::optim(par = init_guess[2:1], fn_mle, method='L-BFGS-B', lower=c(1e-03, 0.5), upper=c(2,20))$par)
}
