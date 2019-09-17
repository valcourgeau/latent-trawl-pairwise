
DupuisSimplified <- function(data_u, mult_fac=c(0.3, 3), cl=NULL){
  params <- rep(0, 6)
  params[1:2] <- CustomMarginalMLE(data_u)
  p <- length(which(data_u > 0))/length(data_u)
  params[4] <- (1 - p^{params[1]}) * params[2]/abs(params[1])
  params[5] <- 1.0/params[1]
  params[6] <- params[2]/ abs(params[1]) - params[4]
  
  depth <- 10
  kk <- acf(data_u, lag.max = depth-1, plot=F)
  
  # print(params[1])
  alpha_tmp <- params[5]
  beta_tmp <- params[6]
  kappa_tmp <- params[4]
  # if(params[1] > 0){
  #   alpha_tmp <- params[5]
  #   beta_tmp <- params[6]
  #   kappa_tmp <- params[4]
  # }else{
  #   alpha_tmp <- 2
  #   beta_tmp <- 1
  #   kappa_tmp <- (p^{-1/alpha_tmp} - 1) * beta_tmp
  # }
  
  # cat('alpha', alpha_tmp, '\n')
  # cat('beta', beta_tmp, '\n')
  # cat('kappa', kappa_tmp, '\n')
  # 
  
  n_trials <- 10
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  # rho_tab <- seq(0.01, (3), length.out = n_trials)
  index <- 1
  
  lin_rho <- line(x = c(0, 1:(depth-2)), log(kk$acf[1:(depth-1)] %>% abs))
  rho_tmp <- abs(lin_rho$coefficients[2])
  rho_tab <- seq(log(rho_tmp*mult_fac[1]), log(rho_tmp*mult_fac[2]), length.out = n_trials) %>% exp
  
  # plot(c(0.05, 1:(depth-1)), kk$acf, ylim=c(0,1))
  # print(kk$acf)
  
  for(rho_iter in rho_tab){
    # print(index)
    
    if(!is.null(cl)){
      parallel::clusterExport(cl, c('acf_trawl', 'alpha_tmp', 'beta_tmp', 'kappa_tmp', 'rho_iter'))
      acf_vals <- parallel::parLapply(cl, X = c(0.05, 1:(depth-1)), fun = function(h){
        return(acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                         rho = rho_iter, delta = 0.1, end_seq = 50))})
      acf_vals <- unlist(acf_vals)
    }else{
      acf_vals <- vapply(c(0.05, 1:(depth-1)), function(h){
        acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                  rho = rho_iter, delta = 0.5, end_seq = 50)}, 1)
    }
    
    # plot(x = c(0.05, 1:(depth-1)), log(acf_vals))
    # lines(lin_rho$fitted.values)
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum(((kk$acf - acf_vals)^2)) #+ 
    
    # cat('rho', rho_iter, '\n')
    # cat('MSE', (mse_tab[index]), '\n')
    # cat('MAE', (mae_tab[index]), '\n')
    # 
    # if(index %% 2 == 0){
    #   acf_vals %>% (function(x){lines(c(0.05, 1:(depth-1)), x, col = 'red')})
    # }
    index <- index + 1
  }
  
  params[3] <- rho_tab[which.min(mse_tab)]
  names(params) <- c('xi', 'sigma', 'rho', 'kappa', 'alpha', 'beta')
  return(params)
}
