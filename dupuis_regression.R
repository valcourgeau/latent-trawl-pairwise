
acf_trawl <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp'){
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  res <- 0
  first_mom <- 0
  res_0 <- 0
  beta <- beta
  
  # first_mom <- 0
  
  for(x in seq_kappa){
    x <- x + delta / 2
    # first_mom <- first_mom + (1+x/beta)^{-alpha}*(1+y/beta)^{-alpha}
    for(y in seq_kappa){
      y <- y + delta / 2
      res <- res + (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      res_0 <- res_0 + (1+(x+y)/beta)^{-alpha}
    }
  }
  
  res <- res * delta^2
  res_0 <- res_0 * delta^2
  # first_mom <- first_mom * delta
  first_mom_sq <-  ((1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha - 1))^2
  # print(res/res_0)
  # print(first_mom^2)
  # res <- (res - first_mom^2)^1 # first_mom_sq
  # res_0 <- (res_0 - first_mom^2)^1 # first_mom_sq
  return((res-first_mom_sq)/(res_0-first_mom_sq))
}

acf_trawl_num_approx <- function(h, alpha, beta, kappa, rho, delta=0.5, type='exp'){
  vapply(h, function(h){
    acf_trawl(h, alpha = alpha, beta = beta, kappa = kappa, 
              rho = rho, delta = delta, type = type)}, 1)}


DupuisSimplified <- function(data_u, n_trials=10, acf_depth=15, mult_fac=c(0.3, 3), cl=NULL){
  params <- rep(0, 6)
  custom_mle_results <- CustomMarginalMLE(data = data_u)
  kappa <- GetKappa(data = data_u, params = custom_mle_results, parametrisation = 'standard')
  params[1:3] <- c(custom_mle_results, kappa)
  
  params[5] <- 1/params[1]
  params[6] <- params[2]/abs(params[1]) - params[3]
  
  depth <- acf_depth
  kk <- acf(as.numeric(data_u>0), lag.max = depth-1, plot=F)
  
  # print(params[1])
  alpha_tmp <- params[5]
  beta_tmp <- params[6]
  kappa_tmp <- params[3]
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
  
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  # rho_tab <- seq(0.01, (3), length.out = n_trials)
  index <- 1
  
  lin_rho <- line(x = c(0, 1:(depth-2)), log(kk$acf[1:(depth-1)] %>% abs))
  rho_tmp <- abs(lin_rho$coefficients[2])
  print(rho_tmp)
  rho_tab <- seq(log(rho_tmp*mult_fac[1]), log(rho_tmp*mult_fac[2]), length.out = n_trials) %>% exp
  rho_tab <-seq(log(0.01), log(2.0), length.out = n_trials) %>% exp
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
  
  params[4] <- rho_tab[which.min(mse_tab)]
  names(params) <- c('xi', 'sigma', 'kappa', 'rho', 'alpha', 'beta')
  return(params)
}
