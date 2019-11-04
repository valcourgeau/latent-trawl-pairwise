
acf_trawl <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp', cov=F){
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  trawl_fct <- GetTrawlFunctions(type)
  B1_func <- trawl_fct[[1]]
  B2_func <- trawl_fct[[2]]
  B3_func <- trawl_fct[[3]]
  
  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  
  res <- 0
  first_mom <- 0
  res_0 <- 0
  sum_over_x <- vapply(seq_kappa, FUN = function(x){
    x <- x + delta / 2
    sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
                    y <- y + delta / 2
                    tmp1 <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
                    tmp2 <- 0.5*(1+(x/2+y/2)/beta)^{-alpha}
                    return(c(tmp1, tmp2))
                  },
                  FUN.VALUE = rep(0, 2))
    return(apply(sum_over_y, MARGIN = 1, sum))},
    FUN.VALUE = rep(0, 2))
  
  final_sum <- apply(sum_over_x, MARGIN = 1, sum)
  final_sum <- final_sum * delta^2
  
  res <- final_sum[1]
  res_0 <- final_sum[2]
  first_mom_sq <-  ((1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha - 1))^2
  # if(h < 1){
  #   cat('res', res, '\n')
  #   cat('res_0', res_0, '\n')
  #   print(first_mom_sq)
  # }
  
  if(cov){
    return(res-first_mom_sq)
  }else{
    return((res-first_mom_sq)/(res_0-first_mom_sq))
  }
}

crossmoment_trawls <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp'){
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  trawl_fct <- GetTrawlFunctions(type)
  B1_func <- trawl_fct[[1]]
  B2_func <- trawl_fct[[2]]
  B3_func <- trawl_fct[[3]]
  
  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))

  res <- 0
  first_mom <- 0
  res_0 <- 0
  beta <- beta
  
  sum_over_x <- vapply(seq_kappa, FUN = function(x){
    x <- x + delta / 2
    sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
      y <- y + delta / 2
      tmp1 <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      return(tmp1)
    },
    FUN.VALUE = rep(0, 1))
    return(sum(sum_over_y))},
    FUN.VALUE = rep(0, 1))
  
  final_sum <- sum(sum_over_x)
  res <- final_sum * delta^2
  return(res)
}

acf_trawl_num_approx <- function(h, alpha, beta, kappa, rho, delta=0.5, type='exp', cov=T){
  vapply(h, function(h){
    acf_trawl(h, alpha = alpha, beta = beta, kappa = kappa, 
              rho = rho, delta = delta, type = type, cov = cov)}, 1)}

acf_trawl_inv <- function(h, alpha, beta, rho, kappa, delta = 0.1, end_seq = 50, type='exp', cov=T){
  seq_kappa <- seq(kappa, kappa+end_seq, by = delta)
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type)
  b_h_minus_0 <- - alpha  * B1_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_minus_h <- - alpha  * B3_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  b_0_h <- - alpha * B2_func(param=rho, h=h)/(B1_func(param=rho, h=h) + B2_func(param=rho, h=h))
  res <- 0
  first_mom <- 0
  res_0 <- 0
  beta <- beta
  
  sum_over_x <- vapply(seq_kappa, FUN = function(x){
    x <- x + delta / 2
    # first_mom <- first_mom + (1+x/beta)^{-alpha}*(1+y/beta)^{-alpha}
    sum_over_y <- vapply(X = seq_kappa, FUN = function(y){
      y <- y + delta / 2
      tmp_11 <- (1+(x-kappa)/beta)^{b_h_minus_0} * (1+(x+y-2*kappa)/beta)^{b_h_minus_0} * (1+(y-kappa)/beta)^{b_h_minus_0}
      tmp_1e <- 2*(1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+(y-kappa)/beta)^{b_0_h}
      tmp_ee <- (1+x/beta)^{b_h_minus_0} * (1+(x+y)/beta)^{b_0_h} * (1+y/beta)^{b_0_minus_h}
      return(c(tmp_11, tmp_1e, tmp_ee))
    },
    FUN.VALUE = rep(0, 3))
    return(c(apply(sum_over_y, MARGIN = 1, sum), (1+(x)/beta)^{-alpha}))},
    FUN.VALUE = rep(0, 4))
  
  # prob_xt_and_xs_pos <- 1 + (1+kappa/beta)^{b_0_minus_h}*(1+2*kappa/beta)^{b_0_h}*(1+kappa/beta)^{b_h_minus_0}
  # prob_xt_and_xs_pos <- prob_xt_and_xs_pos - 2 * (1+kappa/beta)^{-alpha}
  
  final_sum <- apply(sum_over_x, MARGIN = 1, sum) 
  
  res <- (final_sum[1] - final_sum[2] + final_sum[3]) * delta^2
  res_0 <- sum(sum_over_x[4,]^2) * delta
  
  # first_mom_sq <-  ((1+kappa/beta)^{-alpha}*(beta+kappa)/(alpha - 1))^2
  first_mom_sq <- final_sum[4]^2 * delta^2
  print(res)
  print(res_0)
  print(first_mom_sq)
  if(cov){
    return(res-first_mom_sq)
  }else{
    return((res-first_mom_sq)/(res_0-first_mom_sq))
  }
}

acf_trawl_num_approx_inv <- function(h, alpha, beta, kappa, rho, delta=0.5, type='exp', cov=T){
  vapply(h, function(h){
    acf_trawl_inv(h, alpha = alpha, beta = beta, kappa = kappa, 
              rho = rho, delta = delta, type = type, cov = cov)}, 1)}


DupuisSimplified <- function(data_u, n_trials=10, acf_depth=15, mult_fac=c(0.3, 3), cl=NULL){
  params <- rep(0, 6)
  custom_mle_results <- CustomMarginalMLE(data = data_u)
  kappa <- GetKappa(data = data_u, params = custom_mle_results, parametrisation = 'standard')
  params[1:3] <- c(custom_mle_results, kappa)
  params[5:6] <- ParametrisationTranslator(
    params = params,
    parametrisation = 'standard',
    target = 'noven')[1:2]
  print(params)
  depth <- acf_depth
  kk <- acf(as.numeric(data_u>0), lag.max = depth-1, plot=F)
  
  alpha_tmp <- params[5]
  beta_tmp <- params[6]
  kappa_tmp <- params[3]
  
  mae_tab <- rep(0, n_trials)
  mse_tab <- rep(0, n_trials)
  # rho_tab <- seq(0.01, (3), length.out = n_trials)
  index <- 1
  
  lin_rho <- line(x = c(0, 1:(depth-2)), log(kk$acf[1:(depth-1)] %>% abs))
  rho_tmp <- abs(lin_rho$coefficients[2])
  print(rho_tmp)
  rho_tab <- seq(log(rho_tmp*mult_fac[1]), log(rho_tmp*mult_fac[2]), length.out = n_trials) %>% exp
  # rho_tab <-seq(log(0.01), log(2.0), length.out = n_trials) %>% exp
  
  for(rho_iter in rho_tab){
    if(!is.null(cl)){
      parallel::clusterExport(cl, c('acf_trawl', 'alpha_tmp', 'beta_tmp', 'kappa_tmp', 'rho_iter'))
      acf_vals <- parallel::parLapply(cl, X = c(0.05, 1:(depth-1)), fun = function(h){
        return(acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                         rho = rho_iter, delta = 0.5, end_seq = 50))})
      acf_vals <- unlist(acf_vals)
    }else{
      acf_vals <- vapply(c(0.05, 1:(depth-1)), function(h){
        acf_trawl(h, alpha = alpha_tmp, beta = beta_tmp, kappa = kappa_tmp, 
                  rho = rho_iter, delta = 0.5, end_seq = 50)}, 1)
    }
    
    # plot(kk$acf)
    plot(acf_vals, col='red')
    mae_tab[index] <- sum(abs(kk$acf - acf_vals)) 
    mse_tab[index] <- sum(((kk$acf - acf_vals)^2))
    index <- index + 1
  }
  
  plot(exp(rho_tab), mse_tab)
  points(exp(rho_tab), mae_tab)
  
  params[4] <- rho_tab[which.min(mse_tab)]
  names(params) <- c('xi', 'sigma', 'kappa', 'rho', 'alpha', 'beta')
  return(params)
}
