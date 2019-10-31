source('custom_mle.R')
source('utils.R')

library(zeallot)
library(magrittr)

pollution_data <- read.csv('data/clean_pollution_data.csv')

TrawlObjective <- function(data, depth, parametrisation='standard'){
  function(pars){
    kappa <- GetKappa(data = data, params = pars, parametrisation = 'standard')
    pars <- c(pars[1:2], kappa, pars[3:length(pars)])
    noven_pars <- ParametrisationTranslator(params = pars, parametrisation = parametrisation, target = 'noven')
    return(function(trawl_params){
        acf_vals <- vapply(c(0.01, 1:(depth)), function(h){
          crossmoment_trawls(h, alpha = noven_pars[1], beta = noven_pars[2], kappa = noven_pars[3], 
                             rho = trawl_params, delta = 0.5, end_seq = 50)}, 1)
        sample_cross_mom <- acf(data, demean = F, plot = F, type = 'cov', lag.max = depth)$acf
    
        return(sum((acf_vals-sample_cross_mom)))
      }
    )
  }
}
to <- TrawlObjective(data = pollution_data[,2],
               depth = 10,
               parametrisation = 'standard')


GMMObjective <- function(data, depth, omega='id', parametrisation='standard'){
  composite <- CustomLikelihood(data = data,
                                   parametrisation=parametrisation)
  trawl_objective <- TrawlObjective(data = data,
                                    depth = depth,
                                    parametrisation = parametrisation)
  
  return(function(par){
    grad_vec <- c(
      pracma::grad(composite, x0 = par[1:2]),
      pracma::grad(trawl_objective(par), x0 = par[3:length(par)])
    )
    
    if(omega == 'id'){
      omega <- diag(rep(1, length(par)))
    }else{
      
      if(omega == 'centered'){
        omega <- diag(rep(1, length(par)))
      }
    }
    
    return(t(grad_vec) %*% omega %*% grad_vec)
  })
}

custom_mle <- CustomMarginalMLE(pollution_data[,2])
kappa <- GetKappa(pollution_data[,2] ,params = custom_mle, parametrisation = 'standard')
custom_mle_kappa <- c(custom_mle, kappa)
c_mle_kappa_rho <- c(custom_mle_kappa, 0.2)
gmm_obj <- GMMObjective(data = pollution_data[,2], depth = 10)
gmm_obj(c_mle_kappa_rho[-3])

optim(fn = gmm_obj, par = c_mle_kappa_rho[-3], lower=rep(1e-3, 3), upper=rep(3, 2))

kk <- acf(pollution_data[,2], lag.max = depth-1, plot=F)
kk
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
