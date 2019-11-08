TransformationMap <- function(x, params_std, params_trf){
  # from original to trf
  # params_trf <- ParametrisationTranslator(params = params,
  #                                         parametrisation = parametrisation,
  #                                         target = 'transform',
  #                                         target_alpha = target_alpha)
  # params_std <- ParametrisationTranslator(params = params,
  #                                         parametrisation = parametrisation,
  #                                         target = 'standard')
  x[which(x>0.0)] <- eva::qgpd(p = eva::pgpd(q = x[which(x>0.0)],
                                             loc = 0.0,
                                             shape = params_std[1],
                                             scale = params_std[2]),
                               loc = 0.0,
                               shape = params_trf[1],
                               scale = params_trf[2])
  return(x)
}

TransformationMapInverse <- function(x, params_std, params_trf, target_alpha=3){
  # fromt trf to original
  # params_trf <- ParametrisationTranslator(params = params,
  #                                         parametrisation = parametrisation,
  #                                         target = 'transform',
  #                                         target_alpha = target_alpha)
  # params_std <- ParametrisationTranslator(params = params,
  #                                         parametrisation = parametrisation,
  #                                         target = 'standard')
  print(params_trf)
  x[which(x>0.0)] <- eva::qgpd(p = eva::pgpd(q = x[which(x>0.0)],
                                             loc = 0.0,
                                             shape = params_trf[1],
                                             scale = params_trf[2]),
                               loc = 0.0,
                               shape = params_std[1],
                               scale = params_std[2])
  return(x)
}

TransformationJacobian <- function(params_std, params_trf, target_alpha=3){
  #parametrisation should be standard
  # params_std <- ParametrisationTranslator(params = params, parametrisation = parametrisation, target='standard')
  # params_trf <- ParametrisationTranslator(params = params_std, parametrisation = 'standard', target='transform', target_alpha = target_alpha)
  
  jacob_function <- function(z){
    #takes transformed inputs
    z <- as.vector(z)
    z_non_zero <- z[z>0.0]
    upper <- eva::dgpd(z_non_zero, loc = 0, scale = params_std[2], shape = params_std[1])
    inv_z <- TransformationMapInverse(z_non_zero, params_std = params_std, params_trf = params_trf, target_alpha = target_alpha)
    lower <- eva::dgpd(inv_z, loc = 0, scale = params_trf[2], shape = params_trf[1])
    
    final_val <- rep(1.0, length(z))
    final_val[z>0.0] <- pmin(pmax(upper/lower, 1e-9), 1e9)
    
    return(final_val)
  }
}
# 
# set.seed(42)
# gpd_sample_1 <- eva::rgpd(n = 10000, loc = 0.0, scale = (1+10)/3, shape = 1/3.0)
# plot(sort(gpd_sample_1), sort(TransformationJacobian(c(3.0, 1, 10), c(3.0, 1, 10))(gpd_sample_1)))
# 
# set.seed(42)
# gpd_sample <- eva::rgpd(n = 10000, loc = 0.0, scale = 5, shape = -0.2)
# tmp <- TransformationMap(TransformationMapInverse(gpd_sample, params = c(-0.2, 5, 3), parametrisation = 'standard'),
#                                                   , params = c(-0.2, 5, 3), parametrisation = 'standard')
# plot(sort(tmp), sort(gpd_sample))
# line(0:20, 0:20)
# eva::gpdFit(tmp, 0.0)
