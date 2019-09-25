TransformationMap <- function(x, params, target_alpha=3){
  # from original to trf
  target_beta <- 1+params[3]
  
  
}

TransformationMapInverse <- function(x, params, parametrisation){
  # fromt trf to original
  params[]
}

TransformationJacobian <- function(params, parametrisation, target_alpha=3){
  params_std <- ParametrisationTranslator(params = params, parametrisation = parametrisation, target='standard')
  
  jacob_function <- function(z){
    #takes transformed inputs
    
    upper <- eva::dgpd(z, loc = 0, scale = (1+params[3])*abs(target_alpha), shape = 1.0/target_alpha)
    
    inv_z <- TransformationMapInverse(z,)
    lower <- eva::dgpd(inv_z, loc = 0, scale = (1+params[3])*abs(target_alpha), shape = 1.0/target_alpha)
    return(upper/lower)
  }
  
}

