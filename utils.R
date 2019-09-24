GetKappa <- function(data, params, parametrisation='standard'){
  # params in standard parametrisation
  p_non_zero <- mean(as.numeric(data>0))
  print(p_non_zero)
  
  if(parametrisation == 'standard'){
    return(params[2]/abs(params[1])*(1-p_non_zero^{params[1]}))
  }else{
    if(parametrisation == 'noven'){
      return(params[2]*(p_non_zero^{-1/params[1]}-1.0))
    }
  }
}

ParametrisationTranslator <- function(params, parametrisation, target='noven'){
  # from parametrisation to noven
  params_target <- params
  
  if(parametrisation == target){
    return(params)
  }
  if(parametrisation == 'standard' & target=='noven'){
    params_target[1] <- 1/params[1]
    params_target[2] <- params[2]/abs(params[1]) - params[3]
  }else{
    if(parametrisation='noven' & target=='standard'){
      params_target[1] <- 1/params[1]
      params_target[2] <- (params[2] + params[3])/abs(params[1])
    }
  }
  
  return(params_target)
}

