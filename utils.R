GetKappa <- function(data, params, parametrisation='standard'){
  # params in standard parametrisation
  p_non_zero <- mean(as.numeric(data>0))
  print(p_non_zero)
  
  if(parametrisation == 'standard'){
    if(params[1] < 0){
      return((p_non_zero^{-1./3.} - 1))
    }else{
      return(params[2]/abs(params[1])*(1-p_non_zero^{params[1]}))
    }
  }else{
    if(parametrisation == 'noven'){
      if(params[1] < 0){
        return((p_non_zero^{-1./3.} - 1))
      }else{
        return(params[2]*(p_non_zero^{-1/params[1]}-1.0))
      }
    }
  }
}

ParametrisationTranslator <- function(params, parametrisation, target='noven', target_alpha=3.0){
  # from parametrisation to noven
  params_target <- params
  
  if(parametrisation == target){
    return(params)
  }
  
  if(parametrisation == 'standard' & target == 'noven'){
    params_target[1] <- 1/params[1]
    params_target[2] <- params[2]/abs(params[1]) - params[3]
  }else{
    if(parametrisation == 'noven' & target == 'standard'){
      print('o')
      params_target[1] <- 1/params[1]
      params_target[2] <- (params[2] + params[3])/abs(params[1])
      print(params)
      print(params_target)
      print('o')
    }else{
      if(parametrisation == 'standard' & target == 'transform'){
        params_target[1] <- 1.0/target_alpha
        params_target[2] <- (1.0+params[3])*abs(params[1])
      }else{
        if(parametrisation == 'noven' & target == 'transform'){
          params_target[1] <- 1.0/target_alpha
          params_target[2] <- (1.0+params[3])/abs(params[1])
        }
      }
    }
  }
  
  return(params_target)
}

