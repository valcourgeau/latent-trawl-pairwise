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

ParametrisationTranslator <- function(params, parametrisation){
  # from parametrisation to noven
  params_noven <- params
  if(parametrisation == 'standard'){
    params_noven[1] <- 1/params[1]
    params_noven[2] <- params[2]/abs(params[1]) - params[3]
  }
  
  return(params_noven)
}
