library('zeallot')



GridFoundations <- function(n, vanishing_depth, x=1){
  # returns a triangular matrix
  
  n <- n # accomodate for warm-up
  
  n_elems_per_line <- vanishing_depth
  index_i <- unlist(lapply(1:max(2, n_elems_per_line), FUN = function(i){return(1:(n))}))
  index_j <- unlist(lapply(1:max(2,n_elems_per_line), FUN = function(i){return(i:(n+i-1))}))
  return(Matrix::sparseMatrix(i=index_i, j=index_j, use.last.ij = T, x = x))
}

# n <- 14
# vanishing_depth <- 4
# gf <- GridFoundations(n = 14, vanishing_depth = vanishing_depth)
# gf[1:14,]
# k <- 13
# 
# gf[max(c(1, k-vanishing_depth+1)):min(c(n,k)), k:min(c(n,vanishing_depth+k-1))]

TrawlSlicing <- function(n, vanishing_depth, trawl_parameter, type='exp'){
  # returns matrix with
  # [[S(1,1), S(2,1), S(3,1), ..., S(vanishing_depth,1)],
  #  [S(2,2), S(3,2), S(4,2), ..., S(vanishing_depth+1,2)],
  #  ...,
  #  [S(n,n), S(3,2), S(4,2), ..., S(vanishing_depth+n,n)],
  # ]
  
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type=type)
  one_split <- (-B2_func(param = trawl_parameter, h = 0:vanishing_depth)) %>% diff
  slices <- matrix(rep(one_split, n), ncol = vanishing_depth, byrow = T)
  area_A <- B2_func(trawl_parameter, 0.0)
  slices[2:nrow(slices),] <- slices[2:nrow(slices),] - rep(1, n-1) %o% c(one_split[2:length(one_split)], 0.0)
  
  return(slices/area_A) # divide by \mu^{leb}(A)
}

# tmp <- TrawlSlicing(n = 5, vanishing_depth = 10, trawl_parameter = 0.3)

GammaGrid <- function(alpha, beta, n, vanishing_depth, trawl_parameter, type='exp'){
  n_diags <- max(2, vanishing_depth)
  n_non_zero_elements <- n_diags*n
  
  gamma_shapes <- GenerateShapes(alpha = alpha,
                                 n = n,
                                 vanishing_depth = vanishing_depth,
                                 trawl_parameter = trawl_parameter,
                                 type = type)
  stopifnot(length(gamma_shapes) == n_non_zero_elements)
  gamma_sim_vals <- rgamma(n = n_non_zero_elements,
              shape = gamma_shapes,
              rate = beta)
  return(GridFoundations(n = n, vanishing_depth = vanishing_depth, x = gamma_sim_vals))
}

GenerateShapes <- function(alpha, n, vanishing_depth, trawl_parameter, type='exp'){
  stopifnot(length(alpha) == 1)
  
  slices <- TrawlSlicing(n = n, vanishing_depth = vanishing_depth, trawl_parameter = trawl_parameter, type = type) # already standardised
  return(alpha*as.vector(slices))
}

BlockIndex <- function(n_block, n, vanishing_depth){
  return(
    list(
      block_i=max(c(1, n_block-vanishing_depth+1)):min(c(n,n_block)),
      block_j=n_block:min(c(n,vanishing_depth+n_block-1)))
    )
}

GammaOrchestra <- function(scaled_gamma_grid, parallel=T){
  n <- dim(scaled_gamma_grid)[1]
  vanishing_depth <- dim(scaled_gamma_grid)[2] - n + 1
 
  if(parallel){
    cores <- detectCores(logical = TRUE)
    cl <- makeCluster(cores)
    clusterExport(cl, c('BlockIndex', 'n', 'vanishing_depth'))
    clusterEvalQ(cl, library(Matrix))
    tmp <- parallel::parLapply(cl = cl,
                       X = n:1,
                       fun = function(i){
                         blck_ind <- BlockIndex(i, n = n, vanishing_depth = vanishing_depth)
                         return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
                       })
    parallel::stopCluster(cl)
    return(unlist(tmp))
  }else{
    return(vapply(1:n, FUN = function(i){
      blck_ind <- BlockIndex(i, n = n, vanishing_depth = vanishing_depth)
      return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
    }, FUN.VALUE = 1.0))
  }
}

n <- 10
vanishing_depth <- 3
gf <- GridFoundations(n, vanishing_depth)

tmp <- GammaGrid(alpha = 3, beta = 20, n = 2500, vanishing_depth = 30)
tmp
plot(GammaOrchestra(tmp))
# 
# c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type='exp')
# B2_func(0.3, 0:10)
# 
# (-B2_func(param = 0.8, h = 0:k)) %>% diff

GetVanishingCoverage <- function(trawl_parameter, vanishing_depth, type='exp', get_value=F){
  c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type=type)
  max_val <- B2_func(trawl_parameter, h = 0.0)

  coverage <- abs((-B2_func(param = trawl_parameter, h = c(0,vanishing_depth))) %>% diff %>% sum)/abs(max_val)
  if(get_value){
    return(coverage)
  }else{
    cat('Coverage:', round(coverage*100.0, 2), '%\n')
  }
}

# GetVanishingCoverage(trawl_parameter = 0.3, vanishing_depth = 10)

TrawlSimulation <- function(alpha, beta, n, vanishing_depth, trawl_parameter, type, parallel=F){
  gamma_grid <- GammaGrid(alpha = alpha,
                          beta = beta,
                          n = n,
                          vanishing_depth = vanishing_depth,
                          trawl_parameter = trawl_parameter,
                          type = type)
  return(GammaOrchestra(gamma_grid, parallel = parallel))
}

profvis::profvis({
  set.seed(42)
  trawl_sim <- TrawlSimulation(3,3,n=5000,vanishing_depth=30,0.3,'exp',parallel = F)
})
profvis::profvis({
  set.seed(42)
  trawl_sim <- TrawlSimulation(3,3,n=5000,vanishing_depth=30,0.3,'exp',parallel = T)
})

trawl_sim <- TrawlSimulation(3,3,n=10000,vanishing_depth=30,0.3,'exp')
trawl_sim %>% plot
trawl_sim %>% acf
lines(1:15, exp(-0.3*1:15))

plot(density(trawl_sim))
lines(0:50/10, dgamma(0:50/10, shape = 3, rate = 3))

ExceedancesSimulation <- function(params, parametrisation='standard', n, vanishing_depth, type, parallel=F){
  #TODO have a function that does exactly that: given a parametrisation type, changes it to noven
  
  params_trawl <- ParametrisationTranslator(params = params, parametrisation = parametrisation)
  trawl_parameter <- params[4:length(params)]
  alpha <- params_trawl[1]
  beta <- params_trawl[2]
  kappa <- params_trawl[3]
  probability_zero_model <- 1-(1+kappa/beta)^{-alpha}
  cat('probability_zero_model', probability_zero_model, '\n')
  trawl_simulation <- TrawlSimulation(alpha = alpha,
                                      beta = beta,
                                      trawl_parameter = trawl_parameter,
                                      n = n,
                                      vanishing_depth = vanishing_depth,
                                      type = type,
                                      parallel = parallel)
  trawl_simulation_unif <- TrawlSimulation(alpha = alpha,
                                      beta = beta,
                                      trawl_parameter = trawl_parameter,
                                      n = n,
                                      vanishing_depth = vanishing_depth,
                                      type = type,
                                      parallel = parallel)
  GetVanishingCoverage(trawl_parameter = trawl_parameter, vanishing_depth = vanishing_depth, type=type)
  plot(trawl_simulation)
  
  corr_uniform <- pgamma(trawl_simulation_unif, shape=alpha, rate=beta)
  acf(corr_uniform)
  probabilities_zero <- 1-exp(-kappa*trawl_simulation)
  print(summary(probabilities_zero))
  uniform_samples <- runif(n = n, min = 0, max = 1.0)
  uniform_samples <- corr_uniform
  exceedances <- apply(cbind(probabilities_zero, uniform_samples), MARGIN = 1,
                       FUN = function(p_and_u){
                        if(p_and_u[1] >= p_and_u[2]){
                          return(0.0)
                        }else{
                          return(eva::rgpd(1, loc = 0.0, scale = (beta+kappa)/alpha, shape = 1/alpha))
                        }
  })
  
  return(exceedances)
}

example_params_noven <- c(6.33, 20.12, 12.18, 0.1)
o <- ExceedancesSimulation(params = example_params_noven, parametrisation='noven', n = 5000, vanishing_depth = 30, type = 'exp')
# TODO Check trawl last one generated
(o>0) %>% mean
acf(o)
acf(as.numeric(o>0))
plot(o)
is_pos <- as.numeric(o>0)
flags_chunk_2 <- zoo::rollapply(is_pos, width=2, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0
abline(v=which(flags_chunk_2>0))
flags_chunk_3 <- zoo::rollapply(is_pos, width=3, FUN=function(chunk){sum(chunk)*(chunk[1])*prod(chunk)}) > 0
abline(v=which(flags_chunk_3>0), col='red')
acf(o)
acf(is_pos)
