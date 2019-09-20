library('zeallot')

GridFoundations <- function(n, vanishing_depth, x=1){
  # returns a triangular matrix
  i <- rep(1:n, vanishing_depth)
  index_i <- unlist(lapply(1:max(2,(2*vanishing_depth-1)), FUN = function(i){return(1:(n-i+1))}))
  index_j <- unlist(lapply(1:max(2,(2*vanishing_depth-1)), FUN = function(i){return(i:(n))}))
  return(Matrix::sparseMatrix(i=index_i, j=index_j, use.last.ij = T, x = x))
}
n <- 14
vanishing_depth <- 4
gf <- GridFoundations(n = 14, vanishing_depth = 5)
k <- 13

gf[max(c(1, k-vanishing_depth+1)):min(c(n,k)), k:min(c(n,vanishing_depth+k-1))]

GammaGrid <- function(alpha, beta, n, vanishing_depth){
  n_diags <- max(2, 2*vanishing_depth-1)
  n_non_zero_elements <- n_diags*(n + n-n_diags+1)/2
  x <- rgamma(n = n_non_zero_elements, shape = alpha, rate = beta)
  return(GridFoundations(n = n, vanishing_depth = vanishing_depth, x = x))
}

GenerateShapes <- function(n, trawl_parameters, type='exp'){
  trawl_funcs %<-% GetTrawlFunctions(type=type)
  B1_func <- trawl_funcs[1]
  B2_func <- trawl_funcs[2]
  
}

BlockIndex <- function(n_block, n, vanishing_depth){
  return(
    list(
      block_i=max(c(1, n_block-vanishing_depth+1)):min(c(n,n_block)),
      block_j=n_block:min(c(n,vanishing_depth+n_block-1)))
    )
}

GammaOrchestra <- function(scaled_gamma_grid){
  n <- dim(scaled_gamma_grid)[1]
  vanishing_depth <- dim(scaled_gamma_grid)[2]
  return(vapply(1:n, FUN = function(i){
    blck_ind <- BlockIndex(i, n = n, vanishing_depth = vanishing_depth)
    # print(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j])
    return(sum(scaled_gamma_grid[blck_ind$block_i, blck_ind$block_j]))
    }, FUN.VALUE = 1.0))
}

n <- 10
vanishing_depth <- 3
gf <- GridFoundations(n, vanishing_depth)
gf
summary(gf)
gf[1,2]

tmp <- GammaGrid(alpha = 3, beta = 20, n = 2500, vanishing_depth = 30)
tmp
plot(GammaOrchestra(tmp))

c(B1_func, B2_func, B3_func) %<-% GetTrawlFunctions(type='exp')
B2_func(0.3, 1:10)
