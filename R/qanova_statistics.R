teststat_int <- function(data, grn, ord, m, zalpha_star, t , sig_help, C, C_mat, k ,n_vec, p, nadat, ...){
  
  data1 <- lapply(1:k, function(x){data[,nadat[1]][ grn == x]})
  
  ord <- ord
  m <- m
  sig_help <- sig_help
  n_vec_inv <- 1/n_vec
  sqrt2piinv <- 1/sqrt(2*pi)

 
  # Determine all needed quantiles (including the ones for the variance/density estimation)
  quant <- lapply(1:k, function(x){
    matrix(quantile(data1[[x]], probs = c(p, ord[1, (1:m) + m*(x -1) ]/n_vec[x], ord[2, 1:m + m*(x -1) ]/n_vec[x] ), type = 1), ncol = 3, byrow = FALSE)
    # First col represent the quantiles of interest
    # Second col consists of the (upper) quantiles needed for the interval estimation
    # Third col consists of the (lower) quantiles for the interval estimation
  })

  #####Variances
  # interval estimator
  finv_mat <- lapply( 1:k, function(x){
    sig <-  (quant[[x]][,2] - quant[[x]][,3]) / 2/ ( 1/sqrt(n_vec[x]) + zalpha_star[(1:m) + m*(x -1)]) 
    out <- matrix( rep( sig/ sqrt(p*(1-p)), m ), ncol = m)
  })
 
  Sigma <- lapply(1:k, function(x){
    finv_mat[[x]] * sig_help * t( finv_mat[[x]]) # Factor n of Sigma and the factor n^{1/2} in the statistic cancel out each other(quadratic form!)
    # Factors n_i in Sigma and in the density estimation cancel out each other
  })
  

  
  Cq <- sapply(1:k, function(x){ 
    C[[x]] %*% matrix( quant[[x]][,1] , ncol=1) 
  })
  Cq <- matrix( rowSums(matrix(Cq, ncol = k)), ncol = 1)
  
  HSigH <- sapply(1:k, function(x){
    ( C[[x]] %*% Sigma[[x]] %*% t(C[[x]]))
  })
  HSigH <- matrix( rowSums(matrix(HSigH, ncol = k)), ncol = length(C[[1]][,1]))
  
  
  stat <- t(Cq) %*% MASS::ginv( HSigH) %*% Cq
  
  out <- stat
  
}

  
teststat_boot <- function(data, grn, ord, m,  t , sig_help, C, C_mat, k ,n_vec, p, nadat, ...){
  
  data1 <- lapply(1:k, function(x){data[,nadat[1]][ grn == x]})
  
  ord <- ord
  m <- m
  sig_help <- sig_help
  n_vec_inv <- 1/n_vec
  sqrt2piinv <- 1/sqrt(2*pi)
  
  quant <- lapply(1:k, function(x){
    matrix(quantile(data1[[x]], probs = c(p, ord[1, (1:m) + m*(x -1) ]/n_vec[x], ord[2, 1:m + m*(x -1) ]/n_vec[x] ), type = 1), ncol = 3, byrow = FALSE)
    # First col represent the quantiles of interest
    # Second col consists of the (upper) quantiles needed for the interval estimation
    # Third col consists of the (lower) quantiles for the interval estimation
  })
  
  
  
  # Bootstrap Estimation ala Efron, we copied the formula from RPT::boot.var
  finv_mat_boot <- lapply( 1:k, function(xk){
    n_boot <- n_vec[xk]
    Prob <- matrix(sapply(seq(1:n_boot), function(x){ pbinom(t[,xk], size = n_boot, prob = (x - 1)/n_boot)}) - sapply(seq(1:n_boot), function(x){ pbinom(t[,xk], size = n_boot, prob = x/n_boot)}), ncol = n_boot)
    sig_sq <- sapply(1:m, function(x){
      sum( Prob[x,] * (sort(data1[[xk]]) - quant[[xk]][x,1] )^2 )
    })
    matrix( rep( sqrt(sig_sq)/ sqrt(p*(1-p)), m ), ncol = m)
  })
  
  
  Sigma_boot <- lapply(1:k, function(x){
    finv_mat_boot[[x]] * sig_help * t( finv_mat_boot[[x]]) 
    # Factor n of Sigma and the factor n^{1/2} in the statistic cancel out each other(quadratic form!)
    # Factors n_i in Sigma and in the density estimation cancel out each other
  })
  
  ####Our statistic
  Cq <- sapply(1:k, function(x){ 
    C[[x]] %*% matrix( quant[[x]][,1] , ncol=1) 
  })
  Cq <- matrix( rowSums(matrix(Cq, ncol = k)), ncol = 1)


  HSigH_boot <- sapply(1:k, function(x){
    ( C[[x]] %*% Sigma_boot[[x]] %*% t(C[[x]]))
  })
  HSigH_boot <- matrix( rowSums(matrix(HSigH_boot, ncol = k)), ncol = length(C[[1]][,1]))
  
  
  stat_boot <- t(Cq) %*% MASS::ginv ( HSigH_boot) %*% Cq
  
  out <- stat_boot
  
}


teststat_kernel <- function(data, grn, ord, m,  t , sig_help, C, C_mat, k ,n_vec, p, nadat, ...){
  
  data1 <- lapply(1:k, function(x){data[,nadat[1]][ grn == x]})
  
  
  ord <- ord
  m <- m
  sig_help <- sig_help
  n_vec_inv <- 1/n_vec
  sqrt2piinv <- 1/sqrt(2*pi)
  
  
  # Determine all needed quantiles (including the ones for the variance/density estimation)
  quant <- lapply(1:k, function(x){
    matrix(quantile(data1[[x]], probs = c(p, ord[1, (1:m) + m*(x -1) ]/n_vec[x], ord[2, 1:m + m*(x -1) ]/n_vec[x] ), type = 1), ncol = 3, byrow = FALSE)
    # First col represent the quantiles of interest
    # Second col consists of the (upper) quantiles needed for the interval estimation
    # Third col consists of the (lower) quantiles for the interval estimation
  })
  
  
  #Gaussian kernel estimation with thumb role for the bandwith
  finv_mat_kern <- lapply( 1:k, function(x){
    h <- bw.nrd0(data1[[x]])
    out <- NULL
    for( i in 1:m){
      out <-c(out, sqrt2piinv/h*sum(exp( -(quant[[x]][i,1] - data1[[x]])^2/2/h^2 )) )
    }
    out <- matrix( rep( 1/out, m ), ncol = m)
  })
  Sigma_kern <- lapply(1:k, function(x){
    n_vec[x]*finv_mat_kern[[x]] * sig_help * t( finv_mat_kern[[x]]) # Factor n of Sigma and the factor n^{1/2} in the statistic cancel out each other(quadratic form!)
  })
  
  
 
  
  ####Our statistic
  Cq <- sapply(1:k, function(x){ 
    C[[x]] %*% matrix( quant[[x]][,1] , ncol=1) 
  })
  Cq <- matrix( rowSums(matrix(Cq, ncol = k)), ncol = 1)
 
  HSigH_kern <- sapply(1:k, function(x){
    ( C[[x]] %*% Sigma_kern[[x]] %*% t(C[[x]]))
  })
  HSigH_kern <- matrix( rowSums(matrix(HSigH_kern, ncol = k)), ncol = length(C[[1]][,1]))
 
  stat_kern <- t(Cq) %*% MASS::ginv ( HSigH_kern) %*% Cq
 
  out <- stat_kern
  
}

