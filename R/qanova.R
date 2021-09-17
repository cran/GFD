#' QANOVA: Quantile-based analyis-of-variance
#'
#' The function \code{qanova} calculates the Wald-type statistic based on the
#' quantiles and/or their linear combinations, e.g. the interquartile range. Respective p-values
#' are obtained by a \eqn{\chi^2}-approximation and a permutation approach, respectively.
#' @param formula A model \code{formula} object. The left hand side contains the response variable and the right
#'  hand side contains the factor variables of interest. An interaction term must be
#'  specified.
#' @param data A data.frame, list or environment containing the variables in formula
#' and the censoring status
#' indicator. Default option is \code{NULL}.
#' @param quantiles A vector of probabilties corresponding to the quantiles of interest.
#' By default is c(0.5), i.e. just the median is included.
#' @param lin_mat A matrix specifying which linear combination of the quantiles should be
#' included for the analysis. By default (\code{NULL}) the identity matrix is chosen, i.e.
#' all chosen quantiles are considered simulatenously.
#' @param var_method Method for the variance estimation of the sample quantiles. The default
#' ("interval") is the interval-based estimator of Price and Bonett (2001). Additionally, the bootstrap
#' method ("boot") of Efron (1987) or a kernel density approach ("kernel") can be chosen.
#' @param nperm The number of permutations used for calculating the permuted p-value.
#'   The default option is 1999.
#' @param var_level A number between 0 and 1 specifying the confidence level for the
#'  interval variance estimation method; the default value is 0.95.
#' @param nested.levels.unique A logical specifying whether the levels of the nested
#' factor(s) are labeled uniquely or not.
#'  Default is FALSE, i.e., the levels of the nested factor are the same for each
#'  level of the main factor.
#' @details
#' The \code{qanova} function calculates the Wald-type statistic based on quantiles
#' and linear combinations of them for general factorial designs. The procedure is fully
#' nonparametric and no specific assumption of the underlying distribution is required.
#' In particular, heteroscedastic settings can be studied. The analysis can be based on a single
#' quantile (e.g. the median, default choice), a linear combination of quantiles (e.g. the
#' interquartile range, set \code{quantiles=c(0.25,0.75)} and \code{lin_mat = matrix(c(-1,1),ncol=2)})
#' or on several (combinations of) quantiles simulatenously.
#'
#'   The \code{qanova} function returns the test statistic as well as two
#'   corresponding p-values: the first is based on a \eqn{chi^2} approximation and
#'   the second one is based on a permutation procedure.
#'
#'  @return A \code{qanova} object containing the following components:
#'  \itemize{
#'  \item{pvalues_stat}{The p-values obtained by \eqn{\chi^2}-approximation}
#'  \item{pvalues_per}{The p-values of the permutation approach}
#'  \item{statistics}{The value of the qanova along with degrees of freedom of the
#'  central chi-square distribution and p-value, as well as the p-value of the
#'   permutation procedure.}
#'  \item{nperm}{The number of permutations used for calculating the permuted p-value.}
#'}
#'
#' @examples
#' QANOVA(weightgain ~ source*type, data = HSAUR::weightgain,var_method = "interval", nperm =199)
#' 
#' @references Ditzhaus, M., Fried, R. and Pauly, M. (2021). QANOVA: Quantile-based Permutation Methods For General
#' Factorial Designs. TEST (to appear, ArXiv preprint arXiv:1912.09146).
#' Efron, B. (1979). Bootstrap methods: Another look at the jackknife. Ann. Statist., 7:1-26.
#' Price, R. and Bonett, D. (2001). Estimating the variance of the sample median. J. Stat. Comput.
#' Simul, 68:295-305.
#' 
#' @author Philipp Steinhauer
#'
#' @importFrom  MASS ginv
#' @import plyr
#' 
#' @export


QANOVA <- function(formula, data = NULL, quantiles = c(0.5),lin_mat = NULL,
                   var_method = "interval", nperm = 1999, var_level = 0.95,  nested.levels.unique = FALSE){

  if(is.null(lin_mat) == TRUE){
    lin_mat = diag(rep(1,length(quantiles)))
  }

  if(ncol(lin_mat) != length(quantiles)){
    stop("The length of quantiles and the number of columns of lin_mat need to be the same.")
  }

  ## Check arguments

  if( !(all( (quantiles < 1)) & all( (quantiles > 0)) ) ){
    stop("The entries of quantiles need to be values between 0 and 1.")
  }




  input_list <- list(formula = formula, data = data, nperm = nperm,
                     var_method = var_method,lin_mat = lin_mat, p = quantiles)
  dat <- model.frame(formula, data)
  subject <- 1:nrow(dat)
  dat2 <- data.frame(dat, subject = subject)
  nf <- ncol(dat) - 1
  nadat <- names(dat)
  nadat2 <- nadat[-1]
  fl <- NA

  for (aa in 1:nf) {
    fl[aa] <- nlevels(as.factor(dat[, aa + 1]))
  }
  levels <- list()
  for (jj in 1:nf) {
    levels[[jj]] <- levels(as.factor(dat[, jj + 1]))
  }
  lev_names <- expand.grid(levels)
  if (nf == 1) {
    dat2 <- dat2[order(dat2[, 2]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure
    hypo_matrices <- list(diag(fl) - matrix(1/fl, ncol = fl, nrow = fl))
    
    k <- length(n)
    m <- length(quantiles)
    n_vec <- n
    
    gr <- rep(1:k, each = m)
    grn <- rep(1:k, n_vec)
    
    C_mat <- function(x){
      t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
    }
    
    
    hypo_matrices <- lapply(hypo_matrices, FUN = function(X) X %x% lin_mat)
    
    rank_C <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))
    hypo_matrices <-  lapply(hypo_matrices,C_mat)
    
    
    C_split <- function(Y){
      List <- list()
      for(i in 1:k){
        List[[i]] <- matrix(Y[ ,gr == i], ncol = m, byrow = FALSE)
      }
      List
    }
    
    C_List <-  lapply(hypo_matrices,C_split)
    
    
    mat1 <- matrix( c(rep(quantiles, each = m), rep(quantiles, m)), byrow = FALSE, ncol = 2)
    sig_help <- matrix( apply( mat1, 1, function(x){ min(x)- x[1]*x[2]}), ncol = m)
    
    alpha0 <- 1-var_level
    zalpha <- qnorm(1-alpha0/2)
    mat_gr_qu <- matrix( c(rep(1:k, each = m), rep(1:m, k )), ncol =2, byrow = FALSE) # First col represent the group, second col the quantile
    ord <- apply(mat_gr_qu, 1, function(x){
      c( min(n_vec[x[1]], floor( n_vec[x[1]]*quantiles[x[2]] + zalpha * sqrt( n_vec[x[1]]*quantiles[x[2]]*(1-quantiles[x[2]]) ))),
         max(1,floor( n_vec[x[1]]*quantiles[x[2]] - zalpha * sqrt( n_vec[x[1]]*quantiles[x[2]]*(1-quantiles[x[2]]) ))))})
    t <- matrix(sapply(1:k, function(x){ ceiling(n_vec[x] * quantiles )  - 1 }), ncol = k)
    
    
    if (var_method == "interval"){
      
      alpha_star <- sapply(1:length(ord[1,]), function(x){
        if( n_vec[mat_gr_qu[x,1]] > 1000){
          out <- alpha0
        }else{
          int <-  ( ord[2,x] + 1) : (ord[1,x] - 1)
          out <- 1 - sum( choose(n_vec[mat_gr_qu[x,1]], int) * quantiles[mat_gr_qu[x,2]]^(int) * (1-quantiles[mat_gr_qu[x,2]])^( n_vec[mat_gr_qu[x,1]]  -int ))
        }
      })
      zalpha_star <- qnorm(1-alpha_star/2)
      
      
      stat <-  sapply(1:length(hypo_matrices), function(x){
        teststat_int(data = dat2, grn = grn, ord = ord, m = m,
                     zalpha_star = zalpha_star, t = t, sig_help = sig_help,
                     C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                     nadat = nadat)
      })
      
      stat_perm <- sapply(1:length(hypo_matrices), function(x){
        replicate(nperm, teststat_int(data = dat2, grn = sample(grn), ord = ord, m = m,
                                      zalpha_star = zalpha_star, t = t, sig_help = sig_help,
                                      C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                      nadat = nadat))}
      )
      
      
      
    }
    
    
    if (var_method == "boot"){
      
      stat <-  sapply(1:length(hypo_matrices), function(x){
        teststat_boot(data = dat2, grn = grn, ord = ord, m = m, t = t, sig_help = sig_help,
                      C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                      nadat = nadat)
      })
      stat_perm <- sapply(1:length(hypo_matrices), function(x){
        replicate(nperm, teststat_boot(data = dat2, grn = sample(grn), ord = ord, m = m, t = t, sig_help = sig_help,
                                       C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                       nadat = nadat))}
      )
      
      
    }
    
    if (var_method == "kernel"){
      
      stat <-  sapply(1:length(hypo_matrices), function(x){
        teststat_kernel(data = dat2, grn = grn, ord = ord, m = m, t = t, sig_help = sig_help,
                        C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                        nadat = nadat)
        
      })
      stat_perm <- sapply(1:length(hypo_matrices), function(x){
        replicate(nperm, teststat_kernel(data = dat2, grn = sample(grn), ord = ord, m = m, t = t, sig_help = sig_help,
                                         C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                         nadat = nadat))}
      )
      
    }
    pvalues_stat <- sapply(1:length(stat), function(x) (1- pchisq(stat[x], df = rank_C[x])))
    pvalues_per <- sapply(1:length(stat), function(x) ( mean(stat[x] < stat_perm[,x], na.rm = TRUE )))
    
    output <- list()
    output$input <- input_list
    output$pvalues_stat  <- pvalues_stat
    output$pvalues_per <-  pvalues_per
    output$nperm <-nperm
    
    output$statistics <- cbind(stat,rank_C,round(pvalues_stat,4),round(pvalues_per,4))
    rownames(output$statistics) <- fac_names
    colnames(output$statistics) <- c("Test statistic","df","p-value", "p-value perm")
    



  }
  else {
    lev_names <- lev_names[do.call(order, lev_names[, 1:nf]), ]
    dat2 <- dat2[do.call(order, dat2[, 2:(nf + 1)]), ]
    response <- dat2[, 1]
    nr_hypo <- attr(terms(formula), "factors")
    fac_names <- colnames(nr_hypo)
    fac_names_original <- fac_names
    perm_names <- t(attr(terms(formula), "factors")[-1, ])
    ###

    n <- plyr::ddply(dat2, nadat2, plyr::summarise, Measure = length(subject),
                     .drop = F)$Measure

    if (length(fac_names) != nf && 2 %in% nr_hypo) {
      stop("A model involving both nested and crossed factors is\n           not implemented!")
    }
    if (length(fac_names) == nf && nf >= 4) {
      stop("Four- and higher way nested designs are\n           not implemented!")
    }
    if (length(fac_names) == nf) {
      TYPE <- "nested"
      if (nested.levels.unique) {
        n <- n[n != 0]
        blev <- list()
        lev_names <- list()
        for (ii in 1:length(levels[[1]])) {
          blev[[ii]] <- levels(as.factor(dat[, 3][dat[,
                                                      2] == levels[[1]][ii]]))
          lev_names[[ii]] <- rep(levels[[1]][ii], length(blev[[ii]]))
        }
        if (nf == 2) {
          lev_names <- as.factor(unlist(lev_names))
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev)
        }
        else {
          lev_names <- lapply(lev_names, rep, length(levels[[3]])/length(levels[[2]]))
          lev_names <- lapply(lev_names, sort)
          lev_names <- as.factor(unlist(lev_names))
          blev <- lapply(blev, rep, length(levels[[3]])/length(levels[[2]]))
          blev <- lapply(blev, sort)
          blev <- as.factor(unlist(blev))
          lev_names <- cbind.data.frame(lev_names, blev,
                                        as.factor(levels[[3]]))
        }
        if (nf == 2) {
          fl[2] <- fl[2]/fl[1]
        }
        else if (nf == 3) {
          fl[3] <- fl[3]/fl[2]
          fl[2] <- fl[2]/fl[1]
        }
      }
      hypo_matrices <- HN(fl)
    }
    else {
      TYPE <- "crossed"
      hypo_matrices <- HC(fl, perm_names, fac_names)[[1]]
      fac_names <- HC(fl, perm_names, fac_names)[[2]]
    }
    if (length(fac_names) != length(hypo_matrices)) {
      stop("Something is wrong: Perhaps a missing interaction term in formula?")
    }
    if (TYPE == "nested" & 0 %in% n & nested.levels.unique ==
        FALSE) {
      stop("The levels of the nested factor are probably labeled uniquely,\n           but nested.levels.unique is not set to TRUE.")
    }
    if (sum(n<5) > 0) {
      stop("There is at least one factor-level combination\n           with less than 5 observations!")
    }
    k <- length(n)
    m <- length(quantiles)
    n_vec <- n

    gr <- rep(1:k, each = m)
    grn <- rep(1:k, n_vec)

    C_mat <- function(x){
      t(x) %*% MASS::ginv( x %*% t(x) ) %*% x
    }


   hypo_matrices <- lapply(hypo_matrices, FUN = function(X) X %x% lin_mat)

   rank_C <- unlist(lapply(hypo_matrices, function(x) qr(x)$rank))
    hypo_matrices <-  lapply(hypo_matrices,C_mat)



    C_split <- function(Y){
      List <- list()
      for(i in 1:k){
       List[[i]] <- matrix(Y[ ,gr == i], ncol = m, byrow = FALSE)
      }
      List
    }

    C_List <-  lapply(hypo_matrices,C_split)




    mat1 <- matrix( c(rep(quantiles, each = m), rep(quantiles, m)), byrow = FALSE, ncol = 2)
    sig_help <- matrix( apply( mat1, 1, function(x){ min(x)- x[1]*x[2]}), ncol = m)

    alpha0 <- 1-var_level
    zalpha <- qnorm(1-alpha0/2)
    mat_gr_qu <- matrix( c(rep(1:k, each = m), rep(1:m, k )), ncol =2, byrow = FALSE) # First col represent the group, second col the quantile
    ord <- apply(mat_gr_qu, 1, function(x){
      c( min(n_vec[x[1]], floor( n_vec[x[1]]*quantiles[x[2]] + zalpha * sqrt( n_vec[x[1]]*quantiles[x[2]]*(1-quantiles[x[2]]) ))),
         max(1,floor( n_vec[x[1]]*quantiles[x[2]] - zalpha * sqrt( n_vec[x[1]]*quantiles[x[2]]*(1-quantiles[x[2]]) ))))})
    t <- matrix(sapply(1:k, function(x){ ceiling(n_vec[x] * quantiles )  - 1 }), ncol = k)



    if (var_method == "interval"){

      alpha_star <- sapply(1:length(ord[1,]), function(x){
        if( n_vec[mat_gr_qu[x,1]] > 1000){
          out <- alpha0
        }else{
          int <-  ( ord[2,x] + 1) : (ord[1,x] - 1)
          out <- 1 - sum( choose(n_vec[mat_gr_qu[x,1]], int) * quantiles[mat_gr_qu[x,2]]^(int) * (1-quantiles[mat_gr_qu[x,2]])^( n_vec[mat_gr_qu[x,1]]  -int ))
        }
      })
      zalpha_star <- qnorm(1-alpha_star/2)


      stat <-  sapply(1:length(hypo_matrices), function(x){
      teststat_int(data = dat2, grn = grn, ord = ord, m = m,
                   zalpha_star = zalpha_star, t = t, sig_help = sig_help,
                   C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                   nadat = nadat)
      })

       stat_perm <- sapply(1:length(hypo_matrices), function(x){
         replicate(nperm, teststat_int(data = dat2, grn = sample(grn), ord = ord, m = m,
                                             zalpha_star = zalpha_star, t = t, sig_help = sig_help,
                                             C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                             nadat = nadat))}
         )



    }


    if (var_method == "boot"){

      stat <-  sapply(1:length(hypo_matrices), function(x){
        teststat_boot(data = dat2, grn = grn, ord = ord, m = m, t = t, sig_help = sig_help,
                     C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                     nadat = nadat)
      })
      stat_perm <- sapply(1:length(hypo_matrices), function(x){
        replicate(nperm, teststat_boot(data = dat2, grn = sample(grn), ord = ord, m = m, t = t, sig_help = sig_help,
                                       C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                       nadat = nadat))}
      )
    

    }

    if (var_method == "kernel"){

      stat <-  sapply(1:length(hypo_matrices), function(x){
        teststat_kernel(data = dat2, grn = grn, ord = ord, m = m, t = t, sig_help = sig_help,
                      C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                      nadat = nadat)
        
      })
      stat_perm <- sapply(1:length(hypo_matrices), function(x){
        replicate(nperm, teststat_kernel(data = dat2, grn = sample(grn), ord = ord, m = m, t = t, sig_help = sig_help,
                                       C = C_List[[x]], C_mat = hypo_matrices[[x]],k = k, n_vec = n_vec, p = quantiles,
                                       nadat = nadat))}
      )
    
    }
    
    pvalues_stat <- sapply(1:length(stat), function(x) (1- pchisq(stat[x], df = rank_C[x])))
    pvalues_per <- sapply(1:length(stat), function(x) ( mean(stat[x] < stat_perm[,x], na.rm = TRUE )))
    
    output <- list()
    output$input <- input_list
    output$pvalues_stat  <- pvalues_stat
    output$pvalues_per <-  pvalues_per
    output$nperm <-nperm
    
    output$statistics <- cbind(stat,rank_C,round(pvalues_stat,4),round(pvalues_per,4))
    rownames(output$statistics) <- fac_names
    colnames(output$statistics) <- c("Test statistic","df","p-value", "p-value perm")

  }

  class(output) <- "qanova"
  return(output)
  


}




