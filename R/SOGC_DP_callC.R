#' Bayesian Gaussian Dirichlet Process mixture model incorporating Single-layer Overlapping Group structure (SOGC_dp)
#'
#' This function applies SOGC_dp clustering model to a feature matrix
#'
#' @param Y n by p feature matrix, where n is the sample size and p is feature size
#' @param U1 p by m1 matrix denoting feature group membership
#' @param K number of clusters
#' @param center TRUE/FALSE, if each column of Y should be centered
#' @param burnInIter number of iterations as burn-in period
#' @param keepIter number of iterations to be stored
#' @param maxIter maximum number of iterations
#' @param print_int print progress message every print_int iterations
#' @param debug TRUE/FALSE, indicating debug mode
#' @param init_mu K by p matrix, indicate the initial mu
#' @param init_z vector of length n, indicating cluster label for initial
#' @param fix_z TRUE/FALSE, if TRUE, do not update cluster label (for debugging)
#' @param fix_mu TRUE/FALSE, if TRUE, do not update cluster mean (for debugging)
#' @param adj_ls TRUE/FALSE, if TRUE, adjust for label switching using package label.switching
#' @param seed integer as random seed
#' @param BernoulliWeighted TRUE/FALSE, if true, adding R_j/D_g for fairness (see paper)
#' @param MH_ind TRUE/FALSE, if true, MH is used to update pi_gk. If overlapping group exists, MH_ind must be specified as TRUE
#' @param pi_j_prop_n integer, specifying n_0 in the proposal distribution of pi_gk for MH updating
#' @param s2 numeric number, s2 will be fixed and not be updated
#' @param types a vector of length p, indicating the types of each feature. Type-specific s2 will be estimated.
#' @param report_dup_mean TRUE/FALSE, if true, the v_jgk will be returned
#'
#' @return A list of MCMC output
#' @export
#'
#' @examples
SOGC_DP <- function(Y, U1, K, center=FALSE,
  burnInIter=1000, keepIter=3000, maxIter=10000,
  print_int=1000, debug=FALSE,
  init_mu=NULL, init_z=NULL, fix_z=FALSE, fix_mu=FALSE,
  adj_ls=TRUE, seed=123, BernoulliWeighted=TRUE, MH_ind=0,
  pi_j_prop_n=10, 
  s2=NULL, types=NULL, report_dup_mean=FALSE){

  set.seed(seed)
  library(Rcpp)
  library(RcppEigen)
  library(label.switching)
  library(coda)
  #library(MASS)
  #library(MCMCpack)

  ## Error handeling
  if(maxIter < keepIter){
    stop("maxIter has to be greater than keepIter.")
  }

  ## get some parameters from data
  N <- nrow(Y)
  P_dup <- P <- ncol(Y)
  if(center==TRUE){
    Y <- scale(Y, center=TRUE, scale=FALSE)
  }else{
    Y <- Y
  }
  G <- ncol(U1)
  pg <- apply(U1, 2, sum)
  g_index <- unlist(lapply(1:P, function(j) which(U1[j,]==1)))
  U1_rowsum <- rowSums(U1)
  if(sum(U1_rowsum==0)>1){
    stop("Feature must belongs to at least one group, make a singleton group for feature that does not belong to any groups")
  }

  # index to denote if feature is shared by more than one groups (in the same order as rows of U1)
  feature_dup_index <- seq(1, nrow(U1))
  nondup_features_id <- which(U1_rowsum == 1)
  dup_features_id <- which(U1_rowsum > 1)
  dup_features_id_num_dup <- U1_rowsum[dup_features_id]
  g_index <- unlist(lapply(1:nrow(U1), function(j) which(U1[j,]==1)))
  overlap_index <- !(all(U1_rowsum==1))
  if(overlap_index){
    feature_dup_index <- unlist(lapply(1:nrow(U1), function(j) rep(j, U1_rowsum[j])))
    P_dup <- length(feature_dup_index)
  }
  # if group overlap MH is used to update pi_g
  if(overlap_index){
    if(MH_ind==0){
      stop("When overlapping groups exist and Bernoulli needs to be corrected, MH_ind have to be set 1. \n")
    }
  }

  # weight for Bernoulli prior
  if(BernoulliWeighted){
    weight_for_Rj <- unlist(lapply(1:length(U1_rowsum),
      function(j) rep(1/U1_rowsum[j], U1_rowsum[j])))
  }else{
    weight_for_Rj <- rep(1, ncol(Y))
  }

  ## Random Initials
  eta1 <- eta2 <- 2
  alpha_dp <- 1
  sigma2 <- rep(1, P)
  pi_g <- rep(0.5, K)
  pi_j <- matrix(0.5, K, G)
  gamma_g <- matrix(rbinom(K*G, 1, 0.5), K, G)
  gamma_g_rep <- gamma_g[,g_index]
  gamma_j <- matrix(rbinom(K*P_dup, 1, 0.5), K, P_dup)
  gamma_j[, pg[g_index]==1] <- 1
  b <- matrix(rnorm(K*P_dup), K, P_dup)
  mu_dup <- gamma_g_rep*gamma_j*b

  # allow different feature type has different s2
  if(is.null(types)){
    types <- rep(1, ncol(Y))
  }
  uni_types <- unique(types)
  if(!is.null(s2)){
    s2_fixed <- TRUE
    s2 <- rep(s2,length(uni_types))
  }else{
    s2_fixed <- FALSE
    s2 <- rep(1, length(uni_types))
  }

  ## Initials of DP
  v <- c(rbeta(K-1, 1, alpha_dp), 1)
  pi_c <- v
  pi_c[1] <- v[1]
  pi_c[2:K] <- sapply(2:K, function(t) prod(1-v[1:(t-1)])*v[t])
  z <- sample(seq(1,K), size=N, replace=TRUE, prob=pi_c)

  get_mu <- function(mu_dup){
    mu <- matrix(NA, K, nrow(U1))
    mu[, nondup_features_id] <- mu_dup[, match(nondup_features_id, feature_dup_index)]
    if(length(dup_features_id) >0){
      for(k in 1:K){
        mu[k, dup_features_id] <- sapply(1:length(dup_features_id), function(x)
          sum(mu_dup[k, which(feature_dup_index==dup_features_id[x])]))
      }
    }
    return(mu)
  }
  get_mu_dup <- function(mu){
    mu_dup <- matrix(NA, K, P_dup)
    mu_dup[, match(nondup_features_id, feature_dup_index)] <- mu[, nondup_features_id]

    if(length(dup_features_id) >0){
      for(k in 1:K){
        for(i in 1:length(dup_features_id)){
          mu_dup[k, which(feature_dup_index==dup_features_id[i])] <-
            rep(mu[k, dup_features_id[i]]/dup_features_id_num_dup[i],
              dup_features_id_num_dup[i])
        }
      }
    }
    return(mu_dup)
  }
  mu <- get_mu(mu_dup)

  ## set informative initial if available
  if(!is.null(init_mu)){
    mu <- init_mu
    mu_dup <- get_mu_dup(mu)
    b[mu_dup!=0] <- mu_dup[mu_dup!=0]
    for(k in 1:K){
      gamma_g[k, ] <- sapply(1:G, function(g)sum(mu_dup[k,g_index==g]!=0) > 0)*1
    }
    gamma_j <- (mu_dup!=0)*1
    gamma_j[, pg[g_index]==1] <- 1
    print("Initialize mu by init_mu")
  }
  if(!is.null(init_z)){
    z <- init_z
    print("Initialize z by init_z")
    init_mu <- matrix(0, K, ncol(Y))
    uni_z <- unique(z)
    for(k_in in uni_z){
      init_mu[k_in, ] <- apply(Y[z==k_in,, drop=FALSE], 2, mean)
    }
    mu <- init_mu
    mu_dup <- get_mu_dup(mu)
    b[mu_dup!=0] <- mu_dup[mu_dup!=0]
    for(k in 1:K){
      gamma_g[k, ] <- sapply(1:G, function(g)sum(mu_dup[k,g_index==g]!=0) > 0)*1
    }
    gamma_j <- (mu_dup!=0)*1
    gamma_j[, pg[g_index]==1] <- 1
    print("then initialize mu")
  }
  nk <- sapply(1:K, function(x) sum(z==x))

  # set up necesary steps for MH of weighting in bernoulli prior due to overlapping groups
  BernoulliWeighted_int <- BernoulliWeighted*1
  get_pi_j_logl <- function(k, g){
    sel_index <- which(g_index==g)
    return(sum(gamma_j[k, sel_index]*log(weight_for_Rj[sel_index]*pi_j[k, g])) +
      sum((1-gamma_j[k, sel_index])*log(1-weight_for_Rj[sel_index]*pi_j[k, g])))
  }
  pi_j_loglikeli <- matrix(NA, K, G)
  for(k in 1:K){
    pi_j_loglikeli[k, ] <- sapply(1:G, function(g) get_pi_j_logl(k, g))
  }

  ## MCMC
  #sourceCpp(paste0(cppfile))
  res <- MCMC_sogc_dp(seed, burnInIter, keepIter, print_int, N, P, P_dup, G,
    K, BernoulliWeighted_int, MH_ind,
    s2, pi_j_prop_n,
    uni_types, types,
    eta1, eta2, alpha_dp, g_index,
    pg, z, nk, feature_dup_index, U1_rowsum,
    v, pi_c, sigma2, weight_for_Rj, gamma_g, gamma_j, pi_g, pi_j, b, mu, Y, debug,
    fix_z, fix_mu, pi_j_loglikeli, s2_fixed)

  # re-organize matrix into 3D array
  MU_store_o <- array(NA, dim=c(keepIter, K, P))
  Prob_Assign_store_o <- array(NA, dim=c(keepIter, N, K))
  for(k in 1:K){
    Prob_Assign_store_o[,,k] <- res$Prob_Assign_store[, (((k-1)*N+1):(k*N))]
    MU_store_o[, k, ] <- res$MU_store[,(((k-1)*P+1):(k*P))]
  }
  res$Prob_Assign_store <- Prob_Assign_store_o
  res$MU_store <- MU_store_o
  rm(Prob_Assign_store_o, MU_store_o)

  # check/correct for label switching
  if(adj_ls){
    reswitch <- label.switching::stephens(res$Prob_Assign_store)
    relabels <- reswitch$permutations #iterations*K matrix
    ## relabel cluster paremeters if labels are switched
    res$label_switched_index <- !all(sapply(1:K, function(x) length(unique(relabels[,x]))==1))
    if(res$label_switched_index==TRUE){
      res$MU_store <- label.switching::permute.mcmc(res$MU_store, relabels)$output
    }
  }

  # check convergence
  res$sample_index <- sample(seq(1, dim(res$MU_store)[3]), 5)
  mu_sub <- coda::as.mcmc(res$MU_store[, 1, res$sample_index])
  res$geweke_z <- coda::geweke.diag(mu_sub, frac1=0.1, frac2=0.5)$z
  res$convergence <- all(abs(res$geweke_z) < 1.96, na.rm=TRUE)
  totalIter <- keepIter

  # repeat if not converge
  repeat{
    if(res$convergence | totalIter >= maxIter){
      break
    }
    cat("MCMC not converge. Running more iterations ... \n")

    res <- MCMC_sogc_dp(seed, burnInIter=0, keepIter, print_int,
      N, P, P_dup, G,
      K, BernoulliWeighted_int, MH_ind,
      s2=res$S2_store[, keepIter], pi_j_prop_n,
      uni_types, types,
      eta1, eta2, alpha_dp, g_index, pg,
      z=res$Z_store[keepIter,],
      nk=res$NK_store[keepIter,], feature_dup_index, U1_rowsum,
      v=res$v, pi_c=res$PI_C_store[keepIter,],
      sigma2=res$SIGMA2_store[keepIter,],
      weight_for_Rj,
      gamma_g=matrix(res$GAMMA1_store[keepIter,], K, G, byrow=TRUE),
      gamma_j=matrix(res$GAMMA2_store[keepIter,], K, P_dup, byrow=TRUE),
      pi_g=res$pi_g_store[keepIter,],
      pi_j=matrix(res$pi_j_store[keepIter,], K, G, byrow=TRUE),
      b=matrix(res$B_store[keepIter,], K, P_dup, byrow=TRUE),
      mu=res$MU_store[keepIter,,],
      Y, debug, fix_z, fix_mu,
      pi_j_loglikeli=res$pi_j_loglikeli, s2_fixed)

    # re-organize matrix into 3D array (not perform necessary steps for convergence checkup)
    MU_store_o <- array(NA, dim=c(keepIter, K, P))
    Prob_Assign_store_o <- array(NA, dim=c(keepIter, N, K))
    for(k in 1:K){
      Prob_Assign_store_o[,,k] <- res$Prob_Assign_store[, (((k-1)*N+1):(k*N))]
      MU_store_o[, k, ] <- res$MU_store[,(((k-1)*P+1):(k*P))]
    }
    res$Prob_Assign_store <- Prob_Assign_store_o
    res$MU_store <- MU_store_o
    rm(Prob_Assign_store_o, MU_store_o)

    # check/correct for label switching (not perform necessary steps for convergence checkup)
    if(adj_ls){
      reswitch <- label.switching::stephens(res$Prob_Assign_store)
      relabels <- reswitch$permutations #iterations*K matrix
      ## relabel cluster paremeters if labels are switched
      res$label_switched_index <- !all(sapply(1:K, function(x) length(unique(relabels[,x]))==1))
      if(res$label_switched_index==TRUE){
        res$MU_store <- label.switching::permute.mcmc(res$MU_store, relabels)$output
      }
    }

    # check convergence
    res$sample_index <- sample(seq(1, dim(res$MU_store)[3]), 5)
    mu_sub <- coda::as.mcmc(res$MU_store[, 1, res$sample_index])
    res$geweke_z <- coda::geweke.diag(mu_sub, frac1=0.1, frac2=0.5)$z
    res$convergence <- all(abs(res$geweke_z) < 1.96, na.rm=TRUE)
    totalIter <- totalIter + keepIter
  }
  res$totalIter <- totalIter

  ## re-organize other matrices into 3-D arrays
  if(report_dup_mean){
    B_store_o <- GAMMA2_store_o <- array(NA, dim=c(keepIter, K, P_dup))
    pi_j_store_o <- GAMMA1_store_o <- array(NA, dim=c(keepIter, K, G))
    for(k in 1:K){
      GAMMA1_store_o[, k, ] <- res$GAMMA1_store[,(((k-1)*G+1):(k*G))]
      GAMMA2_store_o[, k, ] <- res$GAMMA2_store[,(((k-1)*P_dup+1):(k*P_dup))]
      B_store_o[, k, ] <- res$B_store[,(((k-1)*P_dup+1):(k*P_dup))]
      pi_j_store_o[, k, ] <- res$pi_j_store[,(((k-1)*G+1):(k*G))]
    }
    res$GAMMA1_store <- GAMMA1_store_o
    res$GAMMA2_store <- GAMMA2_store_o
    res$B_store <- B_store_o
    res$pi_j_store <- pi_j_store_o
    rm(GAMMA1_store_o, GAMMA2_store_o, B_store_o, pi_j_store_o)
  }

  if(adj_ls & res$label_switched_index==TRUE){
    res$PI_C_store <- t(sapply(1:nrow(res$PI_C_store), function(x)
      res$PI_C_store[x, match(relabels[x,], seq(1, K))]))
    res$Z_store <- t(sapply(1:nrow(res$Z_store), function(x)
      relabels[x, res$Z_store[x, ]]))
    if(report_dup_mean){
      res$GAMMA1_store <- label.switching::permute.mcmc(res$GAMMA1_store, relabels)$output
      res$GAMMA2_store <- label.switching::permute.mcmc(res$GAMMA2_store, relabels)$output
      res$B_store <- label.switching::permute.mcmc(res$B_store, relabels)$output
      res$pi_j_store <- label.switching::permute.mcmc(res$pi_j_store, relabels)$output
    }
  }

  res$g_index <- g_index

  ###### feature selection probability ######
  MU_nonzero <- res$MU_store!=0
  Overall_selection_matrix <- matrix(NA, dim(MU_nonzero)[1], dim(MU_nonzero)[3])
  for(i in 1:nrow(Overall_selection_matrix)){
    temp_z <- unique(res$Z_store[i,])
    if(length(temp_z)==1){
      Overall_selection_matrix[i,] <- MU_nonzero[i, temp_z,]
    }else{
      Overall_selection_matrix[i,] <- apply(MU_nonzero[i,temp_z,],2,max)
    }
  }
  res$feature_select_prob <- apply(Overall_selection_matrix, 2, mean)

  ##### dup feature selection probability ####
  if(report_dup_mean){
    res$MU_dup_store <- res$GAMMA1_store[,,g_index]*res$GAMMA2_store*res$B_store
    MU_nonzero <- res$MU_dup_store!=0
    Overall_selection_matrix <- matrix(NA, dim(MU_nonzero)[1], dim(MU_nonzero)[3])
    for(i in 1:nrow(Overall_selection_matrix)){
      temp_z <- unique(res$Z_store[i,])
      if(length(temp_z)==1){
        Overall_selection_matrix[i,] <- MU_nonzero[i, temp_z,]
      }else{
        Overall_selection_matrix[i,] <- apply(MU_nonzero[i,temp_z,],2,max)
      }
    }
    res$feature_dup_select_prob <- apply(Overall_selection_matrix, 2, mean)
  }

  # MAP cluster labels
  get_map_label <- function(x) {
    tb <- table(x)
    return(names(tb)[which.max(tb)])
  }
  label_map <- sapply(1:ncol(res$Z_store), function(x)
    get_map_label(res$Z_store[,x]))
  res$label_map <- as.numeric(label_map)
  res$feature_dup_index <- feature_dup_index

  return(res)
}















