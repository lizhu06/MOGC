#' Bayesian Gaussian Dirichlet Process mixture model incorporating Multi-layer Overlapping Group structure (MOGC_dp)
#'
#' This function applies MOGC clustering model to a feature matrix
#'
#' @param Y n by p feature matrix, where n is the sample size and p is feature size
#' @param U1 p by m1 matrix denoting feature group membership
#' @param U2 m1 by m2 matrix denoting level-2 group membership matrix of level-1 groups
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
#' @param pi_g_prop_n integer, specifying n_0 in the proposal distribution of pi_lk for MH updating
#' @param pi_j_prop_n integer, specifying n_0 in the proposal distribution of pi_glk for MH updating
#' @param s2 numeric number, s2 will be fixed and not be updated
#' @param types a vector of length p, indicating the types of each feature. Type-specific s2 will be estimated.
#' @param report_dup_mean TRUE/FALSE, if true, the v_jgk will be returned
#'
#' @return A list of MCMC output
#' @export
#'
#' @examples
MOGC_DP <- function(Y, U1, U2, K, center=FALSE,
  burnInIter=1000, keepIter=3000, maxIter=10000,
  print_int=1000, debug=FALSE,
  init_mu=NULL, init_z=NULL, fix_z=FALSE, fix_mu=FALSE,
  adj_ls=TRUE, seed=123,
  BernoulliWeighted=TRUE, MH_ind=0, pi_g_prop_n=10, pi_j_prop_n=10,
  s2=NULL, types=NULL,
  report_dup_mean=FALSE){

  set.seed(seed)
  library(Rcpp)
  library(RcppEigen)
  library(label.switching)
  library(coda)

  ## Error handeling
  if(maxIter < keepIter){
    stop("maxIter has to be greater than keepIter.")
  }

  ## get some parameters from data
  N <- nrow(Y)
  P <- ncol(Y)
  if(center==TRUE){
    Y <- scale(Y, center=TRUE, scale=FALSE)
  }else{
    Y <- Y
  }
  U1_rowsum <- rowSums(U1)
  U2_rowsum <- rowSums(U2)
  if(sum(U1_rowsum==0)>1){
    stop("Feature must belongs to at least one group, make a singleton group for feature that does not belong sto any groups")
  }
  if(sum(U2_rowsum==0)>1){
    stop("Level-1 group must belongs to at least one level-2 group, make a singleton group for level-1 group that does not belong sto any level-2 groups")
  }
  overlap_index <- !(all(U1_rowsum==1) & all(U2_rowsum==1))
  if(overlap_index){
    cat("Data include overlapping groups. \n")
  }
  if(overlap_index){
    if(MH_ind==0){
      stop("When overlapping groups exist and Bernoulli needs to be corrected, MH_ind have to be set 1. \n")
    }
  }

  # dup features in the order of level-1 groups
  feature_dup_index_U1 <- g_index_U1 <- NULL
  for(j in 1:nrow(U1)){
    feature_dup_index_U1 <- c(feature_dup_index_U1,
      rep(j, U1_rowsum[j]))
    g_index_U1 <- c(g_index_U1, which(U1[j,]==1))
  }

  # duplicate level-1 group (in the order of level-1 groups)
  feature_dup_index <- g_index <- l_index <- l_index_g <- NULL
  for(g in 1:nrow(U2)){
    g_counter <- length(unique(g_index))
    feature_dup_index <- c(feature_dup_index,
      rep(feature_dup_index_U1[which(g_index_U1==g)],
      times=U2_rowsum[g]))
    g_index <- c(g_index, rep(c((g_counter+1):(g_counter+U2_rowsum[g])),
      each=sum(g_index_U1==g)))
    l_index <- c(l_index, rep(which(U2[g, ]==1), each=sum(g_index_U1==g)))
    l_index_g <- c(l_index_g, which(U2[g, ]==1))
  }
  P_dup <- length(feature_dup_index)

  Rj <- sapply(1:length(feature_dup_index), function(x)
    sum(feature_dup_index==feature_dup_index[x]))

  m1 <- length(unique(g_index))
  m2 <- length(unique(l_index))
  pg <- sapply(1:m1, function(k) sum(g_index==g))
  pl <- sapply(1:m2, function(l) sum(l_index==l))
  gl <- sapply(1:m2, function(l) sum(l_index_g==l))

  # weight for Bernoulli prior
  if(BernoulliWeighted){
    weight_for_Dg <- (1/U2_rowsum)
    weight_for_Dg <- rep(weight_for_Dg, times=U2_rowsum)
    weight_for_Tj <- (1/U1_rowsum)[feature_dup_index]
  }else{
    weight_for_Dg <- rep(1, length(U2_rowsum))
    weight_for_Dg <- rep(weight_for_Dg, times=U2_rowsum)
    weight_for_Tj <- rep(1, length(Rj))
  }

  ## Random Initials
  eta1 <- eta2 <- 2
  alpha_dp <- 1
  sigma2 <- rep(1, P)
  pi_k <- rep(1/K, K)   # mixture probability
  alpha0 <- rep(1, K)
  pi_l <- rep(0.5, K)   # pi for gamma_l
  pi_g <- matrix(0.5, K, m2) # pi for gamma_g
  pi_j <- matrix(0.5, K, m1) # pi for gamma_j

  # indicator of L
  gamma_l <- matrix(rbinom(K*m2, 1, 0.5), K, m2)
  gamma_l_rep <- matrix(NA, K, P_dup)
  for(k in 1:K){
    #gamma_l_rep[k, ] <- rep(gamma_l[k,], times=pl)
    gamma_l_rep[k, ] <- gamma_l[k,l_index]
  }
  # indicator of G
  gamma_g <- matrix(rbinom(K*m1, 1, 0.5), K, m1)
  gamma_g_rep <- matrix(NA, K, P_dup)
  for(k in 1:K){
    #gamma_g_rep[k, ] <- rep(gamma_g[k,], times=pg)
    gamma_g_rep[k, ] <- gamma_g[k,g_index]
  }
  # indicator of j
  gamma_j <- matrix(rbinom(K*P_dup, 1, 0.5), K, P_dup)
  b <- matrix(rnorm(K*P_dup), K, P_dup)

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
  z <- sample(seq(1,K), size=N, replace=TRUE, prob=rep(1/K, K))
  nk <- sapply(1:K, function(x) sum(z==x))
  mu_dup <- gamma_l_rep*gamma_g_rep*gamma_j*b

  get_mu <- function(mu_dup){
    mu <- mu_dup
    dup_feat_index <- feature_dup_index[duplicated(feature_dup_index)]
    if(length(dup_feat_index)>0){
      remove_index <- NULL
      for(jj in 1:length(dup_feat_index)){
        keep_index <- which(feature_dup_index==dup_feat_index[jj])[1]
        remove_index <- c(remove_index,
          which(feature_dup_index==dup_feat_index[jj])[-1])
        mu[, keep_index] <- apply(mu_dup[,
            which(feature_dup_index==dup_feat_index[jj])], 1, sum)
      }
      mu <- mu[, -remove_index]
      mu <- mu[, match(seq(1,P), feature_dup_index[-remove_index])]
    }else{
      mu <- mu[, match(seq(1,P), feature_dup_index)]
    }
    return(mu)
  }

  get_mu_dup <- function(mu){
    mu_dup <- mu[, feature_dup_index]
    dup_feat_index <- feature_dup_index[duplicated(feature_dup_index)]
    if(length(dup_feat_index)>0){
      for(jj in 1:length(dup_feat_index)){
        sel_index <- which(feature_dup_index==dup_feat_index[jj])
        mu_dup[, sel_index] <- mu_dup[,sel_index]/length(sel_index)
      }
    }
    return(mu_dup)
  }
  mu <- get_mu(mu_dup)

  ## Initials of DP
  v <- c(rbeta(K-1, 1, alpha_dp), 1)
  pi_k <- v
  pi_k[1] <- v[1]
  pi_k[2:K] <- sapply(2:K, function(t) prod(1-v[1:(t-1)])*v[t])
  z <- sample(seq(1,K), size=N, replace=TRUE, prob=pi_k)

  ## set informative initial if available
  if(!is.null(init_mu)){
    mu <- init_mu
    mu_dup <- get_mu_dup(mu)
    b[mu_dup!=0] <- mu_dup[mu_dup!=0]
    for(k in 1:K){
      gamma_l[k, ] <- sapply(1:m2, function(l)sum(mu_dup[k,l_index==l]!=0) > 0)*1
      gamma_g[k, ] <- sapply(1:m1, function(g)sum(mu_dup[k,g_index==g]!=0) > 0)*1
    }
    gamma_j <- (mu_dup!=0)*1
    print("Initialize mu by init_mu")
  }
  if(!is.null(init_z)){
    print("Initialize z by init_z")
    z <- init_z
    if(is.null(init_mu)){
      init_mu <- t(sapply(1:K, function(k) apply(Y[z==k,], 2, mean)))
      init_mu[is.na(init_mu)] <- 0
      mu <- init_mu
      mu_dup <- get_mu_dup(mu)
      b[mu_dup!=0] <- mu_dup[mu_dup!=0]
      for(k in 1:K){
        gamma_l[k, ] <- sapply(1:m2, function(l)sum(mu_dup[k,l_index==l]!=0) > 0)*1
        gamma_g[k, ] <- sapply(1:m1, function(g)sum(mu_dup[k,g_index==g]!=0) > 0)*1
      }
      gamma_j <- (mu_dup!=0)*1
      print("then initialize mu")
    }
  }
  nk <- sapply(1:K, function(x) sum(z==x))

  # set up necesary steps for MH of weighting in bernoulli prior due to overlapping groups
  BernoulliWeighted_int <- BernoulliWeighted*1
  get_pi_g_logl <- function(k, l){
    sel_index <- which(l_index_g==l)
    single_logl <- sum(gamma_g[k, sel_index]*log(weight_for_Dg[sel_index]*pi_g[k, l])) +
      sum((1-gamma_g[k, sel_index])*log(1-weight_for_Dg[sel_index]*pi_g[k, l]))
    return(single_logl)
  }
  pi_g_loglikeli <- matrix(NA, K, m2)
  for(k in 1:K){
    pi_g_loglikeli[k,] <- sapply(1:m2, function(l) get_pi_g_logl(k, l))
  }
  get_pi_j_logl <- function(k, g){
    sel_index <- which(g_index==k)
    single_logl <- sum(gamma_j[k, sel_index]*log(weight_for_Tj[sel_index]*pi_j[k, g])) +
    sum((1-gamma_j[k, sel_index])*log(1-weight_for_Tj[sel_index]*pi_j[k, g]))
    return(single_logl)
  }
  pi_j_loglikeli <- matrix(NA, K, m1)
  for(k in 1:K){
    pi_j_loglikeli[k, ] <- sapply(1:m1, function(g) get_pi_j_logl(k, g))
  }

  ## MCMC
  #sourceCpp(paste0(cppfile))
  res <- MCMC_mogc_dp(seed, burnInIter, keepIter, print_int, N, P, P_dup, m1, m2, K,
    BernoulliWeighted_int, MH_ind,
    s2, pi_g_prop_n, pi_j_prop_n,
    uni_types, types,
    eta1, eta2, alpha_dp,
    g_index, l_index, l_index_g,
    pg, pl, gl, z, nk, feature_dup_index,
    v, pi_j, pi_g, pi_l, sigma2, pi_k, weight_for_Dg, weight_for_Tj,
    gamma_j, gamma_g, gamma_l, b, mu, Y, debug, fix_z, fix_mu,
    pi_g_loglikeli, pi_j_loglikeli, s2_fixed)

  ## re-organozie matrix into 3-D array
  MU_store_o <- array(NA, dim=c(keepIter, K, P))
  Prob_Assign_store_o <- array(NA, dim=c(keepIter, N, K))
  for(k in 1:K){
    Prob_Assign_store_o[,,k] <- res$Prob_Assign_store[, (((k-1)*N+1):(k*N))]
    MU_store_o[, k, ] <- res$MU_store[,(((k-1)*P+1):(k*P))]
  }
  res$Prob_Assign_store <- Prob_Assign_store_o
  res$MU_store <- MU_store_o
  rm(Prob_Assign_store_o, MU_store_o)

  # check/correct for label switching (not correct mu for convergence check)
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

    res <- MCMC_mogc_dp(seed, burnInIter=0, keepIter, print_int, N, P, P_dup,
      m1, m2, K,
      BernoulliWeighted_int, MH_ind,
      s2=res$S2_store[, keepIter], pi_g_prop_n, pi_j_prop_n,
      uni_types, types,
      eta1, eta2, alpha_dp,
      g_index, l_index, l_index_g,
      pg, pl, gl, z=res$Z_store[keepIter,],
      nk=res$NK_store[keepIter,], feature_dup_index,
      v, pi_j=matrix(res$PI_J_store[keepIter,], K, m1, byrow=TRUE),
      pi_g=matrix(res$PI_G_store[keepIter,], K, m2, byrow=TRUE),
      pi_l=res$PI_L_store[keepIter,], sigma2=res$SIGMA2_store[keepIter,],
      pi_k=res$PI_K_store[keepIter,], weight_for_Dg, weight_for_Tj,
      gamma_j=matrix(res$GAMMA_J_store[keepIter,], K, P_dup, byrow=TRUE),
      gamma_g=matrix(res$GAMMA_G_store[keepIter,], K, m1, byrow=TRUE),
      gamma_l=matrix(res$GAMMA_L_store[keepIter,], K, m2, byrow=TRUE),
      b=matrix(res$B_store[keepIter,], K, P_dup, byrow=TRUE),
      mu=res$MU_store[keepIter,,],
      Y, debug, fix_z, fix_mu,
      pi_g_loglikeli=res$pi_g_loglikeli,
      pi_j_loglikeli=res$pi_j_loglikeli, s2_fixed)

    ## re-organozie matrix into 3-D array
    MU_store_o <- array(NA, dim=c(keepIter, K, P))
    Prob_Assign_store_o <- array(NA, dim=c(keepIter, N, K))
    for(k in 1:K){
      Prob_Assign_store_o[,,k] <- res$Prob_Assign_store[, (((k-1)*N+1):(k*N))]
      MU_store_o[, k, ] <- res$MU_store[,(((k-1)*P+1):(k*P))]
    }
    res$Prob_Assign_store <- Prob_Assign_store_o
    res$MU_store <- MU_store_o
    rm(Prob_Assign_store_o, MU_store_o)

    # check/correct for label switching (not correct mu for convergence check)
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
    B_store_o <- GAMMA_J_store_o <- array(NA, dim=c(keepIter, K, P_dup))
    PI_J_store_o <- GAMMA_G_store_o <- array(NA, dim=c(keepIter, K, m1))
    PI_G_store_o <- GAMMA_L_store_o <- array(NA, dim=c(keepIter, K, m2))
    for(k in 1:K){
      GAMMA_L_store_o[, k, ] <- res$GAMMA_L_store[,(((k-1)*m2+1):(k*m2))]
      GAMMA_G_store_o[, k, ] <- res$GAMMA_G_store[,(((k-1)*m1+1):(k*m1))]
      GAMMA_J_store_o[, k, ] <- res$GAMMA_J_store[,(((k-1)*P_dup+1):(k*P_dup))]
      B_store_o[, k, ] <- res$B_store[,(((k-1)*P_dup+1):(k*P_dup))]
      PI_G_store_o[, k, ] <- res$PI_G_store[,(((k-1)*m2+1):(k*m2))]
      PI_J_store_o[, k, ] <- res$PI_J_store[,(((k-1)*m1+1):(k*m1))]
    }
    rm(GAMMA_L_store_o, GAMMA_G_store_o, GAMMA_J_store_o,
      PI_G_store_o, PI_J_store_o)
  }

  ## label switching using label.switching r-package, stephens method
  if(adj_ls & res$label_switched_index==TRUE){
    res$PI_K_store <- t(sapply(1:nrow(res$PI_K_store), function(x)
      res$PI_K_store[x, match(relabels[x,], seq(1, K))]))
    res$Z_store <- t(sapply(1:nrow(res$Z_store), function(x)
      relabels[x, res$Z_store[x, ]]))
    res$MU_store <- label.switching::permute.mcmc(res$MU_store, relabels)$output
    if(report_dup_mean){
      res$GAMMA_L_store <- label.switching::permute.mcmc(res$GAMMA_L_store, relabels)$output
      res$GAMMA_G_store <- label.switching::permute.mcmc(res$GAMMA_G_store, relabels)$output
      res$GAMMA_J_store <- label.switching::permute.mcmc(res$GAMMA_J_store, relabels)$output
      res$B_store <- label.switching::permute.mcmc(res$B_store, relabels)$output
      res$PI_G_store <- label.switching::permute.mcmc(res$PI_G_store, relabels)$output
      res$PI_J_store <- label.switching::permute.mcmc(res$PI_J_store, relabels)$output
    }
  }
  res$l_index <- l_index
  res$g_index <- g_index

  ###### Overall feature AUC ######
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
    res$MU_dup_store <- res$GAMMA_L_store[,,l_index] *
      res$GAMMA_G_store[,,g_index] *
      res$GAMMA_J_store * res$B_store
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











