#' Simulate data with MOGC structure
#'
#' @param seed
#' @param K
#' @param numSamplesPerK
#' @param U2
#' @param l_index
#' @param numFeatInLevel1Group
#' @param strongGroupSignal
#' @param weakGroupSignalRatio
#' @param strongGroupPerc
#' @param weakGroupPerc
#' @param noiseMean
#' @param sample_sigma
#' @param sigma_noise
#' @param percConfounder
#' @param rho1
#' @param rho2
#' @param random_mu
#'
#' @return
#' @export
#'
#' @examples
simuMOGC <- function(seed=15213, K=3, numSamplesPerK,
  U2, l_index, numFeatInLevel1Group=3,
  strongGroupSignal, weakGroupSignalRatio=0.5, strongGroupPerc, weakGroupPerc,
  noiseMean=0, sample_sigma, sigma_noise,
  percConfounder=0, rho1, rho2, random_mu=FALSE){

  set.seed(seed)

  ## prepare samples for the subtypes
  nall = numSamplesPerK
  cumnall = cumsum(nall)
  n = sum(nall)
  clusterIndex = NULL
  label = numeric(n)
  for(i in 1:length(nall)){
    if(i==1){
      clusterIndex[[i]] = 1:cumnall[i]
    } else {
      clusterIndex[[i]] = (cumnall[i-1]+1):cumnall[i]
    }
    ## label the samples
    label[clusterIndex[[i]]] = i
  }

  ## generate mu (level-2 group by K)
  m1 <- nrow(U2)
  P <- m1 * numFeatInLevel1Group
  m2 <- ncol(U2)
  num_level1Groups_per_l <- apply(U2, 2, sum)
  mu_m2ByK <- matrix(0, nrow=m2, ncol=K) # number of groups * number of clusters matrix
  mu_per_Level2group_list <- list()
  for(l in 1:m2){
    if(random_mu){
      index <- sample(seq(1,K), replace=FALSE) # randomly permute
    }else{
      index <- seq(1,K)
    }
    if(l_index[l] == 1){
      mu_m2ByK[l, ] <- c(-strongGroupSignal,0,strongGroupSignal)[index]
    }
    mu_per_Level2group_list[[l]] <- matrix(0, num_level1Groups_per_l[l], K)
    strong_gene_number <- floor(num_level1Groups_per_l[l]*strongGroupPerc)
    weak_gene_number <- floor(num_level1Groups_per_l[l]*weakGroupPerc)
    noise_gene_number <- num_level1Groups_per_l[l] - strong_gene_number - weak_gene_number
    for(k in 1:K){
      mu_per_Level2group_list[[l]][,k] <- c(
        rep(mu_m2ByK[l,k], strong_gene_number),
        rep(mu_m2ByK[l,k] * weakGroupSignalRatio, weak_gene_number),
        rep(0, noise_gene_number))
    }
  }

  ## mu (level-1 group by K),include strong/weak genes inside each group
  mu_m1ByK <- matrix(0, nrow=m1, ncol=K)
  for(g in 1:m1){
    lgroup_index <- which(U2[g, ]==1)
    index_inside_each_lgroup <- sapply(1:length(lgroup_index),
      function(x) sum(U2[1:g, lgroup_index[x]]))
    if(length(lgroup_index) == 1){
      mu_m1ByK[g, ] <- mu_per_Level2group_list[[lgroup_index]][index_inside_each_lgroup, ]
    }else{
      mu_m1ByK[g, ] <- apply(sapply(1:length(lgroup_index), function(x)
        mu_per_Level2group_list[[lgroup_index[x]]][index_inside_each_lgroup[x],]), 1, sum) # cbind as columns
    }
  }

  ## extend mu to gene level (p by k)
  level1_group_index <- rep(seq(1, m1), each=numFeatInLevel1Group)
  mu_pByK <- mu_m1ByK[level1_group_index, ]

  geneLabel <- apply(mu_pByK!=0, 1, max)

  ## sample data (P by n)
  data <- mu_pByK[, label] + matrix(rnorm(P*n), P, n)

  # gene predictive label
  result <- list(data=t(data),label=label,
     mu_m2ByK=mu_m2ByK, mu_m1ByK=mu_m1ByK, mu_pByK=mu_pByK, geneLabel=geneLabel)

  return(result)
}
