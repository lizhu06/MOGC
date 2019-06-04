#' Simulate data with SOGC structure
#'
#' @param seed
#' @param K
#' @param numSamplesPerK
#' @param U1
#' @param g_index
#' @param strongSignal
#' @param weakSignalRatio
#' @param strongSignalPerc
#' @param weakSignalPerc
#' @param noiseMean
#' @param sample_sigma
#' @param sigma_noise
#' @param percConfounder
#' @param rho
#' @param random_mu
#'
#' @return
#' @export
#'
#' @examples
simuSOGC <- function(seed=15213, K=3, numSamplesPerK,
  U1, g_index, strongSignal, weakSignalRatio=0.5,
  strongSignalPerc, weakSignalPerc,
  noiseMean=0, sample_sigma, sigma_noise,
  percConfounder=0, rho=0, random_mu=FALSE){

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

  ## generate mu at group level (group by K)
  P <- nrow(U1)
  m1 <- ncol(U1)
  num_genes_per_group <- apply(U1, 2, sum)
  mu_m1ByK <- matrix(0, nrow=m1, ncol=K) # number of groups * number of clusters matrix
  mu_per_group_list <- list()
  for(g in 1:m1){

    if(random_mu){
      index <- sample(seq(1,K), replace=FALSE) # randomly permute
    }else{
      index <- seq(1,K)
    }
    if(g_index[g] == 1){
      mu_m1ByK[g, ] <- c(-strongSignal,0,strongSignal)[index]
    }
    mu_per_group_list[[g]] <- matrix(0, num_genes_per_group[g], K)
    strong_gene_number <- floor(num_genes_per_group[g]*strongSignalPerc)
    weak_gene_number <- floor(num_genes_per_group[g]*weakSignalPerc)
    noise_gene_number <- num_genes_per_group[g] - strong_gene_number - weak_gene_number
    for(k in 1:K){
      mu_per_group_list[[g]][,k] <- c(
        rep(mu_m1ByK[g,k], strong_gene_number),
        rep(mu_m1ByK[g,k] * weakSignalRatio, weak_gene_number),
        rep(0, noise_gene_number))
    }
  }

  ## extend mu to gene level (p by k), include strong/weak genes inside each group
  mu_pByK <- matrix(0, nrow=P, ncol=K)
  for(p in 1:P){
    group_index <- which(U1[p, ]==1)
    index_inside_each_group <- sapply(1:length(group_index),
      function(x) sum(U1[1:p,group_index[x]]))
    if(length(group_index) == 1){
      mu_pByK[p, ] <- mu_per_group_list[[group_index]][index_inside_each_group, ]
    }else{
      mu_pByK[p, ] <- apply(sapply(1:length(group_index), function(x)
        mu_per_group_list[[group_index[x]]][index_inside_each_group[x],]), 1, sum) # cbind as columns
    }
  }
  geneLabel <- apply(mu_pByK!=0, 1, max)

  ## sample data (P by n)
  data <- mu_pByK[, label] + matrix(rnorm(P*n), P, n)

  # gene predictive label
  result <- list(data=t(data),label=label,
    mu_m1ByK=mu_m1ByK, mu_pByK=mu_pByK, geneLabel=geneLabel)

  return(result)
}
