// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <R.h>
#include <Rmath.h>
#include <Rcpp.h>

using namespace Eigen;
using namespace Rcpp;
using namespace std;
using Eigen::Map;                 
using Eigen::MatrixXd;             
using Eigen::VectorXd; 
using Eigen::MatrixXi;                  
using Eigen::VectorXi;                   

// [[Rcpp::export]]
void set_seed_mogc_dp(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
void update_z_mogc_dp(int N, int P, int K, Eigen::VectorXi& z,
  Eigen::VectorXd& pi_k, Eigen::VectorXd& sigma2, 
  Eigen::MatrixXd& mu, Eigen::MatrixXd& Y, Eigen::MatrixXd& Prob_Assign, bool debug) {

  //MatrixXi Prob_Assign_store(N, K);
  for(int i = 0; i < N; i++){
    NumericVector logden_c(K);
    NumericVector un_norm_prob_c(K);
    NumericVector norm_prob_c(K);

    for(int k = 0; k < K; k++){
      double logden = 0.0;
      double resi = (-0.5)*log(2*3.1415);
      for(int p = 0; p < P; p++){
        if(!(mu.col(p).sum()==0.0)){
          //logden = logden + ( - pow((Y(i, p) - mu(k, p)), 2) / (2.0 * sigma2[p]));
          logden = logden + (resi - 0.5 * log(sigma2[p]) - 
            pow((Y(i, p) - mu(k, p)), 2) / (2.0 * sigma2[p]));
        }
      }
      logden_c[k] = logden;
    }

    /* select the largest, and set as the offset*/
    double logden_offset = logden_c[0];
    for(int k = 1; k < K; k++){
      if(logden_offset < logden_c[k]){
        logden_offset = logden_c[k];
      }
    }

    /* get un-normalized prob for each k */
    for(int k = 0; k < K; k++){
      un_norm_prob_c[k] = pi_k[k] * exp(logden_c[k] - logden_offset);
      //if(i==0 && debug){
      //  cout << "i=" << i << ", k=" << k << ", logden_c= " << logden_c[k] << ".\n";
      //  cout << "i=" << i << ", k=" << k << ", un_norm_prob_c= " << un_norm_prob_c[k] << ".\n";
      //}
    }

    /* normalzied the probability*/
    double sum_prob = 0.0;
    for(int k = 0; k < K; k++){
      sum_prob = sum_prob + un_norm_prob_c[k];
    }

    /* Add small numbers to avoid denominator being zero*/
    if(sum_prob == 0.0){
      sum_prob = sum_prob + 3*DBL_MIN;
      for(int k = 0; k < K; k++){
        norm_prob_c[k] = un_norm_prob_c[k] + DBL_MIN;
      }
    }

    for(int k = 0; k < K; k++){
      norm_prob_c[k] = un_norm_prob_c[k] / sum_prob;
      Prob_Assign(i, k) = norm_prob_c[k];
    }

    IntegerVector ans(K);
    rmultinom(1, norm_prob_c.begin(), K, ans.begin());
    for(int k = 0; k < K; k++){
      if(ans[k] == 1){
        z[i] = k + 1;
      }
    }

  }
  //return(Prob_Assign_store);
}

// [[Rcpp::export]]
Eigen::VectorXi get_nk_mogc_dp(int N, int K,  Eigen::VectorXi& z){
  Eigen::VectorXi nk = Eigen::VectorXi::Zero(K);
  for(int i = 0; i < N; i++){
    nk[z[i]-1] = nk[z[i]-1] + 1;
  }
  return(nk);
}

// [[Rcpp::export]]
void update_v_mogc_dp(int K, double alpha_dp, Eigen::VectorXi& nk, Eigen::VectorXd& v){
 for(int k=0; k<K-1; k++){
  double nk_tail=0;
  for(int k2=k+1; k2<K; k2++){
    nk_tail = nk_tail + nk[k2];
  }
  v[k] = R::rbeta(1.0+nk[k], alpha_dp+nk_tail);
 }
 v[K-1] = 1;
}

// [[Rcpp::export]]
void update_pi_k_mogc_dp(int K, Eigen::VectorXd& v, Eigen::VectorXd& pi_k){
  pi_k[0] = v[0];
  for(int k=0; k<K; k++){
    double pi_k_prod = 1-v[0];
    for(int k2=0; k2<k; k2++){
      pi_k_prod = pi_k_prod*(1-v[k2]);
    }
    pi_k[k] = pi_k_prod*v[k];
  }
}

// [[Rcpp::export]]
double update_alpha_dp_mogc_dp(int K, double alpha_dp, double eta1, double eta2, Eigen::VectorXd& v){
  double log_sum_v = 0.0;
  for(int k=0; k<K-1; k++){
    log_sum_v = log_sum_v + log(1-v[k]);
  }
  log_sum_v = log_sum_v + pow(10.0, -10.0);
  alpha_dp = R::rgamma(K + eta1 - 1, eta2 - log_sum_v);
  return(alpha_dp);
}

// [[Rcpp::export]]
Eigen::MatrixXd update_mu_mogc_dp(int P, int P_dup, int K, Eigen::VectorXi& g_index, Eigen::VectorXi& l_index,
  Eigen::VectorXi& feature_dup_index,
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l, 
  Eigen::MatrixXd& b){

  Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(K, P);
  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int p = feature_dup_index[pd] - 1;
      int g = g_index[pd] - 1;
      int l = l_index[pd] - 1;
      mu(k, p) = mu(k, p) + gamma_l(k, l)*
        gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
    }
  }
  return(mu);
}

// [[Rcpp::export]]
void update_gamma_l_mogc_dp(int P_dup, int N, int K, int m2, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& l_index, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index,
  Eigen::VectorXd& pi_l, Eigen::VectorXd& sigma2,
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l, 
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  /*cout<< "gamma_1 " << gamma_1.row(0) << ".\n";*/
  for(int k = 0; k < K; k++){
    for(int l = 0; l < m2; l++){

      int old_gamma_l = gamma_l(k,l);
      double inside = 0.0;
      for(int pd = 0; pd < P_dup; pd++){
        int p = feature_dup_index[pd] - 1;
        int g = g_index[pd] - 1;
        if(l_index[pd] == l + 1){
          double sumy = 0.0;
          for(int i = 0; i < N; i++){
            if(z[i] == k+1){
              sumy = sumy + Y(i, p) - mu(k, p) + gamma_l(k, l)*
                gamma_g(k, g)*gamma_j(k,pd)*b(k,pd);
            }
          }
          inside = inside + (nk[k] * pow(gamma_g(k,g) * 
            gamma_j(k, pd) * b(k, pd), 2) - 2.0 * gamma_g(k,g) * 
            gamma_j(k, pd) * b(k, pd) * sumy) / (2.0 * sigma2[p]);
        }
      }
      //cout<< "k=" << k << ", g=" << g << ", pi1(k)= " << pi1[k] << ".\n";
      //cout<< "inside= " << inside << ".\n";
      double p_gamma_l = 1.0 / (1.0 + (1.0 - pi_l[k])/(pi_l[k]) * 
        exp(inside));
      //cout<< "p_gamma_1 " << p_gamma1 << ".\n";
      gamma_l(k, l) = R::rbinom(1, p_gamma_l);

       // if gamma_l change, update mu
      if(gamma_l(k,l) != old_gamma_l){
       for(int pd = 0; pd < P_dup; pd++){
         int p = feature_dup_index[pd] - 1;
         int g = g_index[pd] - 1;
         int l_2 = l_index[pd] - 1;
         if(l_2==l){
           mu(k, p) = mu(k, p) - old_gamma_l *
             gamma_g(k, g)*gamma_j(k, pd)*b(k, pd) + 
             gamma_l(k, l)*gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
         }
       }
      }

      //cout<< "gamma_1 " << gamma_1(k, g) << ".\n";
    }
  }
  /*cout<< "gamma_1 " << gamma_1.row(0) << ".\n";*/
}
// [[Rcpp::export]]
void update_gamma_g_mogc_dp(int P_dup, int N, int K, int m1, int BernoulliWeighted_int,
  Eigen::VectorXi& g_index, Eigen::VectorXi& l_index, Eigen::VectorXi& l_index_g, 
  Eigen::VectorXi& nk, Eigen::VectorXi& z,  Eigen::VectorXi& feature_dup_index,
  Eigen::VectorXd& sigma2,
  Eigen::VectorXd& weight_for_Dg,
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l, 
  Eigen::MatrixXd& pi_g, Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  /*cout<< "gamma_1 " << gamma_1.row(0) << ".\n";*/
  for(int k = 0; k < K; k++){
    for(int g = 0; g < m1; g++){
      int old_gamma_g = gamma_g(k,g);
      double inside = 0.0;
      for(int pd = 0; pd < P_dup; pd++){
        int p = feature_dup_index[pd] - 1;
        int l = l_index[pd] - 1;
        if(g_index[pd] == g + 1){
          double sumy = 0.0;
          for(int i = 0; i < N; i++){
            if(z[i] == k+1){
              sumy = sumy + Y(i, p)- mu(k, p) + gamma_l(k, l)*
                gamma_g(k, g)*gamma_j(k,pd)*b(k,pd);
            }
          }
          inside = inside + (nk[k] * pow(gamma_l(k,l)*gamma_j(k, pd)*b(k, pd), 2) 
            - 2.0 * gamma_l(k,l)*gamma_j(k, pd) * b(k, pd) * sumy) / (2.0 * sigma2[p]);
        }
      }
      double pi_g_weighted = pi_g(k, l_index_g[g] - 1);
      if(BernoulliWeighted_int == 1){
        pi_g_weighted = pi_g(k,l_index_g[g]-1)*weight_for_Dg[g];
      }
      double p_gamma_g = 1.0 / (1.0 + (1.0 - pi_g_weighted)/pi_g_weighted * 
        exp(inside));
      gamma_g(k, g) = R::rbinom(1, p_gamma_g);

       // if gamma_g change, update mu
      if(gamma_g(k,g) != old_gamma_g){
       for(int pd = 0; pd < P_dup; pd++){
         int p = feature_dup_index[pd] - 1;
         int g_2 = g_index[pd] - 1;
         int l = l_index[pd] - 1;
         if(g_2==g){
           mu(k, p) = mu(k, p) - old_gamma_g *
             gamma_l(k, l)*gamma_j(k, pd)*b(k, pd) + 
             gamma_l(k, l)*gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
         }
       }
      }
    }
  }
}


// [[Rcpp::export]]
void update_gamma_j_mogc_dp(int P_dup, int N, int K, int BernoulliWeighted_int, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& l_index, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index,
  Eigen::VectorXd& sigma2, Eigen::VectorXd& weight_for_Tj,
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l, 
  Eigen::MatrixXd& pi_j, Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int old_gamma_j = gamma_j(k,pd);
      int p = feature_dup_index[pd] - 1;
      int g = g_index[pd] - 1;
      int l = l_index[pd] - 1;
      double sumy = 0.0;
        for(int i = 0; i < N; i++){
          if(z[i] == k + 1){
            sumy = sumy + Y(i, p)- mu(k, p) + gamma_l(k, l)*
                gamma_g(k, g)*gamma_j(k,pd)*b(k,pd);
          }
        }
        double inside = (nk[k] * pow(gamma_l(k, l)*gamma_g(k, g) 
          * b(k, pd), 2) - 2.0 * gamma_l(k, l)*gamma_g(k, g) * 
          b(k, pd) * sumy)/(2.0 * sigma2[p]);

        double pi_j_weighted = pi_j(k, g);
        if(BernoulliWeighted_int == 1){
          pi_j_weighted = pi_j(k, g)*weight_for_Tj[pd];
        }
        double p_gamma_j = 1.0/(1.0 + (1.0 - pi_j_weighted)/
          pi_j_weighted * exp(inside));
        gamma_j(k, pd) = R::rbinom(1, p_gamma_j);
      //}
       // if gamma_g change, update mu
      if(gamma_j(k,pd) != old_gamma_j){
        mu(k, p) = mu(k, p) - old_gamma_j *
          gamma_l(k, l)*gamma_g(k, g)*b(k, pd) + 
          gamma_l(k, l)*gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
      }
    }
  }
  //return(gamma_2);
}

// [[Rcpp::export]]
void update_b_mogc_dp(int P_dup, int N, int K,  
  Eigen::VectorXd& s2,
  Eigen::VectorXi& g_index, Eigen::VectorXi& l_index, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index, Eigen::VectorXi& types,
  Eigen::VectorXd& sigma2,
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l,
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){
 
  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int p = feature_dup_index[pd] - 1;
      int g = g_index[pd] - 1;
      int l = l_index[pd] - 1;

      int a_type = types[p] - 1;

      if(gamma_l(k, l) * gamma_g(k,g) * gamma_j(k, pd) == 0){
          b(k,pd) = R::rnorm(0.0, sqrt(s2[a_type]));
      }else{
        double var_b = sigma2[p] * s2[a_type] / (s2[a_type] * nk[k] + sigma2[p]);
        double sumy = 0.0;
        double mu_remove_this_b = mu(k, p) - b(k,pd);
        for(int i = 0; i < N; i++){
          if(z[i] == k + 1){
            sumy = sumy + Y(i, p)- mu(k, p) + gamma_l(k, l)*
                gamma_g(k, g)*gamma_j(k,pd)*b(k,pd);
          }
        }
        double mean_b = sumy * var_b / sigma2[p];
        b(k, pd) =  R::rnorm(mean_b, sqrt(var_b));
        mu(k, p) = mu_remove_this_b + b(k, pd);
      }
    }
  }
  //return(b);
}

// [[Rcpp::export]]
void update_sigma2_mogc_dp(int P, int N, Eigen::VectorXi& z, Eigen::VectorXd& sigma2, 
  Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  //VectorXd sigma2(P); 
  for(int p = 0; p < P; p++){
    double res = 0.0;
    for(int i = 0; i < N; i++){
      res = res + pow(Y(i, p) - mu(z[i]-1, p), 2);
    }
    sigma2[p] = 1.0 / R::rgamma(N/2.0, 2.0/res);
  }
  //return(sigma2);
}

// [[Rcpp::export]]
void update_pi_l_mogc_dp(int K, int m2, 
  Eigen::VectorXd& pi_l, Eigen::MatrixXi& gamma_l){

  for(int k = 0; k < K; k++){
    double sum_gamma_l = gamma_l.row(k).sum();
    pi_l[k] = R::rbeta(sum_gamma_l + 1.0, m2 - sum_gamma_l + 1.0);
  }
  //return(pi1);
}

// [[Rcpp::export]]
void update_pi_g_mogc_dp(int P, int m2, int m1, int K, int BernoulliWeighted_int,
  int MH_ind, double pi_g_prop_n, 
  Eigen::VectorXi& gl, Eigen::VectorXi& l_index_g, Eigen::VectorXd& weight_for_Dg,
  Eigen::MatrixXi& gamma_g,
  Eigen::MatrixXd& pi_g, Eigen::MatrixXd& pi_g_loglikeli){

  for(int k = 0; k < K; k++){
    for(int l = 0; l < m2; l++){
      double fl = l + 1;
      if(MH_ind==1){
        double loglikeli_new=0.0;
        double pi_g_new = R::rbeta(pi_g_prop_n * pi_g(k, l), 
          pi_g_prop_n * (1.0 - pi_g(k, l)));
        for(int g = 0; g < m1; g++){
          if(l_index_g[g] == fl){
            if(gamma_g(k,g)==1){
              loglikeli_new = loglikeli_new + 
                log(pi_g_new * weight_for_Dg[g]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(1-(pi_g_new * weight_for_Dg[g]));
            }
          }
        }
        double acc_p = min(0.0, loglikeli_new - pi_g_loglikeli(k,l) + 
          R::dbeta(pi_g(k,l), pi_g_prop_n*pi_g_new, pi_g_prop_n*(1-pi_g_new), 1) -
          R::dbeta(pi_g_new, pi_g_prop_n*pi_g(k,l), pi_g_prop_n*(1-pi_g(k,l)), 1));
        double u_uni = R::runif(0, 1);
        if(log(u_uni) < acc_p){
          pi_g(k,l) = pi_g_new;
          pi_g_loglikeli(k,l) = loglikeli_new;
        }
      }else{
        double sum_gamma_g = 0.0;
        for(int g = 0; g < m1; g++){
          if(l_index_g[g] == l + 1){
            sum_gamma_g = sum_gamma_g + gamma_g(k, g);
          }
        }
        pi_g(k, l) = R::rbeta(sum_gamma_g + 1.0, gl[l] - sum_gamma_g + 1.0);
      }
    }
  }
  //return(pi2);
}

// [[Rcpp::export]]
void update_pi_j_mogc_dp(int P_dup, int m1, int K, int BernoulliWeighted_int,
  int MH_ind, double pi_j_prop_n, 
  Eigen::VectorXi& pg, Eigen::VectorXi& g_index, Eigen::VectorXd& weight_for_Tj,
  Eigen::MatrixXi& gamma_j,
  Eigen::MatrixXd& pi_j, Eigen::MatrixXd& pi_j_loglikeli){

  for(int k = 0; k < K; k++){
    for(int g = 0; g < m1; g++){
      double fg = g + 1;
      if(MH_ind == 1){
        double loglikeli_new=0.0;
        double pi_j_new = R::rbeta(pi_j_prop_n * pi_j(k, g), 
          pi_j_prop_n * (1.0 - pi_j(k, g)));

        for(int pd = 0; pd < P_dup; pd++){
          if(g_index[pd] == fg){
            if(gamma_j(k,pd)==1){
              loglikeli_new = loglikeli_new + 
                log(pi_j_new * weight_for_Tj[pd]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(1-(pi_j_new * weight_for_Tj[pd]));
            }
          }
        }

        double acc_p = min(0.0, loglikeli_new - pi_j_loglikeli(k,g) + 
          R::dbeta(pi_j(k,g), pi_j_prop_n*pi_j_new, pi_j_prop_n*(1-pi_j_new), 1) -
          R::dbeta(pi_j_new, pi_j_prop_n*pi_j(k,g), pi_j_prop_n*(1-pi_j(k,g)), 1));
        double u_uni = R::runif(0, 1);
        if(log(u_uni) < acc_p){
          pi_j(k,g) = pi_j_new;
          pi_j_loglikeli(k,g) = loglikeli_new;
        }

      }else{
        double sum_gamma_j = 0.0;
        for(int pd = 0; pd < P_dup; pd++){
          if(g_index[pd] == g + 1){
            sum_gamma_j = sum_gamma_j + gamma_j(k, pd);
          }
        }
        pi_j(k, g) = R::rbeta(sum_gamma_j + 1.0, pg[g] - sum_gamma_j + 1.0);
      }
    }
  }
  //return(pi2);
}

// [[Rcpp::export]]
void update_s2_mogc_dp(int P_dup, int K, int T, Eigen::VectorXi& uni_types, 
  Eigen::VectorXi& types, Eigen::VectorXi& feature_dup_index, 
  Eigen::VectorXd& s2, Eigen::MatrixXd& b){

//  int T = uni_types.size(); 

  for(int t = 0; t < T; t++){
    double sum_t = 0.0;
    double sum_b2 = 0.0;
    for(int k = 0; k < K; k++){
      for(int pd = 0; pd < P_dup; pd++){
        int p = feature_dup_index[pd] - 1;
        if(types[p] == uni_types[t]){
          sum_t = sum_t + 1;
          sum_b2 = sum_b2 + pow(b(k, pd), 2);
        }
      }
    }
    double arg1 = sum_t / 2.0 ;
    double arg2 = sum_b2 / 2.0 ;
    s2[t] = 1.0 / R::rgamma(arg1, 1.0/arg2);
  }
}


// [[Rcpp::export]]
List MCMC_mogc_dp(int seed, int burnInIter, int keepIter, int print_int, int N, 
  int P, int P_dup, int m1, int m2, int K, 
  int BernoulliWeighted_int, int MH_ind, 
  Eigen::VectorXd& s2, double pi_g_prop_n, double pi_j_prop_n,
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types,
  double eta1, double eta2, double alpha_dp, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& l_index, Eigen::VectorXi& l_index_g,
  Eigen::VectorXi& pg, Eigen::VectorXi& pl, Eigen::VectorXi& gl,
  Eigen::VectorXi& z, Eigen::VectorXi& nk, 
  Eigen::VectorXi& feature_dup_index,
  Eigen::VectorXd& v, 
  Eigen::MatrixXd& pi_j, Eigen::MatrixXd& pi_g, Eigen::VectorXd& pi_l,
  Eigen::VectorXd& sigma2, 
  Eigen::VectorXd& pi_k, Eigen::VectorXd& weight_for_Dg, Eigen::VectorXd& weight_for_Tj, 
  Eigen::MatrixXi& gamma_j, Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_l,
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, 
  Eigen::MatrixXd& Y, bool debug, bool fix_z, bool fix_mu,
  Eigen::MatrixXd& pi_g_loglikeli, Eigen::MatrixXd& pi_j_loglikeli,
  bool s2_fixed){

  int NSIM = burnInIter + keepIter;

  /* store matrix */
  Eigen::MatrixXi Z_store(keepIter, N);
  int nrow = K * P;
  Eigen::MatrixXd MU_store(keepIter, nrow);
  Eigen::MatrixXd PI_K_store(keepIter, K);
  Eigen::MatrixXd Prob_Assign(N, K);
  nrow = N * K;
  Eigen::MatrixXd Prob_Assign_store(keepIter, nrow);

  int num_t = uni_types.size();

  Eigen::MatrixXd S2_store(num_t, keepIter);
  Eigen::VectorXd ALPHA_DP_store(keepIter);
  Eigen::MatrixXi NK_store(keepIter, K);
  Eigen::MatrixXi GAMMA_J_store(keepIter, P_dup*K);
  Eigen::MatrixXi GAMMA_G_store(keepIter, m1*K);
  Eigen::MatrixXi GAMMA_L_store(keepIter, m2*K);
  Eigen::MatrixXd B_store(keepIter, P_dup*K);
  Eigen::MatrixXd SIGMA2_store(keepIter, P);
  Eigen::MatrixXd PI_L_store(keepIter, K);
  Eigen::MatrixXd PI_G_store(keepIter, K*m2);
  Eigen::MatrixXd PI_J_store(keepIter, K*m1);

  set_seed_mogc_dp(seed);
  
  /* MCMC */
  for(int nsim = 0; nsim < NSIM; nsim++){

    if(!fix_mu){  
      update_gamma_l_mogc_dp(P_dup, N, K, m2, g_index, l_index, nk, z,
      feature_dup_index,  pi_l, 
        sigma2, gamma_j, gamma_g, gamma_l, b, mu, Y);
      
      update_gamma_g_mogc_dp(P_dup, N, K, m1, BernoulliWeighted_int,
        g_index, l_index, l_index_g, 
        nk, z, feature_dup_index, 
        sigma2, weight_for_Dg, gamma_j, gamma_g, gamma_l, 
        pi_g, b, mu, Y);
      
      update_gamma_j_mogc_dp(P_dup, N, K, BernoulliWeighted_int,
        g_index, l_index, nk, z,  feature_dup_index, 
        sigma2, weight_for_Tj, gamma_j, gamma_g, gamma_l, pi_j, b, mu, Y);

      /* update b */
      update_b_mogc_dp(P_dup, N, K, s2, g_index, l_index,
        nk, z,  feature_dup_index, types,  sigma2, gamma_j, gamma_g,gamma_l, b, mu, Y);

      /* update mu*/
      //mu = update_mu(P, P_dup, K, g_index, l_index, feature_dup_index, gamma_j, gamma_g, gamma_l, 
      //  b);
    }

    /* update sigma2 */
    update_sigma2_mogc_dp(P, N, z, sigma2, mu, Y);

    /* update  pi1*/
    update_pi_l_mogc_dp(K, m2, pi_l, gamma_l);

    /* update pi2*/
    update_pi_g_mogc_dp(P, m2, m1, K, BernoulliWeighted_int, MH_ind,
    pi_g_prop_n, gl, l_index_g, weight_for_Dg, gamma_g, pi_g,
    pi_g_loglikeli);

    update_pi_j_mogc_dp(P_dup, m1, K, BernoulliWeighted_int, MH_ind, 
    pi_j_prop_n, pg, g_index, weight_for_Tj, gamma_j, pi_j,
    pi_j_loglikeli);

    /* update s2*/
    if(!s2_fixed){
      update_s2_mogc_dp(P_dup, K, num_t, uni_types, types, feature_dup_index, s2, b);
    }

    /* update cluster index */
    if(!fix_z){
      update_z_mogc_dp(N, P, K, z, pi_k, sigma2, mu, Y, Prob_Assign, debug);
    }

    /* get nk */
    nk = get_nk_mogc_dp(N, K, z);

    /* update pi_k */
    update_v_mogc_dp(K, alpha_dp, nk, v);
    update_pi_k_mogc_dp(K, v, pi_k);

    alpha_dp=update_alpha_dp_mogc_dp(K, alpha_dp, eta1, eta2, v);

    /* print log */
    int f_nsim = nsim + 1;
    if(f_nsim % print_int == 0){
      cout<< "Loop counter value is " << f_nsim << ".\n";
      cout<< "s2 is " << s2 << ".\n";
    }

    /* store results after burn-in*/
    if(nsim > burnInIter - 1){
      int store_index = nsim - burnInIter;

      Z_store.row(store_index) = z;
      for(int n = 0; n < N; n++){
        for(int k = 0; k < K; k++){
          Prob_Assign_store(store_index, n + N*k) = Prob_Assign(n, k);
        }
      }

      NK_store.row(store_index) = nk;

      PI_K_store.row(store_index) = pi_k;

      for(int k = 0; k < K; k++){
        for(int pd = 0; pd < P_dup; pd++){
          GAMMA_J_store(store_index, pd + P_dup*k) = gamma_j(k, pd);
        }
      }
      for(int k = 0; k < K; k++){
        for(int g = 0; g < m1; g++){
          GAMMA_G_store(store_index, g + m1*k) = gamma_g(k, g);
        }
      }
      for(int k = 0; k < K; k++){
        for(int l = 0; l < m2; l++){
          GAMMA_L_store(store_index, l + m2*k) = gamma_l(k, l);
        }
      }

      for(int k = 0; k < K; k++){
        for(int pd = 0; pd < P_dup; pd++){
          B_store(store_index, pd + P_dup*k) = b(k, pd);
        }
      }

      for(int k = 0; k < K; k++){
        for(int p = 0; p < P; p++){
          int index = p + P * k; 
          MU_store(store_index, index) = mu(k, p);
        }
      }

      SIGMA2_store.row(store_index) = sigma2;

      PI_L_store.row(store_index) = pi_l;

      for(int k = 0; k < K; k++){
        for(int l = 0; l < m2; l++){
          PI_G_store(store_index, l + m2*k) = pi_g(k, l);
        }

        for(int g = 0; g < m1; g++){
          PI_J_store(store_index, g + m1*k) = pi_j(k, g);
        }
      }

      S2_store.col(store_index) = s2;
      ALPHA_DP_store[store_index] = alpha_dp;
    }

  }

  return List::create(Rcpp::Named("Z_store") = Z_store,
                            Rcpp::Named("MU_store") = MU_store,
                            Rcpp::Named("PI_K_store") = PI_K_store,
                            Rcpp::Named("S2_store") = S2_store,
                            Rcpp::Named("ALPHA_DP_store") = ALPHA_DP_store,
                            Rcpp::Named("NK_store") = NK_store,
                            Rcpp::Named("GAMMA_J_store") = GAMMA_J_store,
                            Rcpp::Named("GAMMA_G_store") = GAMMA_G_store,
                            Rcpp::Named("GAMMA_L_store") = GAMMA_L_store,
                            Rcpp::Named("B_store") = B_store,
                            Rcpp::Named("SIGMA2_store") = SIGMA2_store,
                            Rcpp::Named("PI_J_store") = PI_J_store,
                            Rcpp::Named("PI_G_store") = PI_G_store,
                            Rcpp::Named("PI_L_store") = PI_L_store,
                            Rcpp::Named("Prob_Assign_store") = Prob_Assign_store,
                            Rcpp::Named("pi_g_loglikeli") = pi_g_loglikeli,
                            Rcpp::Named("pi_j_loglikeli") = pi_j_loglikeli
                            ); 
}












