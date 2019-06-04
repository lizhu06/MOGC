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
void set_seed_sogc_dp(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export]]
void update_z_sogc_dp(int N, int P, int K, Eigen::VectorXi& z,
  Eigen::VectorXd& pi_c, Eigen::VectorXd& sigma2, 
  Eigen::MatrixXd& mu, Eigen::MatrixXd& Y, Eigen::MatrixXd& Prob_Assign, bool debug) {

  //MatrixXi Prob_Assign_store(N, K);
  double resi = (-0.5)*log(2*3.1415);
  for(int i = 0; i < N; i++){
    NumericVector logden_c(K);
    NumericVector un_norm_prob_c(K);
    NumericVector norm_prob_c(K);

    for(int k = 0; k < K; k++){
      double logden = 0.0;
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
      un_norm_prob_c[k] = pi_c[k] * exp(logden_c[k] - logden_offset);
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
Eigen::VectorXi get_nk_sogc_dp(int N, int K,  Eigen::VectorXi& z){
  Eigen::VectorXi nk = Eigen::VectorXi::Zero(K);
  for(int i = 0; i < N; i++){
    nk[z[i]-1] = nk[z[i]-1] + 1;
  }
  return(nk);
}

// [[Rcpp::export]]
void update_v_sogc_dp(int K, double alpha_dp, Eigen::VectorXi& nk, Eigen::VectorXd& v){
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
void update_pi_c_sogc_dp(int K, Eigen::VectorXd& v, Eigen::VectorXd& pi_c){
  pi_c[0] = v[0];
  for(int k=0; k<K; k++){
    double pi_c_prod = 1-v[0];
    for(int k2=0; k2<k; k2++){
      pi_c_prod = pi_c_prod*(1-v[k2]);
    }
    pi_c[k] = pi_c_prod*v[k];
  }
}

// [[Rcpp::export]]
double update_alpha_dp_sogc_dp(int K, double alpha_dp, double eta1, double eta2, Eigen::VectorXd& v){
  double log_sum_v = 0.0;
  for(int k=0; k<K-1; k++){
    log_sum_v = log_sum_v + log(1-v[k]);
  }
  log_sum_v = log_sum_v + pow(10.0, -10.0);
  alpha_dp = R::rgamma(K + eta1 - 1, eta2 - log_sum_v);
  return(alpha_dp);
}

// [[Rcpp::export]]
Eigen::MatrixXd update_mu_sogc_dp(int P, int P_dup, int K, Eigen::VectorXi& g_index,
  Eigen::VectorXi& feature_dup_index,
  Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_j, 
  Eigen::MatrixXd& b){

  Eigen::MatrixXd mu = Eigen::MatrixXd::Zero(K, P);
  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int p = feature_dup_index[pd] - 1;
      int g = g_index[pd] - 1;
      //cout<< "mu " << mu(k, p) << ".\n";
      mu(k, p) = mu(k, p) + gamma_g(k, g) * gamma_j(k, pd) *b(k, pd);
    }
  }
  return(mu);
}

// [[Rcpp::export]]
void update_gamma_g_sogc_dp(int P_dup, int N, int K, int G, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index, Eigen::VectorXi& U1_rowsum, 
  Eigen::VectorXd& pi_g, Eigen::VectorXd& sigma2,
  Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_j, 
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  /*cout<< "gamma_g " << gamma_g.row(0) << ".\n";*/
  for(int k = 0; k < K; k++){
    for(int g = 0; g < G; g++){
      int old_gamma_g = gamma_g(k,g);
      double inside = 0.0;
      for(int pd = 0; pd < P_dup; pd++){
        int p = feature_dup_index[pd] - 1; 
        if(g_index[pd] == g + 1){
          double sumy = 0.0;
          for(int i = 0; i < N; i++){
            if(z[i] == k+1){
              if(U1_rowsum[p]==1){
                sumy = sumy + Y(i, p);
              }else{
                sumy = sumy + Y(i, p) - mu(k, p) + gamma_g(k, g) * gamma_j(k, pd) * b(k, pd);
              }
            }
          }
          inside = inside + (nk[k] * pow(gamma_j(k, pd) * b(k, pd), 2) 
            - 2.0 * gamma_j(k, pd) * b(k, pd) * sumy) / (2.0 * sigma2[p]);
        }
      }
      //cout<< "k=" << k << ", g=" << g << ", pi_g(k)= " << pi_g[k] << ".\n";
      //cout<< "inside= " << inside << ".\n";
      double p_gamma_g = 1.0 / (1.0 + (1.0 - pi_g[k])/(pi_g[k]) * 
        exp(inside));
      //cout<< "p_gamma_g " << p_gamma_g << ".\n";
      gamma_g(k, g) = R::rbinom(1, p_gamma_g);
      //cout<< "gamma_g " << gamma_g(k, g) << ".\n";

       // if gamma_g change, update mu
      if(gamma_g(k,g) != old_gamma_g){
       for(int pd = 0; pd < P_dup; pd++){
         int p = feature_dup_index[pd] - 1;
         int g_2 = g_index[pd] - 1;
         if(g_2==g){
           mu(k, p) = mu(k, p) - old_gamma_g *
             gamma_j(k, pd)*b(k, pd) + 
             gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
         }
       }
      }

    }
  }
  /*cout<< "gamma_g " << gamma_g.row(0) << ".\n";*/
  //return(gamma_g);
}


// [[Rcpp::export]]
void update_gamma_j_sogc_dp(int P_dup, int N, int K, int BernoulliWeighted_int, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& pg, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index, Eigen::VectorXi& U1_rowsum, 
  Eigen::VectorXd& sigma2, Eigen::VectorXd& weight_for_Rj,
  Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_j, 
  Eigen::MatrixXd& pi_j, Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){

  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int old_gamma_j = gamma_j(k,pd);
      int p = feature_dup_index[pd] - 1; 
      int g = g_index[pd] - 1;
      //int group_index = g_index[p];
      //if(pg[group_index-1] > 1){  // singleton still has two gammas
        double sumy = 0.0;
        for(int i = 0; i < N; i++){
          if(z[i] == k + 1){
            if(U1_rowsum[p] == 1){
              sumy = sumy + Y(i, p);
            }else{
              sumy = sumy + Y(i, p) - mu(k, p) + gamma_g(k, g) * gamma_j(k, pd) * b(k, pd);
            }
          }
        }
        double inside = (nk[k] * pow(gamma_g(k, g) * b(k, pd), 2) - 
          2.0 * gamma_g(k, g) * 
          b(k, pd) * sumy)/(2.0 * sigma2[p]);

        double pi_j_weighted = pi_j(k, g);
        if(BernoulliWeighted_int == 1){
          pi_j_weighted = pi_j(k, g)*weight_for_Rj[pd];
        }

        double p_gamma_j = 1.0/(1.0 + (1.0 - pi_j_weighted)/
          pi_j_weighted * exp(inside));

        gamma_j(k, pd) = R::rbinom(1, p_gamma_j);

         // if gamma_g change, update mu
        if(gamma_j(k,pd) != old_gamma_j){
          mu(k, p) = mu(k, p) - old_gamma_j *
            gamma_g(k, g)*b(k, pd) + 
            gamma_g(k, g)*gamma_j(k, pd)*b(k, pd);
        }
      //}
    }
  }
  //return(gamma_j);
}


// [[Rcpp::export]]
void update_b_sogc_dp(int P_dup, int N, int K,  
  Eigen::VectorXd& s2,
  Eigen::VectorXi& g_index, Eigen::VectorXi& nk, Eigen::VectorXi& z, 
  Eigen::VectorXi& feature_dup_index, Eigen::VectorXi& U1_rowsum, Eigen::VectorXi& types, 
  Eigen::VectorXd& sigma2,
  Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_j,
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y){
 
  for(int k = 0; k < K; k++){
    for(int pd = 0; pd < P_dup; pd++){
      int p = feature_dup_index[pd] - 1; 
      int g = g_index[pd] - 1;

      int a_type = types[p] - 1;

      if(gamma_g(k, g) * gamma_j(k, pd) == 0){
        b(k,pd) = R::rnorm(0.0, sqrt(s2[a_type]));
      }else{

        double var_b = sigma2[p] * s2[a_type] / (s2[a_type] * nk[k] + sigma2[p]);
        double sumy = 0.0;
        double mu_remove_this_b = mu(k, p) - b(k,pd);
        for(int i = 0; i < N; i++){
          if(z[i] == k + 1){
            if(U1_rowsum[p]==1){
              sumy = sumy + Y(i, p);
            }else{
              sumy = sumy + Y(i, p) - mu(k, p) + gamma_g(k, g) * gamma_j(k, pd) * b(k, pd);
            }
          }
        }
        double mean_b = sumy * var_b / sigma2[p];
        b(k, pd) =  R::rnorm(mean_b, sqrt(var_b));
        mu(k, p) = mu_remove_this_b + b(k, pd);
        //cout<< "b=" << b(k, pd)   << ".\n";
      }
    }
  }
  //return(b);
}

// [[Rcpp::export]]
void update_sigma2_sogc_dp(int P, int N, Eigen::VectorXi& z, Eigen::VectorXd& sigma2, 
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
void update_pi_g_sogc_dp(int K, int G, 
  Eigen::VectorXd& pi_g, Eigen::MatrixXi& gamma_g){

  for(int k = 0; k < K; k++){
    double sum_gamma_g = gamma_g.row(k).sum();
    pi_g[k] = R::rbeta(sum_gamma_g + 1.0, G - sum_gamma_g + 1.0);
  }
  //return(pi_g);
}

// [[Rcpp::export]]
void update_pi_j_sogc_dp(int P_dup, int G, int K, int BernoulliWeighted_int,
  int MH_ind, double pi_j_prop_n, 
  Eigen::VectorXi& pg, Eigen::VectorXi& g_index, Eigen::VectorXd& weight_for_Rj,
  Eigen::MatrixXi& gamma_j,
  Eigen::MatrixXd& pi_j, Eigen::MatrixXd& pi_j_loglikeli){

  for(int k = 0; k < K; k++){
    for(int g = 0; g < G; g++){
      double fg = g + 1;

      if(MH_ind == 1){
        double loglikeli_new=0.0;
        double pi_j_new = R::rbeta(pi_j_prop_n * pi_j(k, g), 
          pi_j_prop_n * (1.0 - pi_j(k, g)));

        for(int pd = 0; pd < P_dup; pd++){
          if(g_index[pd] == fg){
            if(gamma_j(k,pd)==1){
              loglikeli_new = loglikeli_new + 
                log(pi_j_new * weight_for_Rj[pd]);
            }else{
              loglikeli_new = loglikeli_new + 
                log(1-(pi_j_new * weight_for_Rj[pd]));
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
          if(g_index[pd] == fg){
            sum_gamma_j = sum_gamma_j + gamma_j(k, pd);
          }
        }
        pi_j(k, g) = R::rbeta(sum_gamma_j + 1.0, pg[g] - sum_gamma_j + 1.0);
      }
    }
  }
  //return(pi_j);//
}

// [[Rcpp::export]]
void update_s2_sogc_dp(int P_dup, int K, int T, Eigen::VectorXi& uni_types, 
  Eigen::VectorXi& types, Eigen::VectorXi& feature_dup_index, 
  Eigen::VectorXd& s2, Eigen::MatrixXd& b){

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
List MCMC_sogc_dp(int seed, int burnInIter, int keepIter, int print_int, 
  int N, int P, int P_dup, int G, 
  int K, int BernoulliWeighted_int, int MH_ind, 
  Eigen::VectorXd& s2, double pi_j_prop_n,
  Eigen::VectorXi& uni_types, Eigen::VectorXi& types,
  double eta1, double eta2, double alpha_dp, 
  Eigen::VectorXi& g_index, Eigen::VectorXi& pg, Eigen::VectorXi& z, Eigen::VectorXi& nk, 
  Eigen::VectorXi& feature_dup_index, Eigen::VectorXi& U1_rowsum, 
  Eigen::VectorXd& v, 
  Eigen::VectorXd& pi_c, Eigen::VectorXd& sigma2, Eigen::VectorXd& weight_for_Rj, 
  Eigen::MatrixXi& gamma_g, Eigen::MatrixXi& gamma_j, 
  Eigen::VectorXd& pi_g, Eigen::MatrixXd& pi_j, 
  Eigen::MatrixXd& b, Eigen::MatrixXd& mu, Eigen::MatrixXd& Y, bool debug, 
  bool fix_z, bool fix_mu, Eigen::MatrixXd& pi_j_loglikeli, bool s2_fixed){

  int NSIM = burnInIter + keepIter;
  /* store matrix */
  Eigen::MatrixXi Z_store(keepIter, N);
  int nrow = K * P;
  Eigen::MatrixXd MU_store(keepIter, nrow);
  Eigen::MatrixXd PI_C_store(keepIter, K);
  Eigen::MatrixXd Prob_Assign(N, K);
  nrow = N * K;
  Eigen::MatrixXd Prob_Assign_store(keepIter, nrow);

  int num_t = uni_types.size();
  
  Eigen::MatrixXd S2_store(num_t, keepIter);
  Eigen::MatrixXi NK_store(keepIter, K);
  Eigen::MatrixXi gamma_g_store(keepIter, G*K);
  Eigen::MatrixXi gamma_j_store(keepIter, P_dup*K);
  Eigen::MatrixXd B_store(keepIter, P_dup*K);
  Eigen::MatrixXd SIGMA2_store(keepIter, P);
  Eigen::MatrixXd pi_g_store(keepIter, K);
  Eigen::MatrixXd pi_j_store(keepIter, K*G);

  set_seed_sogc_dp(seed);

  /* MCMC */
  for(int nsim = 0; nsim < NSIM; nsim++){
    
    /* update cluster index */
    if(!fix_z){
      update_z_sogc_dp(N, P, K, z, pi_c, sigma2, mu, Y, Prob_Assign, debug);
    }

    /* get nk */
    nk = get_nk_sogc_dp(N, K, z);

    /* update pi_c */
    update_v_sogc_dp(K, alpha_dp, nk, v);
    update_pi_c_sogc_dp(K, v, pi_c);

    alpha_dp=update_alpha_dp_sogc_dp(K, alpha_dp, eta1, eta2, v);

    if(!fix_mu){
      /* update gamma_g */
      update_gamma_g_sogc_dp(P_dup, N, K, G, g_index, nk, 
        z, feature_dup_index, U1_rowsum, pi_g, sigma2, 
        gamma_g, gamma_j, b, mu, Y);

      /* update gamma_j */
      update_gamma_j_sogc_dp(P_dup, N, K, BernoulliWeighted_int, 
        g_index, pg, nk, z, feature_dup_index, U1_rowsum,
        sigma2, weight_for_Rj, gamma_g, gamma_j, pi_j, b, mu, Y);

      /* update b */
      update_b_sogc_dp(P_dup, N, K, s2, g_index, 
        nk, z, feature_dup_index, U1_rowsum, types,
        sigma2, gamma_g, gamma_j, b, mu, Y);

      /* update mu*/
      //mu = update_mu(P, P_dup, K, g_index, feature_dup_index,
      //  gamma_g, gamma_j, b);

    }

    /* update sigma2 */
    update_sigma2_sogc_dp(P, N, z, sigma2, mu, Y);

    /* update  pi_g*/
    update_pi_g_sogc_dp(K, G, pi_g, gamma_g);

    /* update pi_j*/
    update_pi_j_sogc_dp(P_dup, G, K, BernoulliWeighted_int, MH_ind, pi_j_prop_n,
      pg, g_index, weight_for_Rj, gamma_j, pi_j, pi_j_loglikeli);

    /* update s2*/
    if(!s2_fixed){
      update_s2_sogc_dp(P_dup, K, num_t, uni_types, types, feature_dup_index, s2, b);
    }

    int f_nsim = nsim + 1;
    if(f_nsim%print_int == 0){
      cout<< "Loop counter value is " << f_nsim << ".\n";
      cout<< "s2 is " << s2 << ".\n";
    }

    /* store results after burn-in*/
    if(nsim > burnInIter - 1 ){

      int store_index = nsim - burnInIter;

      Z_store.row(store_index) = z;
      for(int n = 0; n < N; n++){
        for(int k = 0; k < K; k++){
          Prob_Assign_store(store_index, n + N*k) = Prob_Assign(n, k);
        }
      }

      NK_store.row(store_index) = nk;

      PI_C_store.row(store_index) = pi_c;

      for(int k = 0; k < K; k++){
        for(int g = 0; g < G; g++){
          gamma_g_store(store_index, g + G*k) = gamma_g(k, g);
        }
      }
      for(int k = 0; k < K; k++){
        for(int pd = 0; pd < P_dup; pd++){
          gamma_j_store(store_index, pd + P_dup*k) = gamma_j(k, pd);
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

      pi_g_store.row(store_index) = pi_g;

      for(int k = 0; k < K; k++){
        for(int g = 0; g < G; g++){
          pi_j_store(store_index, g + G*k) = pi_j(k, g);
        }
      }

      S2_store.col(store_index) = s2;
    }
  }

  return List::create(Rcpp::Named("Z_store") = Z_store,
                            Rcpp::Named("MU_store") = MU_store,
                            Rcpp::Named("PI_C_store") = PI_C_store,
                            Rcpp::Named("S2_store") = S2_store,
                            Rcpp::Named("NK_store") = NK_store,
                            Rcpp::Named("GAMMA1_store") = gamma_g_store,
                            Rcpp::Named("GAMMA2_store") = gamma_j_store,
                            Rcpp::Named("B_store") = B_store,
                            Rcpp::Named("SIGMA2_store") = SIGMA2_store,
                            Rcpp::Named("pi_g_store") = pi_g_store,
                            Rcpp::Named("pi_j_store") = pi_j_store,
                            Rcpp::Named("Prob_Assign_store") = Prob_Assign_store,
                            Rcpp::Named("pi_j_loglikeli") = pi_j_loglikeli,
                            Rcpp::Named("v") = v
                            ); 
}












