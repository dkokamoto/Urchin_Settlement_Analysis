data { 
  // input information 
  int NO; // number of observations
  int NS; // number of sites
  int NM; // number of months
  int NYM; // number of months
  int NP; // number of variables
  // input data 
  int YSP[NO]; // observations (on the original scale)
  int n[NO]; // subsample count (on the original scale)
  int N[NO]; // total count (on the original scale)
  int OBS_MONTH[NO]; // continuous month-site observation index 
  int OBS_SITE[NO]; // month-site observation index 
  int PRED_MONTH[NM]; // month-site index associated with MA
  matrix[NM*NS*2,(NP)*NS*2] X; // covariate matrix
  int xind[NM,NS,2]; // covariate matrix indices
  vector[NO] D; // number of days 
  vector[NO] NB; // number of brushes 
  cholesky_factor_cov[NYM] D_seas;   // seasonal cholesky covariance factor
}

parameters {
  // process variables 
  vector[(NP)*NS*2] beta;
  
  // mean volatility
  vector<lower=0>[NS] sigma_S_mu[2]; 
  
  // spatial volatility correlation matrix
  cholesky_factor_corr[NS] L_Omega_S[2]; // 
  
  // estimated abundance
  vector<lower= -7>[NS] S[2,NM];

  // seasonal and annual variance
  real<lower = 0> sigma[2,NS];
  
  // iid errors for seasonality
  vector<lower= -7>[NYM] z_s[2,NS];
}

transformed parameters { 
  vector[NYM] LS[2,NS]; // log scale expected sp per brush, per day for every month
  
  vector[NS] LS2[2,NM]; // log scale expected sp per brush, per day for every month
  
  vector[NO] S_exp[2]; // normal scale estimated sp for each observation
  
  vector[2*NM*NS] s_hat;
  
  for(q in 1:2){ 
    for (i in 1:NS) {
      LS[q,i]  = D_seas*z_s[q,i];
    }
  }
  
  s_hat = X*beta;
  
  // put it all together
  for (q in 1:2){
    for(i in 1:NM){
      for(j in 1:NS){
        LS2[q,i,j]=LS[q,j,PRED_MONTH[i]]*sigma[q,j]+s_hat[xind[i,j,q]];
      }
    }
  
    for(i in 1:NO){
      S_exp[q,i] = exp(S[q,OBS_MONTH[i],OBS_SITE[i]])*D[i];
    }
  }
}

model { 
  for(q in 1:2){
    for(i in 1:NS){
      sigma[q,i] ~ cauchy(0, 0.5);
    }
    // process corr matrix 
    L_Omega_S[q] ~ lkj_corr_cholesky(1);

    // seasonal GP error
    for(i in 1:NS){
      z_s[q,i]~ normal(0,1);
    }
    
     // process likelihood
      S[q] ~ multi_normal_cholesky(LS2[q], diag_pre_multiply(sigma_S_mu[q],L_Omega_S[q]));
  
    sigma_S_mu[q] ~ cauchy(0,2.5);
  }
  for(i in 1:(NP*NS*2)){
    if(i<(NS*2+1)){
     beta[i] ~ cauchy(0,10); //prior for the intercept following Gelman 2008
    } else {
      beta[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
    }
  }
  // observation likelihood 
  target+= binomial_lpmf(YSP|n,S_exp[2]./ (S_exp[1]+S_exp[2])); // define the likelihood
  target+= poisson_lpmf(N|S_exp[1]+S_exp[2]);
}

generated quantities {
  corr_matrix[NS] Omega_S[2];
  for(q in 1:2){
    Omega_S[q] = multiply_lower_tri_self_transpose(L_Omega_S[q]);
  }
}
