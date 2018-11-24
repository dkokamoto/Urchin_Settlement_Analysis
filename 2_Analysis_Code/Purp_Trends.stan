data { 
  // input information 
  int NO; // number of observations
  int NS; // number of sites
  int NM; // number of months
  int NYM; // number of months
  int NSUB; // number of observations with subsamples
  // input data 
  //  Beta priors for uncertainty in subsamples
  vector[NSUB] P1; // SP observations (on the original scale)+0.5 
  vector[NSUB] P2; // SP observations - subsample count (on the original scale) + 0.5
  int SUB[NSUB];
    
  int N[NO]; // total count (on the original scale)
  int OBS[NO]; // continuous month-site observation index 
  int PRED_MONTH[NM]; // month-site observation index associated with MA
  vector[NO] D; // number of days 
  cholesky_factor_cov[NYM] D_seas;   // seasonal cholesky covariance factor
  cholesky_factor_cov[NM] D_ann;   // annual cholesky covariance factor 
}

parameters {
  // mean volatility
  vector<lower=1e-7>[NS] sigma_S_mu; 
  
  // spatial volatility correlation matrix
  cholesky_factor_corr[NS] L_Omega_S; // 
  cholesky_factor_corr[NS] L_Omega_Spatial;
  
  // seasonal and annual variance
  vector<lower = 0>[NS] sigma;
  vector<lower = 0>[NS] sigma2;
  
  // mean site settlement
  row_vector<lower= -10, upper= 10>[NS] mu_S;
  
  // iid errors for seasonality
  matrix<lower= -10, upper= 10>[NYM,NS] z_s;
  
  // iid errors for annual trends
  matrix<lower= -10, upper= 10>[NS,NM] w_z;
  
  // iid errors for residuals
  matrix<lower= -4, upper= 4>[NS,NM] e_z;
  
  // proportion of samples that are purps
  vector<lower= 0,upper= 1>[NSUB] theta;
}

transformed parameters { 
  matrix[NYM,NS] LS; // log scale expected sp per brush, per day for every month
  
  matrix[NM,NS] S; // estimated abundance
  
  vector[NO] S_exp; // normal scale estimated sp for each observation
  
  vector[NO] theta2; // fraction that is SP
  
  matrix[NM,NS] w;
  
  // correlated errors 
  matrix[NM,NS] e;
  
  LS  = D_seas*z_s;
  w = D_ann*transpose(L_Omega_Spatial*w_z);

  theta2 = rep_vector(1.0,NO);
  theta2[SUB] = theta;
  // put it all together
  e = transpose(diag_pre_multiply(sigma_S_mu,L_Omega_S)*e_z);
  S= rep_matrix(mu_S,NM)+diag_post_multiply(LS[PRED_MONTH,],sigma2)+diag_post_multiply(w,sigma)+e;
  S_exp = exp(to_vector(S)[OBS]) .*D ./theta2;
}

model { 
  sigma ~ cauchy(0, 0.5);
  sigma2 ~ cauchy(0, 0.5);

  L_Omega_S ~ lkj_corr_cholesky(2);
  L_Omega_Spatial ~ lkj_corr_cholesky(2);

  for(i in 1:NS){
    z_s[i]~ normal(0,1);
    mu_S[i] ~ normal(0,5);
  }
 
  sigma_S_mu ~ cauchy(0,2.5);
    
  // uncertainty on proportions of purples in the sample 
  // jeffrey's prior given subsampling
  theta~beta(P1, P2);
  
  // observation likelihood 
  N~poisson(S_exp);
}