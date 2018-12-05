
data { 
  // input information 
  int NO; // number of observations
  int NS; // number of sites
  int NM; // number of months
  int NYM; // number of months
  int NSUB; // number of observations with subsamples
  int NSP; // number of subsample counts
  int NP; // number of variables
  // input data 
  //  Beta priors for uncertainty in subsamples
  vector[NSP] P1; // SP observations (on the original scale)+0.5 
  vector[NSP] P2; // NON SP observations (on the original scale) + 0.5
  int SUB[NSUB];
  int ID[NO];
  int OBS[NO]; // continuous month-site observation index 
  // input data 
  int N[NO]; // total count (on the original scale)
  int OBS_MONTH[NO]; // continuous month-site observation index associated with MA
  int OBS_SITE[NO]; // month-site observation index associated with MA
  int PRED_MONTH[NM]; // month-site observation index associated with MA
  matrix[NM*NS,NP] X; // covariate matrix
  int xind[NM,NS]; // covariate matrix indices
  vector[NO] D; // number of days 
  cholesky_factor_cov[NYM] D_seas;   // seasonal cholesky covariance factor
  real<lower= 0> sigma_scale;
}

parameters {
  // mean volatility
  vector<lower=0>[NS] sigma_S_mu; 
  
  // spatial volatility correlation matrix
  cholesky_factor_corr[NS] L_Omega_S; // 
  cholesky_factor_corr[NS] L_Omega_Spatial_seas;
  
  // iid errors for seasonality
  matrix[NYM,NS] z_s;

  // seasonal and annual variance
  real<lower = 0> sigma[NS];
  
   // mean site settlement
  row_vector[NS] mu_S;
  
  // iid errors for residuals
  matrix<lower= -4,upper=4>[NS,NM] e_z;
  
  // AR1 parameters
  vector<lower= -0.9,upper= 0.9>[NS] phi;

  // proportion of samples that are purps
  vector<lower= 0,upper= 1>[NSP] theta;
  
  // iid errors for seasonal gaussian process
  vector[NP] beta;
}

transformed parameters { 
  matrix[NYM,NS] LS; // log scale expected sp per brush, per day for every month
  matrix[NM,NS] S; // estimated abundance
  vector[NO] S_exp; // normal scale estimated sp for each observation
  
  vector[NO] N_mu; // expected total count
  vector[NO] q_mu; // expected fraction of N that is SP
  vector[NS] sigma_star;
  
  vector[NM*NS] s_hat;

  // process variables 
  vector[NO] theta2; // fraction that is SP
  
  // correlated errors 
  matrix[NM,NS] e;
  matrix[NM,NS] e_hat;
  
  // regression model 
  s_hat = X*beta;
  
  // seasonal gaussian process model 
  LS  = D_seas*z_s;
  
  // process error
  for(i in 1:NS){
    sigma_star[i] = sigma_S_mu[i]*sqrt(1.0-phi[i]^2.0);
  } 
  e = transpose(diag_pre_multiply(sigma_star,L_Omega_S)*e_z);  
  // put it all together

  for(j in 1:NS){
    e_hat[1,j] = e[1,j];
    S[1,j]=mu_S[j]+LS[PRED_MONTH[1],j]*sigma[j]+s_hat[xind[1,j]]+e_hat[1,j]-5;
    for(i in 2:NM){
      e_hat[i,j] =  phi[j]*e_hat[i-1,j]+e[i,j];
      S[i,j]=mu_S[j]+LS[PRED_MONTH[i],j]*sigma[j]+s_hat[xind[i,j]]+e_hat[i,j]-5;
    }
  }
  
  theta2 = rep_vector(1.0,NO);
  theta2[SUB] = theta[ID[SUB]];
  
  S_exp = exp(to_vector(S)[OBS]) .*D ./theta2;
}

model { 
  sigma ~ normal(0.0, 2);
  
  L_Omega_S ~ lkj_corr_cholesky(2.0);

  for(i in 1:NS){
    z_s[i]~ normal(0.0,1.0);
    e_z[i]~ normal(0.0,1.0);
    mu_S[i] ~ cauchy(0.0,10); 
  }
  
  sigma_S_mu ~ normal(0.0,sigma_scale);

  beta ~ normal(0, 10);

  // jeffrey's prior given subsampling
  theta~beta(P1, P2);
  
  // observation likelihood 
  N~poisson(S_exp);
}
