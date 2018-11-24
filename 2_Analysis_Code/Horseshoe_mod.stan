
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
  int OBS[NO]; // continuous month-site observation index 
  int NP; // number of variables
  // input data 
  int YSP[NO]; // observations (on the original scale)
  int n[NO]; // subsample count (on the original scale)
  int N[NO]; // total count (on the original scale)
  int OBS_MONTH[NO]; // continuous month-site observation index associated with MA
  int OBS_SITE[NO]; // month-site observation index associated with MA
  int PRED_MONTH[NM]; // month-site observation index associated with MA
  matrix[NM*NS,NP] X; // covariate matrix
  int xind[NM,NS]; // covariate matrix indices
  vector[NO] D; // number of days 
  cholesky_factor_cov[NYM] D_seas;   // seasonal cholesky covariance factor
}

parameters {
  // mean volatility
  vector<lower=1e-7>[NS] sigma_S_mu; 
  
  // spatial volatility correlation matrix
  cholesky_factor_corr[NS] L_Omega_S; // 
  cholesky_factor_corr[NS] L_Omega_Spatial_seas;
  
  // iid errors for seasonality
  matrix<lower= -10, upper= 10>[NYM,NS] z_s;

  // seasonal and annual variance
  real<lower = 0> sigma[NS];
  
  // process variables 
  vector[NS] beta0;
  
  // iid errors for residuals
  matrix<lower= -4, upper= 4>[NS,NM] e_z;
  
  // proportion of samples that are purps
  vector<lower= 0,upper= 1>[NSUB] theta;
  
  // iid errors for seasonal gaussian process
  vector [NP] z;
  real < lower =0 > r1_global ;
  real < lower =0 > r2_global ;
  vector < lower =0 >[NP] r1_local ;
  vector < lower =0 >[NP] r2_local ;
}

transformed parameters { 
  matrix[NYM,NS] LS; // log scale expected sp per brush, per day for every month
  matrix[NM,NS] S; // estimated abundance
  vector[NO] S_exp; // normal scale estimated sp for each observation
  
  vector[NO] N_mu; // expected total count
  vector[NO] q_mu; // expected fraction of N that is SP
  
  vector[NM*NS] s_hat;
  
  // horseshoe hyperparameters
  vector<lower=0>[NP] lambda_h;
  real<lower=0> tau_h;
  
  // process variables 
  vector[NP] beta;
  vector[NO] theta2; // fraction that is SP
  
  // correlated errors 
  matrix[NM,NS] e;
  
  // horsehoe
  lambda_h = r1_local .* sqrt( r2_local );
  tau_h = r1_global * sqrt( r2_global );
  beta = z .* lambda_h * tau_h ;
  
  // regression model 
  s_hat = X*beta;
  
  // seasonal gaussian process model 
  LS  = D_seas*z_s;

  // process error
  e = transpose(diag_pre_multiply(sigma_S_mu,L_Omega_S)*e_z);  
  
  // put it all together
  for(i in 1:NM){
    for (j in 1:NS){
      S[i,j]=beta0[j]+s_hat[xind[i,j]]+LS[PRED_MONTH[i],j]*sigma[j]-5+e[i,j];
    }
  }

  theta2 = rep_vector(1.0,NO);
  theta2[SUB] = theta;
  
  S_exp = exp(to_vector(S)[OBS]) .*D ./theta2;
}

model { 
  sigma ~ cauchy(0, 0.5);
  
  L_Omega_S ~ lkj_corr_cholesky(2);

  for(i in 1:NS){
    z_s[i]~ normal(0,1);
    e_z[i] ~normal(0,1);
  }
  
  sigma_S_mu ~ cauchy(0,2.5);
  
  beta0 ~ cauchy(0,10); //prior for the intercept following Gelman 2008
 
 // iid errors for process error

  
  // iid errors for seasonal gaussian process
  z ~ normal (0 , 1);
  r1_local ~ normal (0.0 , 1.0);
  r2_local ~ inv_gamma (0.5, 0.5);
  
  // half - t prior for tau
  r1_global ~ normal (0.0 , 1);
  r2_global ~ inv_gamma (0.5 , 0.5);
  
  // jeffrey's prior given subsampling
  theta~beta(P1, P2);
  
  // observation likelihood 
  N~poisson(S_exp);
}
