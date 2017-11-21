data { 
  // input information 
  int NO; // number of observations
  int NS; // number of sites
  int NM; // number of months
  int NYM; // number of months
  // input data 
  int YSP[NO]; // observations (on the original scale)
  int n[NO]; // subsample count (on the original scale)
  int N[NO]; // total count (on the original scale)
  int OBS_MONTH[NO]; // continuous month-site observation index associated with MA
  int OBS_SITE[NO]; // month-site observation index associated with MA
  int PRED_MONTH[NM]; // month-site observation index associated with MA
  vector[NO] D; // number of days 
  vector[NO] NB; // number of brushes 
  cholesky_factor_cov[NYM] D_seas;   // seasonal cholesky covariance factor
  cholesky_factor_cov[NM] D_ann;   // annual cholesky covariance factor 
}

parameters {
  // mean volatility
  vector<lower=1e-7>[NS] sigma_S_mu[2]; 
  
  // spatial volatility correlation matrix
  cholesky_factor_corr[NS] L_Omega_S[2]; // 
  cholesky_factor_corr[NS] L_Omega_Spatial[2];
  
  // estimated abundance
  vector<lower= -10,upper= 10>[NS] S[2,NM];

  // seasonal and annual variance
  real<lower = 0> sigma[2,NS];
  real<lower = 0> sigma2[2,NS];
  
  // mean site settlement
  vector<lower= -10, upper= 10>[NS] mu_S[2];
  
  // iid errors for seasonality
  vector<lower= -10, upper= 10>[NYM] z_s[2,NS];
  
  // iid errors for annual trends
  vector<lower= -10, upper= 10>[NS] w_z[2,NM];
  vector<lower= -10, upper= 10>[NS] e_z[2,NM];
}

transformed parameters { 
  vector[NYM] LS[2,NS]; // log scale expected sp per brush, per day for every month
  
  vector[NS] LS2[2,NM]; // log scale expected sp per brush, per day for every month
  
  vector[NO] S_exp[2]; // normal scale estimated sp for each observation
  
  vector[NO] N_mu; // expected total count
  vector[NO] q_mu; // expected fraction of N that is SP
  
  vector[NM] w[2,NS];
    
  vector[NS] w_z2[2,NM];
  vector[NM] w_z3[2,NS];
  vector[NS] w_e[2,NM];
  
  for(q in 1:2){ 
    for(i in 1:NM){
      w_z2[q,i] = L_Omega_Spatial[q]*w_z[q,i];
    }
    for(i in 1:NM){
      for(j in 1:NS){
        w_z3[q,j,i]=w_z2[q,i,j];
      }
      w_e[q,i] = diag_pre_multiply(sigma_S_mu[q],L_Omega_S[q])*e_z[q,i];
    }
    for (i in 1:NS) {
      LS[q,i]  = D_seas*z_s[q,i];
      w[q,i] = D_ann*w_z3[q,i];
    }
  }
  
  // put it all together
  for (q in 1:2){
    for(i in 1:NM){
      for (j in 1:NS){
        LS2[q,i,j]=mu_S[q,j]+LS[q,j,PRED_MONTH[i]]*sigma[q,j]+w[q,j,i]*sigma2[q,j]+w_e[q,i,j];
      }
    }
    for(i in 1:NO){
      S_exp[q,i] = exp(S[q,OBS_MONTH[i],OBS_SITE[i]])*NB[i]*D[i];
    }
  }
  
  N_mu = S_exp[1]+S_exp[2];
  q_mu = S_exp[1] ./ N_mu;
}

model { 
  for(q in 1:2){
    for(i in 1:NS){
      sigma[q,i] ~ cauchy(0, 0.5);
      sigma2[q,i] ~ cauchy(0, 0.5);
    }

    L_Omega_S[q] ~ lkj_corr_cholesky(1);
    L_Omega_Spatial[q] ~ lkj_corr_cholesky(1);
  
    for(i in 1:NS){
      z_s[q,i]~ normal(0,1);
    }
    for(i in 1:NM){
      w_z[q,i]~ normal(0,1);
      e_z[q,i]~ normal(0,1);
    }
  
    sigma_S_mu[q] ~ cauchy(0,2.5);
  
    mu_S[q] ~ normal(0,5);
  }
    
  // observation likelihood 
  target+= binomial_lpmf(YSP|n,q_mu); // define the likelihood
  target+= poisson_lpmf(N|N_mu);
}

generated quantities {
  corr_matrix[NS] Omega_S[2];
  corr_matrix[NS] Omega_Spatial[2];
  for(q in 1:2){
    Omega_S[q] = multiply_lower_tri_self_transpose(L_Omega_S[q]);
    Omega_Spatial[q] = multiply_lower_tri_self_transpose(L_Omega_Spatial[q]);
  }
}
