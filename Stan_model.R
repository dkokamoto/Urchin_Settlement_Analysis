model <- 
  'data { 
    // input information 
    int NO; // number of observations
    int NP; // number of parameters
    int NS; // number of sites
    int NM; // number of months
    // input data 
    int YSP[NO]; // observations (on the original scale)
    int n[NO]; // subsample count (on the original scale)
    int N[NO]; // total count (on the original scale)
    int OBS_MONTH[NO]; // month-site observation index associated with MA
    int OBS_SITE[NO]; // month-site observation index associated with MA
    vector[NO] D; // number of days 
    vector[NO] NB; // number of brushes  
    matrix[NS,NP] MM[NM]; // design matrix
  }
    
  parameters {
    // process variables 
    vector<lower= -10,upper= 10>[NP] beta_SP;
    
    vector<lower=0>[NS] L_sigma_SP;  
    cholesky_factor_corr[NS] L_Omega_SP; // 
    
    vector<lower= -10,upper= 10>[NS] SP[NM];
    vector<lower= -10,upper= 10>[NS] SF[NM]; 
  }
  
  transformed parameters { 
    vector[NS] LSP[NM]; // log scale expected sp per brush, per day for every month
    
    vector[NO] SP_exp; // normal scale estimated sp for each observation
    vector[NO] SF_exp; // normal scale estimated sp for each observation
  
    vector[NO] N_mu; // expected total count
    vector[NO] q_mu; // expected fraction of N that is SP

    matrix[NS,NS] L_Sigma_SP;

    L_Sigma_SP <- diag_post_multiply(L_Omega_SP, L_sigma_SP);  

    for (i in 1:NM){
      LSP[i] <- MM[i]*beta_SP;
    }
  
    for(i in 1:NO){
      SP_exp[i] <- exp(SP[OBS_MONTH[i],OBS_SITE[i]])*NB[i]*D[i];
      SF_exp[i] <- exp(SF[OBS_MONTH[i],OBS_SITE[i]])*NB[i]*D[i];
    }
    
      N_mu <- SF_exp+SP_exp; // expected total count is sum of expected total SP & SF
      q_mu <- SP_exp ./ N_mu;// expected proportion is expected number of SP divided by expected total count
    }

  model { 
    L_Omega_SP ~ lkj_corr_cholesky(NS);  //prior for cholesky factor for correlation matrix
    L_sigma_SP ~ cauchy(0, 2.5); //prior for variance
    
    for (i in 1:NM) {
      SP[i] ~ multi_normal_cholesky(LSP[i], L_Sigma_SP);  //multivariate normal prior for each obs
    }

    YSP ~ binomial(n,q_mu); // define the likelihoods
    N ~ poisson(N_mu); 
  }
  generated quantities {
    corr_matrix[NS] Omega;  //produce correlation matrix
    Omega <- multiply_lower_tri_self_transpose(L_Omega_SP);
  }
'
modelAR <- 
  'data { 
    // input information 
    int NO; // number of observations
    int NP; // number of parameters
    int NS; // number of sites
    int NM; // number of months
    // input data 
    int YSP[NO]; // observations (on the original scale)
    int n[NO]; // subsample count (on the original scale)
    int N[NO]; // total count (on the original scale)
    int OBS_MONTH[NO]; // month-site observation index associated with MA
    int OBS_SITE[NO]; // month-site observation index associated with MA
    vector[NO] D; // number of days 
    vector[NO] NB; // number of brushes  
    matrix[NS,NP] MM[NM]; // design matrix
  }
  
  parameters {
    // process variables 
    vector<lower= -10,upper= 10>[NP] beta_SP;
    
    vector<lower=0>[NS] L_sigma_SP;  
    cholesky_factor_corr[NS] L_Omega_SP; // 
      
    vector<lower= -10,upper= 10>[NS] SP[NM];
    vector<lower= -10,upper= 10>[NS] SF[NM]; 
    vector<lower= -0.9, upper= 0.9>[NS] phi;
  }
  
  transformed parameters { 
    vector[NS] LSP[NM]; // log scale expected sp per brush, per day for every month
    
    vector[NO] SP_exp; // normal scale estimated sp for each observation
    vector[NO] SF_exp; // normal scale estimated sp for each observation
    
    vector[NO] N_mu; // expected total count
    vector[NO] q_mu; // expected fraction of N that is SP
    
    vector[NS] lag_res[NM];
    
    matrix[NS,NS] L_Sigma_SP;
    
    L_Sigma_SP <- diag_post_multiply(L_Omega_SP, L_sigma_SP);  
    
    lag_res[1] <- rep_vector(0, NS);
    LSP[1] <- MM[1]*beta_SP;

    for (i in 2:NM){
      lag_res[i] <- phi .*(SP[i-1]-LSP[i-1]);
      LSP[i] <- MM[i]*beta_SP+lag_res[i];
    }

    for(i in 1:NO){
      SP_exp[i] <- exp(SP[OBS_MONTH[i],OBS_SITE[i]])*NB[i]*D[i];
      SF_exp[i] <- exp(SF[OBS_MONTH[i],OBS_SITE[i]])*NB[i]*D[i];
    }
    
    N_mu <- SF_exp+SP_exp;
    q_mu <- SP_exp ./ N_mu;
  }
  
  model { 
    L_Omega_SP ~ lkj_corr_cholesky(NS);
    L_sigma_SP ~ cauchy(0, 2.5);
    
    for (i in 1:NM) {
      SP[i] ~ multi_normal_cholesky(LSP[i], L_Sigma_SP);
    }
    
    YSP ~ binomial(n,q_mu); // define the likelihood
    N ~ poisson(N_mu);
  }
  generated quantities {
    corr_matrix[NS] Omega;
    Omega <- multiply_lower_tri_self_transpose(L_Omega_SP);
  }
'


