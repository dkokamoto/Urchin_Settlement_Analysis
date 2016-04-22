
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
    vector<lower= -10,upper= 10>[NP] beta_SF;
    
    vector<lower=0>[NS] L_sigma_SP;  
    cholesky_factor_corr[NS] L_Omega_SP; // 
    vector<lower=0>[NS] L_sigma_SF;  
    cholesky_factor_corr[NS] L_Omega_SF; // 
      
    vector<lower= -10,upper= 10>[NS] SP[NM];
    vector<lower= -10,upper= 10>[NS] SF[NM]; 
    vector<lower= -0.9, upper= 0.9>[NS] phi_SP;
    vector<lower= -0.9, upper= 0.9>[NS] phi_SF;
  }
  
  transformed parameters { 
    vector[NS] LSP[NM]; // log scale expected sp per brush, per day for every month
    vector[NS] LSF[NM]; // log scale expected sf per brush, per day for every month
    
    vector[NO] SP_exp; // normal scale estimated sp for each observation
    vector[NO] SF_exp; // normal scale estimated sf for each observation
    
    vector[NO] N_mu; // expected total count
    vector[NO] q_mu; // expected fraction of N that is SP
    
    vector[NS] lag_res_SP[NM];
    vector[NS] lag_res_SF[NM];
    
    matrix[NS,NS] L_Sigma_SP;
    matrix[NS,NS] L_Sigma_SF;
    
    L_Sigma_SP <- diag_pre_multiply(L_sigma_SP,L_Omega_SP);  
    L_Sigma_SF <- diag_pre_multiply(L_sigma_SF,L_Omega_SF);  
    
    lag_res_SP[1] <- rep_vector(0, NS);
    lag_res_SF[1] <- rep_vector(0, NS);
    LSP[1] <- MM[1]*beta_SP;
    LSF[1] <- MM[1]*beta_SF;
    
    for (i in 2:NM){
      lag_res_SP[i] <- phi_SP .*(SP[i-1]-LSP[i-1]);
      LSP[i] <- MM[i]*beta_SP+lag_res_SP[i];
      lag_res_SF[i] <- phi_SF .*(SF[i-1]-LSF[i-1]);
      LSF[i] <- MM[i]*beta_SF+lag_res_SF[i];
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
    L_Omega_SF ~ lkj_corr_cholesky(NS);
    L_sigma_SF ~ cauchy(0, 2.5);
    
    for (i in 1:NM) {
      SP[i] ~ multi_normal_cholesky(LSP[i], L_Sigma_SP);
      SF[i] ~ multi_normal_cholesky(LSF[i], L_Sigma_SF);
    }
    
    YSP ~ binomial(n,q_mu); // define the likelihood
    N ~ poisson(N_mu);
  }
  generated quantities {
    corr_matrix[NS] Omega_SP;
    corr_matrix[NS] Omega_SF;
    Omega_SP <- multiply_lower_tri_self_transpose(L_Omega_SP);
    Omega_SF <- multiply_lower_tri_self_transpose(L_Omega_SF);
  }
'


