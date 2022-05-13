// This Stan program defines a model for GPP as a linear function of light and antecedent precip

data {
    int<lower=0> N;            // number of observations
    int<lower=0> nweights;     // number of time lags
    vector[N] GPP_obs;         // observed GPP time series
    vector<lower=0>[N] light;  // light observations normalized to a max of 1
    matrix<lower=0>[N,nweights] P_ante; // matrix of antecedent driver with a column for each memory lag
    vector<lower=0>[N] missing_dat; // vector with locations of missing data
}

transformed data{
    vector[nweights] w_prior;  // a vector of 1's to serve as w's prior
    for(i in 1:nweights){
        w_prior[i] = 1;
    }
}

parameters {
    real<lower=0> beta_light;     // coefficient on light
    real<lower=0,upper=1> beta_precip; // coefficient on antecedent driver
    simplex[nweights] w;          // weights on each time lag
    real<lower=0> sigma_obs;      // standard deviation on observation error
}

transformed parameters {
    vector<lower=0>[N] ant;
    ant = P_ante * w;
}

model {
    vector[N] GPP; // underlying GPP state

    GPP = (beta_light * light + beta_precip * ant);
    for(i in 1:N){
        if(missing_dat[i] == 1){
            target += normal_lpdf(GPP_obs[i] | GPP[i], sigma_obs);
        }
    }

    // Priors
    beta_light ~ normal(1, 1);
    beta_precip ~ normal(1, 1);
    w ~ dirichlet(w_prior);
    sigma_obs ~ normal(0,0.5);
}

