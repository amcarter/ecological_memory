# Functions for simulating metabolism data
# A Carter
# 4/2022

# Antecedent ####
calc_antecedent_drivers <- function(driver, nweights, scale = 1){
    driver <- zoo::na.approx(driver, na.rm = F)
    if(sum(is.na(driver)) > 0) {
        print('the antecedent driver variable should not start or end with NA')
        return()
    }

    # calc antecedent interval lengths
    index <- seq(1,nweights, by = 1)
    intervals = scale ^ index

    ant <- data.frame(matrix(0, ncol = nweights, nrow = length(driver)))
    colnames(ant) <- paste0('int_', seq(1:nweights))

    N <- 0
    for (i in 1:nweights){
        start <- N + 1
        N <- N + intervals[i]
        for(j in start:N){
            pd <- c(rep(NA, j), driver[1:(length(driver)-j)])
            ant[, i] <- ant[, i] + pd
        }
        ant[, i] <- ant[, i]/intervals[i]
        ant[1:N, i] <- 0
    }

    return(ant)
}


calc_antecedent_driver_vi <- function(GPP, nweights, w = NULL,
                                   intervals = NULL, scale = NULL){

    GPP <- zoo::na.approx(GPP, na.rm = F)
    if(sum(is.na(GPP)) > 0) {
        print('GPP should not start or end with NA')
        return()
    }

    if(is.null(w)) w = rep(1/nweights, nweights)
    if(is.null(intervals)){
        index <- seq(0,(nweights-1), by = 1)
        intervals = scale * 2 ^ index
    }

    n_pre <- sum(intervals)
    Pant <- GPP <- c(rep(GPP[1], n_pre), GPP)
    ndays <- length(GPP)

    for (i in (n_pre+1):ndays){
        Pvec <- numeric(nweights)
        N <- 0
        for(j in 1:nweights){
            n <- intervals[j]
            pp <- numeric(n)
            for(k in 1:n){
                pp[k] <- GPP[i-N-k]
            }
            Pvec[j] <- w[j] * mean(pp)
            N <- N + n
        }
        Pant[i]<-sum(Pvec)
    }

    Pant <- Pant[(n_pre+1):ndays]

    return(Pant)
}

# Simulate detrital carbon and detrital based respiration ####
calc_rate_coef <- function(temp_C,
                           K_20 = 0.01,# respiration rate coefficient at 20C
                           E_a = 0.63 # activation energy for heterotrophic respiration (eV)
                           ){
    k_b = 8.6173*10^-5 # boltzmann's constant in eV/K
    K = K_20 * exp(-E_a/k_b *(1/(temp_C+273) - 1/293))

    return(K)
}

# Simulation functions for detrital C in rivers: ####

# Shear stress and storm sloughing ####
calc_tau_gm2 <- function(depth, slope){
    tau = 10^3 * depth * slope

    return(tau)
}

calc_disturbance <- function(tau, tau0){
    tau_max = max(tau, na.rm = T)
    x2 = case_when(tau < tau0 ~ 0,
                   TRUE ~ 1/(tau_max - tau0))
    disturb <- ((tau - tau0)*x2)^2

    return(disturb)
}

# litterfall ####
calc_litter_oneyear <- function(annual_litter, LAI){
    LAI <- zoo::na.approx(LAI, na.rm = F)
    diffLAI <- c(0, diff(LAI))
    diffLAI <- case_when(diffLAI < 0 ~ diffLAI,
                         TRUE ~ 0)
    litter <- annual_litter * diffLAI/sum(diffLAI, na.rm = T)

    return(litter)
}

calc_litter_from_LAI <- function(dd){
    if(!all(c('date', 'LAI') %in% colnames(dd))){
        print('Your dataframe must have columns called \'date\' and \'LAI\'')
        return()
    }
    if(!'year' %in% colnames(dd)){
        dd$year = lubridate::year(dd$date)
    }
    tmp <- data.frame()
    for(y in unique(dd$year)){
        yy <- dd[dd$year == y,]
        yy$litter = calc_litter_oneyear(100*max(yy$LAI), yy$LAI)
        tmp <- bind_rows(tmp, yy)
    }

    return(tmp$litter)
}


# Stan Wrappers ####
fill_stan_dat <- function(stan_dat){

    for(i in 1:length(stan_dat)){
        if(is.data.frame(stan_dat[[i]])) next
        stan_dat[[i]] <- zoo::na.approx(stan_dat[[i]], na.rm = F)
        if(sum(is.na(stan_dat[[i]])) != 0){
            print('Your data cannot start or end with NA')
            return()
        }
    }

    return(stan_dat)
}

plot_post_sim <- function(fit, pars, vals, title = NULL, xlim = NULL){
    dd <- data.frame(x = vals, y = length(vals):1)
    p <- rstan::plot(fit, show_density = T, fill_color = 'grey',
                     pars = pars) +
        geom_point(data = dd, aes(x = x, y = y), size = 3,
                   shape = 17, col = 'brown3')+
        xlim(xlim)+
        ggtitle(title)
    return(p)
}
