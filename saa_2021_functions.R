# Do a maximum likelihood fit for a truncated Gaussian mixture model using
# parallel tempering. This function is a candidate for inclustion in the R
# package baydem.
temper_trunc_gauss_mix <- function(K,
                                   phi_m,
                                   sig_m,
                                   tau_min,
                                   tau_max,
                                   dtau,
                                   calib_df, num_restart=10,
                                   num_cores=NA) {

  # Create the sampling grid, tau
  tau <- seq(tau_min,tau_max,dtau)

  # Calculate the measurement matrix, M
  M <- bd_calc_meas_matrix(tau, phi_m, sig_m, calib_df)

  # Iterate over restarts. To generate the starting vectors, assume the
  # following:
  # (1) pi, the mixture proportions, are drawn from a Dirichlet distribution.
  # (2) mu, the means, are restricted to the interval tau_min to tau_max.
  # (3) s, the standard deviations, are restricted to the interval 10*dtau to
  #     tau_max-tau_min

  # Create the lower and upper bounds for optimization
  lower_bounds <- c(rep(0,K-1),rep(tau_min,K),rep(dtau*10,K))
  upper_bounds <- c(rep(1,K-1),rep(tau_max,K),rep(tau_max-tau_min,K))

  # Re-wrap the negative log-likelihood (the objective function) to enforce
  # the bounds on the parameter vector.
  neg_log_lik_rewrap <- function(th_reduced,M,tau,
                                 lower_bounds,upper_bounds) {
    if (any(th_reduced < lower_bounds)) {
      return(Inf)
    }

    if (any(th_reduced > upper_bounds)) {
      return(Inf)
    }

    return(calc_trunc_gauss_mix_neg_log_lik(th_reduced,M,tau))
  }

  # For the parallel tempering, use a temperature vector with temperatures
  # between sqrt(10) and 1/sqrt(10)
  temp_vect <- rev(10^seq(-0.5,0.5,len=10))

  # For the largest temperature, use the following for the standard deviation of
  # the proposal distribution:
  # pi_k  0.5
  # mu_k  (tau_max-tau_min)/2
  #  s_k  (tau_max-tau_min)/4
  tau_range <- tau_max - tau_min
  scale_vect0 <- c(rep(0.5,K-1),rep(tau_range/2,K),rep(tau_range/4,K))

  # For other temperatures, rescale the standard deviation used for the proposal
  # distribution such that the the smallest temperature has a scale 0.01 that
  # of the largest.
  scale_matrix <- matrix(NA,3*K-1,length(temp_vect))
  rescale_vect <- seq(1,.01,len=length(temp_vect))
  for(tt in 1:length(temp_vect)) {
    scale_matrix[,tt] <- scale_vect0 * rescale_vect[tt]
  }

  # Create a matrix of starting parameter vectors for the restart
  TH0_reduced <- matrix(NA,num_restart,3*K-1)
  for(n in 1:num_restart) {
    th0 <- c(MCMCpack::rdirichlet(1,rep(1,K)),
             sort(runif(K,tau_min,tau_max)),
             runif(K,10*dtau,tau_max-tau_min))
    # For optimization, the parameter vector must be reduced in length by 1
    # since the mixture proportions must sum to 1. Therefore, remove the first
    # entry of th0
    TH0_reduced[n,] <- th0[-1]
  }

  if (is.na(num_cores)) {
    # If num_cores is NA, iterate over restarts using a regular for loop
    # A list of fits
    fit_list <- list()
    for (n in 1:num_restart) {
     # Do the actual fit for this restart
     fit_list[[n]] <- enneal::par_temper(TH0_reduced[n,],
                                         #calc_trunc_gauss_mix_neg_log_lik,
                                         neg_log_lik_rewrap,
                                         samps_per_cyc=10,
                                         temp_vect=temp_vect,
                                         prop_scale=scale_matrix,
                                         num_cyc=1000,
                                         M=M, tau=tau,
                                         lower_bounds=lower_bounds,
                                         upper_bounds=upper_bounds)
    }
  } else {
    # For speed, do not use :: to call functions. This means that doParallel
    # (for registerDoParallel) and foreach (for %dopar%) must have been
    # pre-loaded.
    registerDoParallel(num_cores)
    # Must explicitly export functions so that dopar works on Windows
    function_exports <- c("calc_trunc_gauss_mix_neg_log_lik",
                          "is_th_reduced_valid")
    fit_list <-
      foreach(n=1:num_restart,.packages=c('baydem'),
              .export=function_exports) %dopar% {
        enneal::par_temper(TH0_reduced[n,],
                           neg_log_lik_rewrap,
                           samps_per_cyc=10,
                           temp_vect=temp_vect,
                           prop_scale=scale_matrix,
                           num_cyc=1000,
                           M=M, tau=tau,
                           lower_bounds=lower_bounds,
                           upper_bounds=upper_bounds)
      }
  }

  # Identify the restart with the best value of the objective function
  best_chain_in_fit <- function(fit) {
    ind_min <- which.min(unlist(lapply(fit$chains,
                                       function(chain){chain$eta_best})))
    best_chain <- fit$chains[[ind_min]]
  }

  best_chains <- lapply(fit_list,function(fit){best_chain_in_fit(fit)})

  ind_best <- which.min(unlist(lapply(best_chains,function(fit){fit$eta_best})))
  th_best  <- best_chains[[ind_best]]$theta_best
  th_best <- c(1-sum(th_best[1:(K-1)]),th_best)
  eta_best <- best_chains[[ind_best]]$eta_best
  #ind_best <- which.min(unlist(lapply(best_chain_in_fit(fit_list),function(fit){fit$eta_best})))
#  th_best  <- fit_list[[ind_best]]$theta_best
#  th_best <- c(1-sum(th_best[1:(K-1)]),th_best)
#  eta_best <- fit_list[[ind_best]]$eta_best

  # Calculate the probability density for the best fit parameter
  f  <- bd_calc_gauss_mix_pdf(th_best,tau,tau_min,tau_max)

  # Calculate the Bayesian information criterion (BIC) and Akaike information
  # criterion (AIC). General wisdom holds that BIC is preferred if the true
  # model is in the set of models tested and AIC is preferred if it may not be.
  N <- length(phi_m)
  bic <- log(N)*(3*K-1) + 2*eta_best
  aic <-      2*(3*K-1) + 2*eta_best

  # Return a list with:
  # th          The best fit paramater vector
  # neg_log_lik The best value of the objective function
  # tau         The sampling grid
  # f           The probability density
  return(list(th=th_best,neg_log_lik=eta_best,tau=tau,f=f,bic=bic,aic=aic))
}

# This function is a candidate for inclustion in the R package baydem.
calc_trunc_gauss_mix_neg_log_lik <- function(th_reduced,M,tau,sig_min=0) {
  # If the parameter vector is invalid, return infinity
  if (!is_th_reduced_valid(th_reduced,sig_min)) {
    return(Inf)
  }

  tau_min <- min(tau)
  tau_max <- max(tau)

  # Add an undersore to pi, the mixture proportions, since pi is 3.14... in R
  K <- (length(th_reduced) + 1)/3
  pi_ <- th_reduced[1:(K-1)]

  th <- c(1-sum(pi_),th_reduced)
  # Calculate v, the vector of densities
  v <- bd_calc_gauss_mix_pdf(th,tau,tau_min,tau_max)
  h <- M %*% v
  return(-sum(log(h)))
}

# This function is a candidate for inclustion in the R package baydem.
is_th_reduced_valid <- function(th_reduced,sig_min=0) {
  # TODO: Consider adding a check on the length of th_reduced
  K <- (length(th_reduced) + 1)/3

  # Add an undersore to pi, the mixture proportions, since pi is 3.14149... in R
  pi_ <- th_reduced[1:(K-1)]
  pi_ <- c(sum(pi_),pi_) # mixture proportions must sum to 1

  # Mixture propotions must all lie between 0 and 1
  if (any(pi_ < 0)) {
    return(FALSE)
  }

  if (any(pi_ > 1)) {
    return(FALSE)
  }

  # The means must be ordered
  mu <- th_reduced[K:(2*K-1)]
  if (is.unsorted(mu)) {
    return(FALSE)
  }

  # The standard deviations must be strictly positive (or greater than sig_min)
  sig <- th_reduced[(2*K):(3*K-1)]
  if(any(sig <= sig_min)) {
    return(FALSE)
  }
  return(TRUE)
}

# The following function is based on:
# https://github.com/andrewcparnell/Bchron/blob/master/R/BchronDensityFast.R
#
# The only difference is that in the call to mclust::densityMclust a one
# dimensional variable/unequal variance model is used rather than letting the
# model be chosen via the BIC. This ensures that this calculation most closely
# mirrors the inference done with baydem.
BchronDensityFast_modified <-
  function(ages, ageSds, calCurves, pathToCalCurves = system.file("data", package = "Bchron"), dfs = rep(100, length(ages)), samples = 2000, G = 30) {
    if (length(ages) != length(ageSds)) stop("ages and 1-sigma errors must be same length")
    if (length(ages) != length(calCurves)) stop("ages and Calibration curves must be same length")

    # Calibrate ages
    x <- BchronCalibrate(ages = ages, ageSds = ageSds, calCurves = calCurves, pathToCalCurves = pathToCalCurves, dfs = rep(100, length(ages)))

    # Get number of dates
    n <- length(x)

    # Get a huge load of samples from the posteriors here
    thetaBig <- vector(length = n * samples)
    for (i in 1:n) thetaBig[((i - 1) * samples + 1):(i * samples)] <- sample(x[[i]]$ageGrid, size = samples, prob = x[[i]]$densities, replace = TRUE)

    # Now run mclust
    mclustOutput = mclust::densityMclust(data = thetaBig,
                                     G = G,
                                     modelNames = 'V')

    output <- list(out = mclustOutput, calAges = x)
    class(output) <- "BchronDensityRunFast"
    return(output)
  }

do_BchronDensityFast_modified_fit <- function(phi_m,sig_m,tau,K) {
  kappa <- 8033
  trc_m <- -kappa*log(phi_m)
  sig_trc_m <- kappa * sig_m / phi_m

  N <- length(phi_m)
  bchronFit <- BchronDensityFast_modified(ages=round(trc_m),
                                          ageSds=round(sig_trc_m),
                                          calCurves=rep('intcal20',N),G=K)
  f <- rep(0,length(tau))
  for(k in 1:K) {
    f <-
      f + bchronFit[[1]]$parameters$pro[k]*
        dnorm(1950-tau,bchronFit[[1]]$parameters$mean[k],
              sqrt(bchronFit[[1]]$parameters$variance$sigmasq[k]))
  }

  return(list(tau=tau,f=f,bchronFit=bchronFit))
}

# Fit a KDE to the radiocarbon measurements ala Bronk Ramsey (2009).
# Unfortunately, Bronk Ramsey (2009) does not describe the algorithm in
# sufficient detail to be certain of the precise algorithm used and I do not
# have access to the Oxcal code. Therefore, I supplement my reading of Bronk
# Ramsey (2009) with the information in Carleton and Groucutt (2020). This
# yields myy following best guess for how the Oxcal KDE algorithm works.
#
# Bronk Ramsey 2009 -- Methods for summarizing 14C dates
# Carleton and Groucutt 2020 -- Sum things are not what they seem
#
# The algorithm involves doing the following multiple times:
#
# (1) Take a sample of potential dates (in Oxcal, stratigraphic information can
#     can be used at this step, but for my implementation I simply use
#     bchron::BchronCalibrate to make samples)
# (2) Sample the bandwidth modifier, g, which is uniformily distributed on 0 to
#     1.
# (3) Utilize a Metropolis accept-reject step where the predictive likelihood
#     (Equation 5 in Bronk Ramsey (2009)) is used to weight the acceptance
#     probabilities
# This function is a candidate for inclustion in the R package baydem.
fit_kde <- function(phi_m,sig_m,tau,num_samp_t,num_samp_g) {
  kappa <- 8033
  trc_m <- -kappa*log(phi_m)
  sig_trc_m <- kappa * sig_m / phi_m
  N <- length(phi_m)
  calib <- Bchron::BchronCalibrate(ages=round(trc_m),
                                   ageSds=round(sig_trc_m),
                                   calCurves=rep('intcal20',N))

  # Sample the calibrated curves. The variable t_e_matrix has dimensions
  # num_samp_t by N. The samples are in BP (before present).
  t_e_matrix <- Bchron::sampleAges(calib,num_samp_t)

  kde <- rep(0,length(tau))
  tau_BP <- rev(1950 - tau)
  for (s in 1:num_samp_t) {
    # Extract the event times
    t_e <- as.vector(t_e_matrix[s,])
    # Calculate the Silverman bandwidth
    h_S <- (4/3)^(1/5) * sd(t_e) * N^(-1/5)

    for (i in 1:num_samp_g) {
      # Draw the proposal for g, the bandwidth modifier, from the interval 0 to 1
      g_prop <- runif(1)

      # Calculate the bandwidth
      h <- g_prop * h_S

      # Vectorize the kernel calculation to speed up computation. This requires
      # attention to detail since (adopting Bronk Ramsey's notation) j=i is
      # excluded from the sum. To begin, create two matrices:
      #
      # (a) A matrix in which the columns are t_e
      # (b) A matrix in which the rows are t_e
      A <- replicate(N,t_e)
      B <- t(A)

      # Calculate the kernel
      K <- dnorm((B-A)/h)

      # The rows now have the values that need to be summed, except that we must
      # substract K(0) from each row to account for the fact that j cannot equal
      # i.
      neg_log_lik_prop <-
        -(N-2)*log(1/(N-1)/h) -(N-2)/N*sum(log(rowSums(K) - dnorm(0)))

      if(i == 1) {
          # For the first sample, we must accept the "proposal"
          accept <- 1
      } else {
        # Do a Metroplis step
        if (is.na(neg_log_lik_prop)) {
          accept <- FALSE
        } else {
          if (!is.finite(neg_log_lik_prop)) {
            accept <- FALSE
          } else {
            alpha <- exp(-(neg_log_lik_prop-neg_log_lik))
            if(alpha > 1) {
              alpha <- 1
            }
            accept <- runif(1) <= alpha
          }
        }
      }
      if(accept) {
        g <- g_prop
        neg_log_lik <- neg_log_lik_prop
      }
      # Update the density --  vectorized as above, with:
      # (a) A matrix in which the columns are t_e
      # (b) A matrix in which the rows are tau
      A <- replicate(length(tau_BP),t_e)
      B <- t(replicate(length(t_e),tau_BP))

      kde_si <- colSums(dnorm((B-A)/h))/h/num_samp_g/N
      kde <- kde + kde_si/num_samp_t
    }
  }
  return(list(kde=rev(kde),
              tau=tau,
              t_e_matrix=t_e_matrix))
}
