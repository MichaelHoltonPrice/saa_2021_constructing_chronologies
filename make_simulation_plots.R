# Create simulation plots for Michael Holton Price's following 2021 Society for
# American Archaeology (SAA) presentation:
#
# Title:
#   End-to-end Bayesian inference for summarizing sets of radiocarbon dates
# Session:
#  Constructing Chronologies
# 
# This script relies on the R packages baydem and enneal available at:
#
# https://github.com/eehh-stanford/baydem
# https://github.com/MichaelHoltonPrice/enneal
#
# The specific commits of each used for this script are:
#
# 2ab1cf3495e71d0d509aad7b8d806851186f844f
# 6edd8a53177b221dbd7917597b6c9dfd74e7d612
#
# These versions can be installed in R using the following commands:
#
# devtools::install_github("eehh-stanford/baydem",ref="2ab1cf3495e71d0d509aad7b8d806851186f844f")
# devtools::install_github("MichaelHoltonPrice/enneal",ref="6edd8a53177b221dbd7917597b6c9dfd74e7d612")

# Since multiple cores are used to speed up computation, load some necessary
# packages; also load Bchron (if necessary, install these packages first):
library(doParallel)
library(foreach)
library(parallel)
library(Bchron)
library(baydem)

# Clear the workspace
rm(list=ls(all=TRUE))

# Load the functions in saa_2021_functions.R
source("saa_2021_functions.R")

# Set the model hyperparameters
hp <-
  list(
    # Class of fit (Gaussian mixture)
    fitType = "gaussmix",
    # Parameter for the dirichlet draw of the mixture probabilities
    alpha_d = 1,
    # The gamma distribution shape parameter for sigma
    alpha_s = 3,
    # The gamma distribution rate parameter for sigma, yielding a mode of 300
    alpha_r = (3 - 1) / 300,
    # Minimum calendar date (years BC/AD)
    taumin = 600,
    # Maximum calendar date (years BC/AD)
    taumax = 1300,
    # Spacing for the measurement matrix (years)
    dtau = 1,
    # Number of mixtures
    K = 2
  )

# The simulation distribution: a Gaussian mixture with ordering pi1, pi2, mu1,
# mu2, sig1, sig2
th_sim <-
  c(
    pi1 = 0.2,
    pi2 = 0.8,
    mu1 = 775,
    mu2 = 1000,
    sig1 = 35,
    sig2 = 45
  )

# For reproducibility, set the random number seed (seed from random.org)
set.seed(94385)

# Define the error spec for the fraction modern measurements, which assumes
# measurement errors are uniformily distributed on the interval 0.0021 to
# 0.0028 (the uncertainty is in the fraction modern value, not the radiocarbon
# year).

errorSpec <- list(min=.0021,max=.0028)

# Load the intcal20 calibration curve, which will be used for all results.
calibDf <- baydem::bd_load_calib_curve("intcal20")

# Create 10,000 samples of the calendar date (AD) using the simulation
# distribution
date_AD <- baydem::bd_sample_gauss_mix(10000,th_sim,hp$taumin,hp$taumax)

# Simulate the radiocarbon measurement process
rc_meas <- baydem::bd_draw_rc_meas_using_date(date_AD,calibDf,errorSpec,isAD=T)

# Create the data for plotting

# Create the sampling grid (also used for plotting)
tau <- seq(hp$taumin,hp$taumax,by=hp$dtau)

# Calculate the probability density of the target curve
fsim <- baydem::bd_calc_gauss_mix_pdf(th_sim,
                                      tau,
                                      taumin=hp$taumin,
                                      taumax=hp$taumax)

# Create the histogram data of the sample dates
hist_data <- hist(date_AD,plot=F)

# Calculate the measurement matrix
M <- baydem::bd_calc_meas_matrix(tau,rc_meas$phi_m,rc_meas$sig_m,calibDf,
                                 addCalibUnc=T,useTrapez=F)

# Normalize the rows of the measurement matrix. For clarity, use a for loop
# rather than, say, rowSums.
Mnorm <- M
for(i in 1:nrow(Mnorm)) {
  Mnorm[i,] <- M[i,] / (hp$dtau * sum(Mnorm[i,]))
}

# Create the summed densities from the normalized measurement matrix for 100,
# 1000, and 10000 samples.
fspd100   <- colSums(Mnorm[1:100,])/100
fspd1000  <- colSums(Mnorm[1:1000,])/1000
fspd10000 <- colSums(Mnorm[1:10000,])/10000

# Do maximum likelihood fits for N=100,1000, and 10000 and for K=2 and K=3.
# This can take awhile -- a perfect opportunity to do some stretches or grab a
# cup of coffee.
cores_to_use <- detectCores()-2
max_lik_fit100_K2 <- temper_trunc_gauss_mix(2,
                                            rc_meas$phi_m[1:100],
                                            rc_meas$sig_m[1:100],
                                            hp$taumin,
                                            hp$taumax, hp$dtau,
                                            calibDf,num_restart=10,
                                            num_cores=cores_to_use)
max_lik_fit100_K3 <- temper_trunc_gauss_mix(3,
                                            rc_meas$phi_m[1:100],
                                            rc_meas$sig_m[1:100],
                                            hp$taumin,
                                            hp$taumax, hp$dtau,
                                            calibDf,num_restart=10,
                                            num_cores=cores_to_use)
if(max_lik_fit100_K2$bic < max_lik_fit100_K3$bic) {
  print("For N=100, K=2 is preferred")
} else {
  print("For N=100, K=2 is preferred")
}
max_lik_fit1000_K2 <- temper_trunc_gauss_mix(2,
                                             rc_meas$phi_m[1:1000],
                                             rc_meas$sig_m[1:1000],
                                             hp$taumin,
                                             hp$taumax, hp$dtau,
                                             calibDf,num_restart=10,
                                             num_cores=cores_to_use)
max_lik_fit1000_K3 <- temper_trunc_gauss_mix(3,
                                             rc_meas$phi_m[1:1000],
                                             rc_meas$sig_m[1:1000],
                                             hp$taumin,
                                             hp$taumax, hp$dtau,
                                             calibDf,num_restart=10,
                                             num_cores=cores_to_use)
if(max_lik_fit1000_K2$bic < max_lik_fit1000_K3$bic) {
  print("For N=1000, K=2 is preferred")
} else {
  print("For N=1000, K=2 is preferred")
}
max_lik_fit10000_K2 <- temper_trunc_gauss_mix(2,
                                              rc_meas$phi_m[1:10000],
                                              rc_meas$sig_m[1:10000],
                                              hp$taumin,
                                              hp$taumax, hp$dtau,
                                              calibDf,num_restart=10,
                                              num_cores=cores_to_use)
max_lik_fit10000_K3 <- temper_trunc_gauss_mix(3,
                                              rc_meas$phi_m[1:10000],
                                              rc_meas$sig_m[1:10000],
                                              hp$taumin,
                                              hp$taumax, hp$dtau,
                                              calibDf,num_restart=10,
                                              num_cores=cores_to_use)
if(max_lik_fit10000_K2$bic < max_lik_fit10000_K3$bic) {
  print("For N=10000, K=2 is preferred")
} else {
  print("For N=10000, K=2 is preferred")
}

# Do a modified version of the BchronDensityFast fit for N=100,1000, and 10000
# with K=2 This can also take awhile -- more coffee?
bchron_fit100_K2 <- do_BchronDensityFast_modified_fit(rc_meas$phi_m[1:100],
                                                      rc_meas$sig_m[1:100],
                                                      tau,
                                                      2)

bchron_fit1000_K2 <- do_BchronDensityFast_modified_fit(rc_meas$phi_m[1:1000],
                                                       rc_meas$sig_m[1:1000],
                                                       tau,
                                                       2)

bchron_fit10000_K2 <- do_BchronDensityFast_modified_fit(rc_meas$phi_m[1:10000],
                                                        rc_meas$sig_m[1:10000],
                                                        tau,
                                                        2)

# Do KDE fits
#kde_fit100 <- fit_kde(rc_meas$phi_m[1:100],
#                      rc_meas$sig_m[1:100],
#                      tau,
#                      num_samp_g=10,
#                      num_samp_t=100)
#
#kde_fit1000 <- fit_kde(rc_meas$phi_m[1:1000],
#                       rc_meas$sig_m[1:1000],
#                       tau,
#                       num_samp_g=10,
#                       num_samp_t=100)
#
#kde_fit10000 <- fit_kde(rc_meas$phi_m[1:10000],
#                        rc_meas$sig_m[1:10000],
#                        tau,
#                        num_samp_g=10,
#                        num_samp_t=100)

# Determine the maximum density across relevant curves so that all plots can be
# set to the same y-range.
fmax <- max(fsim,hist_data$density,fspd100,fspd1000,fspd10000,
            max_lik_fit100_K2$f,max_lik_fit1000_K2$f,max_lik_fit10000_K2$f,
            bchron_fit100_K2$f,bchron_fit1000_K2$f,bchron_fit10000_K2$f,
            kde_fit100$kde,kde_fit1000$kde)

# Show just the target curve
pdf('sim_target.pdf',width=8,height=6)
  plot(tau,fsim,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="blue")
  legend("topright",
         legend=c("Target"),
         col="blue",lty=1)
dev.off()

# Show the target curve and a histogram of the samples
pdf('sim_target_with_hist.pdf',width=8,height=6)
  hist(date_AD,freq=F,
       xlab='Calendar Year [AD]',ylab='Density',main="",
       xlim=c(hp$taumin,hp$taumax),ylim=c(0,fmax))
  lines(tau,fsim,,lwd=3,col="blue")
  legend("topright",
         legend=c("Target"),
         col="blue",
         lty=1)
dev.off()

# Show the SPD for N=10000
pdf('sim_spd_10000.pdf',width=8,height=6)
  plot(tau,fspd10000,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="black")
  legend("topright",
       legend="SPD",
       col="black",
       lty=1)
dev.off()

# Show the target curve and the SPD for N=10000
pdf('sim_target_spd_10000.pdf',width=8,height=6)
  plot(tau,fspd10000,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="black")
  lines(tau,fsim,,lwd=3,col="blue")
  legend("topright",
       legend=c("Target","SPD"),
       col=c("blue","black"),
       lty=c(1,1))
dev.off()

# Compares fits for N=100
pdf('sim_target_spd_max-lik_bchron_100.pdf',width=8,height=6)
  plot(tau,fspd100,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="black")
  lines(tau,fsim,,lwd=3,col="blue")
  lines(max_lik_fit100_K2$tau,max_lik_fit100_K2$f,lwd=3,col="red")
  lines(bchron_fit100_K2$tau,bchron_fit100_K2$f,lwd=3,col="black",lty=3)
  legend("topright",
       legend=c("Target","SPD","Max-Lik","Bchron"),
       col=c("blue","black","red","black"),
       lty=c(1,1,1,3))
dev.off()

#pdf('sim_target_spd_max-lik_kde_100.pdf',width=8,height=6)
#  plot(tau,fspd100,
#       xlab='Calendar Year [AD]',ylab='Density',
#       ylim=c(0,fmax),type="l",lwd=3,col="black")
#  lines(tau,fsim,,lwd=3,col="blue")
#  lines(max_lik_fit100_K2$tau,max_lik_fit100_K2$f,lwd=3,col="red")
#  lines(kde_fit100$tau,kde_fit100$f,lwd=3,col="black",lty=3)
#  legend("topright",
#         legend=c("Target","SPD","Max-Lik","KDE"),
#         col=c("blue","black","red","black"),
#         lty=c(1,1,1,3))
#dev.off()

# Compares fits for N=1000
pdf('sim_target_spd_max-lik_bchron_1000.pdf',width=8,height=6)
  plot(tau,fspd1000,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="black")
  lines(tau,fsim,,lwd=3,col="blue")
  lines(max_lik_fit1000_K2$tau,max_lik_fit1000_K2$f,lwd=3,col="red")
  lines(bchron_fit1000_K2$tau,bchron_fit1000_K2$f,lwd=3,col="black",lty=3)
  legend("topright",
       legend=c("Target","SPD","Max-Lik","Bchron"),
       col=c("blue","black","red","black"),
       lty=c(1,1,1,3))
dev.off()

#pdf('sim_target_spd_max-lik_kde_1000.pdf',width=8,height=6)
#  plot(tau,fspd1000,
#       xlab='Calendar Year [AD]',ylab='Density',
#       ylim=c(0,fmax),type="l",lwd=3,col="black")
#  lines(tau,fsim,,lwd=3,col="blue")
#  lines(max_lik_fit1000_K2$tau,max_lik_fit1000_K2$f,lwd=3,col="red")
#  lines(kde_fit1000$tau,kde_fit1000$f,lwd=3,col="black",lty=3)
#  legend("topright",
#         legend=c("Target","SPD","Max-Lik","KDE"),
#         col=c("blue","black","red","black"),
#         lty=c(1,1,1,3))
#dev.off()

# Compares fits for N=10000
pdf('sim_target_spd_max-lik_bchron_10000.pdf',width=8,height=6)
  plot(tau,fspd10000,
       xlab='Calendar Year [AD]',ylab='Density',
       ylim=c(0,fmax),type="l",lwd=3,col="black")
  lines(tau,fsim,,lwd=3,col="blue")
  lines(max_lik_fit10000_K2$tau,max_lik_fit10000_K2$f,lwd=3,col="red")
  lines(bchron_fit10000_K2$tau,bchron_fit10000_K2$f,lwd=3,col="black",lty=3)
  legend("topright",
       legend=c("Target","SPD","Max-Lik","Bchron"),
       col=c("blue","black","red","black"),
       lty=c(1,1,1,3))
dev.off()
