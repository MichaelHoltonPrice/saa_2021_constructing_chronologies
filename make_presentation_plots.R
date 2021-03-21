uncalYear <- 1100 # Years BP
uncalSig  <- 20   # 20 year measurement uncertainty

phi_m <- exp(-uncalYear/ 8033)
sig_m <- uncalSig * exp(-uncalYear / 8033) / 8033

tau_min <- 850
tau_max <- 1050
dtau    <- 5
tau <- seq(tau_min,tau_max,by=dtau)
calibDf = baydem::bd_load_calib_curve("intcal13")
M <- baydem::bd_calc_meas_matrix(tau, phi_m, sig_m, calibDf, T, F)

fprior <- rep(1/(tau_max-tau_min),length(tau))
fpost  <- M[1,] / (sum(M[1,]*dtau))
pdf('single_obs_inf_plot1.pdf',width=8,height=6)
  plot(tau,fprior,type='l',lwd=3,ylim=c(0,max(fprior,fpost)),xlab='Calendar Year [AD]',ylab='Density')
dev.off()

pdf('single_obs_inf_plot2.pdf',width=8,height=6)
  baydem::bd_vis_calib_curve(tau_min, tau_max, calibDf, xlab = "Calendar Year [AD]", ylab = "Fraction Modern", invertCol = "gray80")
dev.off()

pdf('single_obs_inf_plot3.pdf',width=8,height=6)
  plot(tau,fprior,type='l',lwd=3,ylim=c(0,max(fprior,fpost)),xlab='Calendar Year [AD]',ylab='Density',col='grey')
  lines(tau,fpost,lwd=3)
dev.off()

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

simData <- readRDS('sim1000.rds')

# Minimum calendar date (AD)
taumin = 600
# Maximum calendar date (AD)
taumax = 1300
# Spacing for the measurement matrix (years)
dtau = 1

tau <- seq(taumin,taumax,by=dtau)
M <- baydem::bd_calc_meas_matrix(tau, simData$phi_m, simData$sig_m, calibDf, addCalibUnc=T, useTrapez=F)

# Normalize densities. For clarity, use a for loop rather than, say, rowSums
for(i in 1:nrow(M)) {
  M[i,] <- M[i,] / (dtau * sum(M[i,]))
}

# Check summed density against rcarbon
mycaldates <- rcarbon::calibrate(simData$trc_m,simData$sig_trc_m, normalised=T,F14C=T)
myspd <- rcarbon::spd(mycaldates, timeRange=c(1350,650))

if(F) { # Set T to check individual calibrations
  for(n in 1:nrow(M)) {
    print(n)
    tau1 <- 1950-mycaldates$grids[[n]]$calBP
    f1   <- mycaldates$grids[[n]]$PrDens
    f1   <- f1/sum(f1) # assumes spacing of tau1 is 1 year
    f2 <- M[n,]
    fmax <- max(f1,f2)
    plot(tau1,f1,type='l',col='green',lwd=3,ylim=c(0,fmax))
    lines(tau,M[n,],lwd=3,col='black')
    readline(prompt="Press [enter] to continue")
  }
}

f1 <- M[1,]
f2 <- M[2,]
f3 <- M[3,]
fall <- colMeans(M)

fmax <- max(f1,f2,f3)

pdf('spd1.pdf',width=8,height=6)
  plot(tau,f1,type='l',lwd=3,ylim=c(0,fmax),xlab='Calendar Year [AD]',ylab='Density')
dev.off()

pdf('spd2.pdf',width=8,height=6)
  plot(tau,f1,type='l',lwd=3,ylim=c(0,fmax),xlab='Calendar Year [AD]',ylab='Density',col='grey')
  lines(tau,f2,lwd=3)
dev.off()

pdf('spd3.pdf',width=8,height=6)
  plot(tau,f1,type='l',lwd=3,ylim=c(0,fmax),xlab='Calendar Year [AD]',ylab='Density',col='grey')
  lines(tau,f2,lwd=3,col='grey')
  lines(tau,f3,lwd=3)
dev.off()

fsim <- th_sim['pi1']*dnorm(tau,th_sim['mu1'],th_sim['sig1']) + th_sim['pi2']*dnorm(tau,th_sim['mu2'],th_sim['sig2'])
fmax <- max(fall,fsim)
pdf('spdall.pdf',width=8,height=6)
  plot(tau,fall,type='l',lwd=3,xlab='Calendar Year [AD]',ylab='Density',ylim=c(0,fmax))
dev.off()


pdf('spdall_sim.pdf',width=8,height=6)
  plot(tau,fall,type='l',lwd=3,xlab='Calendar Year [AD]',ylab='Density',ylim=c(0,fmax))
  lines(tau,fsim,lwd=3,col='blue')
dev.off()

pdf('spdcheck.pdf',width=8,height=6)
  tau_rcarb <- 1950 - myspd$grid[,1]
  dtau_rcarb <- unique(diff(tau_rcarb))
  if(length(dtau_rcarb) != 1) {
    stop('myspd is not regularly spaced')
  }
  f_rcarb <- myspd$grid[,2]
  f_rcarb <- f_rcarb / sum(f_rcarb) / dtau_rcarb
  fmax <- max(f_rcarb,fall,fsim)
  plot(tau_rcarb,f_rcarb,col='green',lwd=3,type='l',xlab='Calendar Date [AD]',ylab='Density',ylim=c(0,fmax))
  lines(tau,fall,lwd=3,col='black')
  lines(tau,fsim,lwd=3,col='blue')
dev.off()


