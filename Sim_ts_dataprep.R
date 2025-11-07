library("purrr")
library("spam")
library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library("plotrix")
library("nloptr")
library("parallel")
library("sgt")
source("Rfns.R")


# =======================================================
# SIMULATION SET-UP
# Synthetic dataset generated from an AR(1)-GARCH(1,1) process with skewed t innovations
# AR-GARCH filter parameters:
# mu=-.05; ar1 = .3 # AR(1) part
# omega=.01; al=.1; be=.85 # GARCH(1,1) parameters
# Innovation distribution parameters:
# nu=5 # shape parameter
# ga=1.5 # skewness parameter
# Burn-in period of 1000 points was used
# The simulated series is saved as "simdat.RDATA"

n=5000 # out-of-sample size to evaluate forecasts

w=500 # moving window size


set.seed(2345)
spec = garchSpec(model = list(mu = -.05, ar = .3, omega = .01, alpha = .1, beta = .85, skew = 1.5, shape = 5),
                 cond.dist = "sstd")
simdat = garchSim(spec, n = n+2*w, n.start = 1000, extended = TRUE)
save(simdat, file = "simdat.RDATA")
# =======================================================

load("simdat.RDATA")
mut = (simdat$garch) - (simdat$sigma) * (simdat$eps)
sigt = simdat$sigma
simdat = simdat$garch

avec=c(.90, .95, .99) # vector of alpha levels for VaR
tvec=c(.96561, .98761, .99855) # vector of tau levels for expectile
nvec=c(.754, .875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

x=tail(simdat,n+2*w) # time series to be used for fitting and forecasting


# =======================================================
# Normal innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1)
# esmodel = NULL # expected shortfall for the assumed model/distribution
# for(i in inu)
#   esmodel = c(esmodel, integrate(function(x) x*dnorm(x), qmodel[i], Inf)$value/(1-VaR.levels[i]))
esmodel = (ddist("norm", qdist("norm",p=nvec,mu=0,sigma=1), mu=0, sigma=1)) / (1 - nvec)
emodel = enorm(asy=tvec)

VaR <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRfhs <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRevt <- matrix(nrow=n+w,ncol=length(VaR.levels))

ES <- matrix(nrow=n+w,ncol=length(nvec))
ESfhs <- matrix(nrow=n+w,ncol=length(nvec))
ESevt <- matrix(nrow=n+w,ncol=length(nvec))

EXP <- matrix(nrow=n+w,ncol=length(tvec))
EXPfhs <- matrix(nrow=n+w,ncol=length(tvec))
EXPevt <- matrix(nrow=n+w,ncol=length(tvec))


# estimated parameters
xi <- vector(mode="numeric", length=n+w)
beta  <- vector(mode="numeric", length=n+w)
fit.par <- matrix(nrow=n+w, ncol=5) # 5 model parameters
mut  <- vector(mode="numeric", length=n+w)
sigt  <- vector(mode="numeric", length=n+w)

# ----------------------------------------------------

for(i in 1:(n+w))
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
  
  VaR[i,]=foc$VaRmodel
  VaRfhs[i,] = foc$VaRfhs
  VaRevt[i,] = foc$VaRevt
  
  EXP[i,]=foc$EXPmodel
  EXPfhs[i,] = foc$EXPfhs
  EXPevt[i,] = foc$EXPevt
  
  ES[i,]=foc$ESmodel
  ESfhs[i,] = foc$ESfhs
  ESevt[i,] = foc$ESevt
  xi[i] = foc$evt.shape
  beta[i] = foc$evt.scale
  mut[i] = foc$mut
  sigt[i] = foc$sigt
  
  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)
save(out, file="Sim3norm3.RDATA")

# =======================================================
# Student t innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")

VaR <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRfhs <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRevt <- matrix(nrow=n+w,ncol=length(VaR.levels))

ES <- matrix(nrow=n+w,ncol=length(nvec))
ESfhs <- matrix(nrow=n+w,ncol=length(nvec))
ESevt <- matrix(nrow=n+w,ncol=length(nvec))

EXP <- matrix(nrow=n+w,ncol=length(tvec))
EXPfhs <- matrix(nrow=n+w,ncol=length(tvec))
EXPevt <- matrix(nrow=n+w,ncol=length(tvec))

xi <- vector(mode="numeric", length=n+w)
beta  <- vector(mode="numeric", length=n+w)
fit.par <- matrix(nrow=n+w, ncol=6) # 6 model parameters
mut  <- vector(mode="numeric", length=n+w)
sigt  <- vector(mode="numeric", length=n+w)


# ----------------------------------------------------
for(i in 1:(n+w))
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  nu=coef(fit)["shape"]
  
  qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu)
  
  # esmodel = NULL # expected shortfall for the assumed model/distribution
  # for(j in inu)
  #   esmodel = c(esmodel, integrate(function(x) x*dt(x, df=nu), qmodel[j], Inf)$value/(1-VaR.levels[j]))
  # esmodel=sqrt((nu-2)/nu)*esmodel # assuming nu>2
  
  esmodel = (dt(qt(nvec, df=nu), df=nu)) / (1-nvec) * (nu + (qt(nvec, df=nu))^2) / (nu-1)
  esmodel = sqrt((nu-2)/nu)*esmodel
  
  emodel = sqrt((nu-2)/nu)*et(asy=tvec, df=nu) # assuming nu>2
  
  foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
  
  VaR[i,]=foc$VaRmodel
  VaRfhs[i,] = foc$VaRfhs
  VaRevt[i,] = foc$VaRevt
  
  EXP[i,]=foc$EXPmodel
  EXPfhs[i,] = foc$EXPfhs
  EXPevt[i,] = foc$EXPevt
  
  ES[i,]=foc$ESmodel
  ESfhs[i,] = foc$ESfhs
  ESevt[i,] = foc$ESevt
  
  xi[i] = foc$evt.shape
  beta[i] = foc$evt.scale
  mut[i] = foc$mut
  sigt[i] = foc$sigt
  
  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)

save(out, file="Sim3std3.RDATA")


# =======================================================
# skewed Student t innovations
# The version of Fernandez & Steel (1998)
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")

VaR <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRfhs <- matrix(nrow=n+w,ncol=length(VaR.levels))
VaRevt <- matrix(nrow=n+w,ncol=length(VaR.levels))

ES <- matrix(nrow=n+w,ncol=length(nvec))
ESfhs <- matrix(nrow=n+w,ncol=length(nvec))
ESevt <- matrix(nrow=n+w,ncol=length(nvec))

EXP <- matrix(nrow=n+w,ncol=length(tvec))
EXPfhs <- matrix(nrow=n+w,ncol=length(tvec))
EXPevt <- matrix(nrow=n+w,ncol=length(tvec))

xi <- vector(mode="numeric", length=n+w)
beta  <- vector(mode="numeric", length=n+w)
fit.par <- matrix(nrow=n+w, ncol=7) # 7 model parameters
mut  <- vector(mode="numeric", length=n+w)
sigt  <- vector(mode="numeric", length=n+w)

for(i in 1:(n+w))
  
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  
  # mean and std dev of skewed t rv
  nu=coef(fit)["shape"]; ga=coef(fit)["skew"]
  
  # calculate VaR of normalized skewed-t
  qmodel = qsgt(VaR.levels,mu=0,sigma=1,lambda=(ga^2-1)/(ga^2+1),p=2,q=nu/2)
  
  esmodel = NULL # expected shortfall for the assumed model/distribution
  # for(j in inu)
  #   esmodel = c(esmodel, integrate(function(x) x*dskt(x, df=nu, gamma=ga), qmodel[j], Inf)$value/(1-VaR.levels[j]))
  # esmodel=(esmodel-m)/s
  
  # parameters of Hansen (1994)
  lam = -(ga^2-1)/(ga^2+1) # transfer from ga of Fernandez and Steel (1998) to lam of Hansen (1994)
  # minus sign is taken because we change from loss to return
  c = gamma((nu+1)/2) / (gamma(nu/2) * sqrt(pi * (nu-2)))
  a = 4*lam*c*(nu-2)/(nu-1)
  b = sqrt(1+3*lam^2-a^2)
  
  # Calculate explicit ES for skewed-t described in Patton et al. (2019)
  for(j in inu){
    # We make some transformation here because we are handling loss data
    if(qmodel[j] >= (a/b)){
      alpha_tilde1 = psgt(b/(1-lam)*(-qmodel[j]+a/b),mu=0,sigma=1,lambda=0,p=2,q=nu/2)
      es_t1 = sqrt((nu-2)/nu) * nu^(nu/2)/(2*(alpha_tilde1)*sqrt(pi))*gamma((nu-1)/2)/gamma(nu/2)*(qt(1-alpha_tilde1, df=nu)^2+nu)^((1-nu)/2) # ES for standardized t distribution
      esmodel = c(esmodel, -(alpha_tilde1)/(1-VaR.levels[j]) * (1-lam) * (-a/b - (1-lam)/b*es_t1))
    }else{
      lam2 = -lam
      ga2 = 1/ga
      nu2 = nu
      c2 = c
      a2 = 4*lam2*c2*(nu2-2)/(nu2-1)
      b2 = b
      alpha_tilde2 = psgt(b2/(1-lam2)*(qsgt(VaR.levels[j],mu=0,sigma=1,lambda=lam2,p=2,q=nu2/2)+a2/b2),mu=0,sigma=1,lambda=0,p=2,q=nu2/2)
      es_t2 = sqrt((nu2-2)/nu2) * nu2^(nu2/2)/(2*(alpha_tilde2)*sqrt(pi))*gamma((nu2-1)/2)/gamma(nu2/2)*(qt(1-alpha_tilde2, df=nu2)^2+nu2)^((1-nu2)/2)
      esmodel = c(esmodel, -(alpha_tilde2)/(1-VaR.levels[j]) * (1-lam2) * (-a2/b2 - (1-lam2)/b2*es_t2))
    }
  }
  
  m = mean.st(shape=nu,skew=ga)
  s = sqrt(var.st(shape=nu,skew=ga))
  
  emodel = (est(asy=tvec,shape=nu,skew=ga)-m)/s
  
  foc=RM.forecasts3(fit=fit,alpha=avec,tau=tvec,nu=nvec,qmodel=qmodel,emodel=emodel,esmodel=esmodel,seed=i)
  
  VaR[i,]=foc$VaRmodel
  VaRfhs[i,] = foc$VaRfhs
  VaRevt[i,] = foc$VaRevt
  
  EXP[i,]=foc$EXPmodel
  EXPfhs[i,] = foc$EXPfhs
  EXPevt[i,] = foc$EXPevt
  
  ES[i,]=foc$ESmodel
  ESfhs[i,] = foc$ESfhs
  ESevt[i,] = foc$ESevt
  
  xi[i] = foc$evt.shape
  beta[i] = foc$evt.scale
  mut[i] = foc$mut
  sigt[i] = foc$sigt
  
  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)

save(out, file="Sim3sstd3.RDATA")


