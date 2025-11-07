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
source("Ecomp-Rfns.R")


# =======================================================
# structural change point b*=2000
# The simulated series is saved as "structural_change.RDATA"

n=4000
w1=500 # time window for forecast
w=500 # time window for grel


# =======================================================
# SIMULATION SET-UP 1
# Synthetic dataset generated from an AR(1)-GARCH(1,1) process with normal innovations
# AR-GARCH filter parameters:
# mu= -0.05; ar1 = 0.1 # AR(1) part
# omega=.3; al=.01, be=.1 + .7 * (t > b*) # GARCH(1,1) parameters
# Burn-in period of 1000 points was used


# set.seed(233)
# spec = garchSpec(model = list(mu = -0.05, ar = 0.1, omega = .3, alpha = .01, beta = 0.1),
#                  cond.dist = "norm")
# simdat1 = garchSim(spec, n = n/2+(w+w1), n.start = 1000, extended = TRUE)
# simdat1 = simdat1$garch

# spec = garchSpec(model = list(mu = -0.05, ar = 0.1, omega = .3, alpha = .01, beta = 0.9),
#                  cond.dist = "norm")
# simdat2 = garchSim(spec, n = n/2, n.start = 1000, extended = TRUE)
# simdat2 = simdat2$garch


# SIMULATION SET-UP 2
# Synthetic dataset generated from an AR(1)-GARCH(1,1) process with skewed-t innovations
# AR-GARCH filter parameters:
# mu= -0.05; ar1 = 0.1 # AR(1) part
# omega=.3; al=.1, be=.5 # GARCH(1,1) parameters
# Innovation distribution parameters:
# nu=6 - 3 * (t > b*) # shape parameter
# ga=1 # skewness parameter
# Burn-in period of 1000 points was used

set.seed(233)
spec = garchSpec(model = list(mu = -0.05, ar = 0.1, omega = .3, alpha = 0.1, beta = 0.5, skew = 1, shape = 6),
                 cond.dist = "sstd")
simdat1 = garchSim(spec, n = n/2+(w+w1), n.start = 1000, extended = TRUE)
simdat1 = simdat1$garch

spec = garchSpec(model = list(mu = -0.05, ar = 0.1, omega = .3, alpha = .1, beta = 0.5, skew = 1, shape = 3),
                 cond.dist = "sstd")
simdat2 = garchSim(spec, n = n/2, n.start = 1000, extended = TRUE)
simdat2 = simdat2$garch


simdat = c(simdat1,simdat2)
plot(1:(n),tail(simdat, n),type="l")
save(simdat, file = "structural_change.RDATA")

# =======================================================

load("structural_change.RDATA")

avec=.9 # vector of alpha levels for VaR
tvec=.99855 # vector of tau levels for expectile
nvec=.975 # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

x=tail(simdat,n+w1+w) # time series to be used for fitting and forecasting

# =======================================================
# Normal innovations
# =======================================================
spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1) # VaR values

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
  fit = ugarchfit(spec, x[i:(i+w1-1)], solver="hybrid")
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
  
  if(i %% 5 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)
save(out, file="SC_norm.RDATA")

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
  fit = ugarchfit(spec, x[i:(i+w1-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  
  # mean and std dev of skewed t rv
  nu=coef(fit)["shape"]; ga=coef(fit)["skew"]
  
  # calculate VaR of normalized skewed-t
  qmodel = qsgt(VaR.levels,mu=0,sigma=1,lambda=(ga^2-1)/(ga^2+1),p=2,q=nu/2)
  
  esmodel = NULL # expected shortfall for the assumed model/distribution
  
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
  
  if(i %% 5 == 0) print(i)
}

out=list(VaR=VaR,VaRfhs=VaRfhs,VaRevt=VaRevt,EXP=EXP,EXPfhs=EXPfhs,EXPevt=EXPevt,ES=ES,ESfhs=ESfhs,ESevt=ESevt,xi=xi, beta=beta,par=fit.par,mut=mut,sigt=sigt)

save(out, file="SC_sstd.RDATA")


load("structural_change.RDATA")


# --------------------------------------------------------
# Forecasts
VaRout<-NULL
EXPout<-NULL
VaRout2<-NULL
ESout<-NULL

load("SC_norm.RDATA")
VaRout = rbind(out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(out$ES[,1],out$ESfhs[,1],out$ESevt[,1])

load("SC_sstd.RDATA")
VaRout = rbind(out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(out$ES[,1],out$ESfhs[,1],out$ESevt[,1])

# Average forecasts and method ranks
nm=dim(VaRout)[1]; n=dim(VaRout)[2]
y = tail(simdat,n)
plot(y, type="l")
method.name <- c("n-FP", "n-FHS", "n-EVT", "st-FP", "st-FHS", "st-EVT")


#---------------------------------------#
#
#      Comparative E-backtest
#
#---------------------------------------#

final_c = 0.1 # tuning parameter c for calculation of betting processes

# Matrices to store score values
smatVaRES1 <- smatVaRES0 <- smatExp2 <- smatExp0 <- smatVaR1 <-smatVaR0 <- matrix(nrow=nm,ncol=n)
#Matrices to store score values' difference
smatdiff <- matrix(nrow = nm*(nm-1), ncol = n)
# List to save gamma (i.e. upper bound)
galist.VaRES <- galist.EXP <- galist.VaR <- list()
for (i in 1:nm) {galist.VaRES[[i]] <- galist.EXP[[i]] <- galist.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save lambda
lambda.VaRES.list <-  lambda.VaR.list <- lambda.EXP.list <- list()
for (i in 1:nm) {lambda.VaRES.list[[i]] <- lambda.VaR.list[[i]] <- lambda.EXP.list[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save e-process
elist.VaRES <- elist.EXP <- elist.VaR <- list()
for (i in 1:nm) {elist.VaRES[[i]] <- elist.EXP[[i]] <- elist.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save e-process with stopping rule
estop.VaRES <- estop.EXP <- estop.VaR <- list()
for (i in 1:nm) 
  estop.VaRES[[i]] <- estop.EXP[[i]] <- estop.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)


# Matrices to store e-values of two-sided calibration tests
for(i in 1:nm)
{
  smatVaR1[i,] <- sfVaR(r=VaRout[i,],x=y,a=avec,h=1)
  smatVaR0[i,] <- sfVaR(r=VaRout[i,],x=y,a=avec,h=0)
  smatExp2[i,] <- sfEXP(r=EXPout[i,],x=y,a=tvec,h=1)
  smatExp0[i,] <- sfEXP(r=EXPout[i,],x=y,a=tvec,h=0)
  smatVaRES1[i,] <- sfVaRES(r1=VaRout2[i,], r2=ESout[i,], x=y, a=nvec, h=1)
  smatVaRES0[i,] <- sfVaRES(r1=VaRout2[i,], r2=ESout[i,], x=y, a=nvec, h=0)
}

# Calculate the difference of VaR-score functions between methods

egrel.mat.VaR.max <- egrel.mat.VaR <- matrix(nrow = nm, ncol = nm) # matrix to store maximum and final e-values

options(warn = 0)

# run to the end

system.time(
  for (i in 1:nm) {
    for (j in (1:nm)[-i]) {
      
      ##'@Notice: please adjust all the score functions including "sdiff.temp" & "sdiff.grel"
      
      sdiff.temp <- smatVaR1[i,] - smatVaR1[j,]
      #sdiff.temp <- smatVaR0[i,] - smatVaR0[j,]
      
      sf.code = 1 #score function = 1 or 0
      
      
      lambda.grel.temp <- vector()
      eprocess.grel.temp <- vector()
      
      e.temp <- 1
      gamma.vec <- vector()
      
      M <- max(abs(y))
      #M <- max(abs(y),abs(VaRout1[c(i,j),]))
      
      # E value calculation
      
      for (l in (1+w):n) {
        
        r = VaRout1[i,l]
        r.star = VaRout1[j,l]
        
        # gamma: upper bound of lambda
        if (sf.code == 1) {
          gamma.temp <- (1-avec)*(r-r.star) - (r<=r.star)*(min(r,-M) - min(r.star,-M)) - (r>r.star)*(min(r,M) - min(r.star,M))
        }else{
          gamma.temp <- (1-avec)*(log(r)-log(r.star)) - (r>r.star)*(log((min(max(r.star,M),r))/(r.star)))
        }
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        
        
        ## GREL method
        
        sdiff.grel <- sfVaR(r=rep(r,w), x=y[(l-w):(l-1)], a=avec, h=sf.code) - sfVaR(r=rep(r.star,w), x=y[(l-w):(l-1)], a=avec, h=sf.code)
        
        
        ## New method
        if(sum(sdiff.grel) == 0) {lambda.grel.temp[l-w] <- 0} else {
          lambda.grel.temp[l-w] <- (sum(sdiff.grel))/(sum(sdiff.grel^2))}
        lambda.grel.temp[l-w] <- min(max(lambda.grel.temp[l-w],0),gamma)
        
        
        e.temp <- e.temp*(1 + lambda.grel.temp[l-w]*sdiff.temp[l])
        eprocess.grel.temp[l-w] <- e.temp
      }
      
      
      lambda.VaR.list[[i]][j,1:(n-w)] <- lambda.grel.temp
      elist.VaR[[i]][j,1:(n-w)] <- eprocess.grel.temp
      galist.VaR[[i]][j,1:(n-w)] <- gamma.vec
      egrel.mat.VaR[i,j] = eprocess.grel.temp[n-w]
      egrel.mat.VaR.max[i,j] = max(eprocess.grel.temp)
    }}
)

# run again adding stopping rule at the structural change point

for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatVaR1[i,] - smatVaR1[j,]
    
    ##Adding the stopping rule in to the calculation
    estop.temp <- vector()
    
    ## restart at middle
    estop.temp[1:((n-w)/2)] <- elist.VaR[[i]][j, 1:((n-w)/2)]
    estop.temp[(n-w)/2+1] <- 1
    for (k in ((n-w)/2+2):(n-w)) {estop.temp[k] <- estop.temp[k-1]*(1 + lambda.VaR.list[[i]][j,k] * sdiff.temp[k+w])}
    
    
    estop.VaR[[i]][j,1:(n-w)] <- estop.temp
  }}


# Plot the result


for (i in 1:(nm-1)) {
  stop.num = 2
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("SC_VaR stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",".VaR"))
      plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i])
      matplot(c(1:(n-w)),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xlab = "number of data", ylab = "e-process (log scale)",
              main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),1]),type = "l",col = "red")
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),2]),type = "l",col = "blue")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),1]),type = "l",col = "red")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),2]),type = "l",col = "blue")
      grid(nx = NA, ny = NULL)
      
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      abline(v=(n-w)/2+1,lty=2)
      legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+"))),lty = 1,col = c("red", "blue"),bty = "n", cex = 1)
      dev.off()
    }}}



##------------------------------------------------------------------------------------------------##

## Expectile

egrel.mat.EXP.max <- egrel.mat.EXP <- matrix(nrow = nm, ncol = nm)

# run to the end

#options(warn = 2)
system.time(
  for (i in 1:nm) {
    for (j in (1:nm)[-i]) {
      
      sdiff.temp <- smatExp2[i,] - smatExp2[j,]
      #sdiff.temp <- smatExp0[i,] - smatExp0[j,]
      
      sf.code = 1 #score function = 1 or 0
      
      
      lambda.grel.temp <- vector()
      eprocess.grel.temp <- vector()
      
      e.temp <- 1
      
      gamma.vec <- vector()
      
      # M <- max(abs(y),abs(EXPout1[c(i,j),]))
      M <- max(abs(y))
      
      # E value calculation
      for (l in (1+w):n) {
        
        r = EXPout1[i,l]
        r.star = EXPout1[j,l]
        
        if (sf.code == 1) {
          # gamma: upper bound of lambda
          gamma.temp <- (((r.star <= r) & (r <= M)) | ((r < r.star) & (r.star < -M)))*tvec*(r^2 - r.star^2 - 2*M*abs(r-r.star)) +
            ((r >= r.star) & (r > M))*( (1-tvec)*(r^2 - r.star^2 - 2*M*(r-r.star)) + (1-2*tvec)*(max(M,r.star) - r.star)^2 ) +
            ((r < r.star) & (r.star >= -M))*( (1-tvec)*(r^2 - r.star^2 + 2*M*(r-r.star)) - (1-2*tvec)*(max(-M,r)-r)^2 )
          
        } else {
          gamma.temp <- (r <= r.star) * (1-tvec) * (log(r/r.star)-M/r+M/r.star)
          - ((r.star < M) & (r > M)) * (1 - 2*tvec) * (log(M/r.star)+1-M/r.star)
          + (r > r.star) * abs(1-tvec-(r <= M)) * (log(r/r.star)+M/r-M/r.star)
        }
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        ## GREL method
        sdiff.grel <- sfEXP(r=rep(r,w), x=y[(l-w):(l-1)], a=tvec, h=sf.code) - sfEXP(r=rep(r.star,w), x=y[(l-w):(l-1)], a=tvec, h=sf.code)
        
        
        lambda.grel.temp[l-w] <- (sum(sdiff.grel))/(sum(sdiff.grel^2))
        lambda.grel.temp[l-w] <- min(max(lambda.grel.temp[l-w],0),gamma)
        
        e.temp <- e.temp*(1 + lambda.grel.temp[l-w]*sdiff.temp[l])
        eprocess.grel.temp[l-w] <- e.temp
        
      }
      
      
      lambda.EXP.list[[i]][j,1:(n-w)]  <- lambda.grel.temp
      galist.EXP[[i]][j,1:(n-w)] <- gamma.vec
      elist.EXP[[i]][j,1:(n-w)] <- eprocess.grel.temp
      egrel.mat.EXP[i,j] = eprocess.grel.temp[n-w]
      egrel.mat.EXP.max[i,j] = max(eprocess.grel.temp)
    }}
)


# run again adding stopping rule at the structural change point


for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatExp2[i,] - smatExp2[j,]
    estop.temp <- vector()
    
    # restart at middle
    estop.temp[1:(n/2)] <- elist.EXP[[i]][j, 1:(n/2)]
    estop.temp[n/2+1] <- 1
    for (k in (n/2+2):n) {estop.temp[k] <- estop.temp[k-1]*(1 + lambda.EXP.list[[i]][j,k] * sdiff.temp[k])}
    
    
    estop.EXP[[i]][j,1:n] <- estop.temp
    
    
  }}


# plot the result

for (i in 1:5) {
  stop.num = 2
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("SC_EXP stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",".EXP"))
      plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i])
      matplot(c(1:(n-w)),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xlab = "number of data", ylab = "e-process (log scale)",
              main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),1]),type = "l",col = "red")
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),2]),type = "l",col = "blue")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),1]),type = "l",col = "red")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),2]),type = "l",col = "blue")
      grid(nx = NA, ny = NULL)
      
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      abline(v=(n-w)/2+1,lty=2)
      legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+"))),lty = 1,col = c("red", "blue"),bty = "n", cex = 1)
      dev.off()
    }}}


##------------------------------------------------------------------------------------------------##



##------------------------------------------------------------------------------------------------##

## (VaR,ES)

egrel.mat.VaRES.max <- egrel.mat.VaRES <- matrix(nrow = nm, ncol = nm)


#run to the end

options(warn = 0)
system.time(
  for (i in 1:nm) {
    for (j in (1:nm)[-i]) {
      
      sdiff.temp <- smatVaRES1[i,] - smatVaRES1[j,]
      #sdiff.temp <- smatVaRES0[i,] - smatVaRES0[j,]
      
      sf.code = 1 #score function = 1 or 0
      
      
      lambda.grel.temp <- vector()
      eprocess.grel.temp <- vector()
      
      e.temp <- 1
      gamma.vec <- vector()
      
      M <- max(abs(y))
      #M <- max(abs(y),abs(ESout1[c(i,j),]),abs(VaRout1b[c(i,j),]))
      
      # E value calculation
      
      for (l in (1+w):n) {
        r = max(ESout1[i,l],0)
        r.star = max(ESout1[j,l],0)
        z = VaRout1b[i,l]
        z.star = VaRout1b[j,l]
        
        ## gamma: upper bound of lambda
        if (sf.code == 1) {
          if(r==0 | r.star==0){gamma.temp<-Inf}
          else{gamma.temp <- (1-nvec)*((r+z)/(2*sqrt(r))-(r.star+z.star)/(2*sqrt(r.star)))+
            (r<=r.star)*((max(z,-M)-z)/(2*sqrt(r))-max((max(min(M,z),-M)-z.star)/(2*sqrt(r.star)),0))+
            (r>r.star)*min((M-min(z,M))/(2*sqrt(r))-(M-min(z.star,M))/(2*sqrt(r.star)),(max(z,-M)-z)/(2*sqrt(r)))}
        } else {
          gamma.temp <- (1-nvec) * (z/r - z.star/r.star +log(r/r.star)) + (r<=r.star) * ((max(z,-M)-z)/r - max((max(min(M,z),-M)-z.star)/r.star,0))
          + (r>r.star) * (min((M-min(z,M))/r-(M-min(z.star,M))/r.star, (max(z,-M)-z)/r))
        }
        
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        
        ## GREL method
        sdiff.grel <- sfVaRES(r1=rep(z,w), r2=rep(r,w), x=y[(l-w):(l-1)], a=nvec, h=sf.code) - sfVaRES(r1=rep(z.star,w), r2=rep(r.star,w), x=y[(l-w):(l-1)], a=nvec, h=sf.code)
        
        
        if(sum(sdiff.grel) == 0) {lambda.grel.temp[l-w] <- 0} else {
          lambda.grel.temp[l-w] <- (sum(sdiff.grel))/(sum(sdiff.grel^2))}
        lambda.grel.temp[l-w] <- min(max(lambda.grel.temp[l-w],0),gamma)
        
        e.temp <- e.temp*(1 + lambda.grel.temp[l-w]*sdiff.temp[l])
        eprocess.grel.temp[l-w] <- e.temp
        
      }
      
      
      lambda.VaRES.list[[i]][j,1:(n-w)] <- lambda.grel.temp
      elist.VaRES[[i]][j,1:(n-w)] <- eprocess.grel.temp
      galist.VaRES[[i]][j,1:(n-w)] <- gamma.vec
      egrel.mat.VaRES[i,j] = eprocess.grel.temp[n-w]
      egrel.mat.VaRES.max[i,j] = max(eprocess.grel.temp)
      
      
    }}
)


# run again adding stopping rule at the structural change point

for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatVaRES1[i,] - smatVaRES1[j,]
    estop.temp <- vector()
    
    
    # restart at middle
    estop.temp[1:(n/2)] <- elist.VaRES[[i]][j, 1:(n/2)]
    estop.temp[n/2+1] <- 1
    for (k in (n/2+2):n) {estop.temp[k] <- estop.temp[k-1]*(1 + lambda.VaRES.list[[i]][j,k] * sdiff.temp[k])}
    
    
    estop.VaRES[[i]][j,1:n] <- estop.temp
    
    
  }}


for (i in 1:5) {
  stop.num = 2
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("SC_VaRES stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",".VaRES"))
      matplot(c(1:(n-w)),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xlab = "number of data", ylab = "e-process (log scale)",
              main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),1]),type = "l",col = "red")
      lines(1:((n-w)/2), log(plotmat[1:((n-w)/2),2]),type = "l",col = "blue")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),1]),type = "l",col = "red")
      lines(((n-w)/2+1):(n-w), log(plotmat[((n-w)/2+1):(n-w),2]),type = "l",col = "blue")
      grid(nx = NA, ny = NULL)
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      abline(v=(n-w)/2+1,lty=2)
      legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+"))),lty = 1,col = c("red", "blue"),bty = "n", cex = 1)
      
      dev.off()
    }}}



