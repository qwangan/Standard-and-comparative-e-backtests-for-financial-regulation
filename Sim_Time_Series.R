library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library("plotrix")

source("Ecomp-Rfns.R")

n=5000 # out-of-sample size

w=500 # moving window size

# Levels for risk measures
avec=c(.90, .95, .99) # vector of alpha levels for VaR
tvec=c(.96561, .98761, .99855) # vector of tau levels for expectile
nvec=c(.754, .875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

n=n+w

# Parameters
nu=5; xi=1.5
mst = mean.st(shape=nu,skew=xi) # mean and sd for optimal forecast computations
sst = sqrt(var.st(shape=nu,skew=xi))

# recovering the mu[t] and sig[t] used in the data generation
# keeping values relevant for the optimal forecast computations
load("simdat.RDATA")
mut = tail((simdat$garch) - (simdat$sigma) * (simdat$eps),n)
sigt = tail(simdat$sigma,n)

simdat = simdat$garch
y=tail(simdat,n)


# --------------------------------------------------------
# Forecasts

#=====================================================
# VaR_alpha
#=====================================================

VaRout1 = getVaR(k=1, levels=avec,mut=mut,sigt=sigt)
VaRout2 = getVaR(k=2, levels=avec,mut=mut,sigt=sigt)
VaRout3 = getVaR(k=3, levels=avec,mut=mut,sigt=sigt)

#=====================================================
# tau-Expectile
#=====================================================

EXPout1 = getEXP(k=1, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
EXPout2 = getEXP(k=2, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
EXPout3 = getEXP(k=3, levels=tvec,mut=mut,sigt=sigt,mst=mst,sst=sst)

#=====================================================
# (VaR_nu, ES_nu)
#=====================================================

out=getVaRES(k=1, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout1b = out$VaR; ESout1 = out$ES

out=getVaRES(k=2, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout2b = out$VaR; ESout2 = out$ES

out=getVaRES(k=3, levels=nvec,mut=mut,sigt=sigt,mst=mst,sst=sst)
VaRout3b = out$VaR; ESout3 = out$ES

#=====================================================
# --------------------------------------------------------
# Volatility estimates

load("Sim3norm3.RDATA")
sigt.norm = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("Sim3std3.RDATA")
sigt.std = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("Sim3sstd3.RDATA")
sigt.sstd = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

sigt.mat = rbind(sigt.norm,sigt.std,sigt.sstd,sigt)


#---------------------------------------#
#
#      Comparative E-backtest
#
#---------------------------------------#
final_c = 0.5 # tuning parameter c for calculation of betting processes

method.name <- c("n-FP", "n-FHS", "n-EVT", "t-FP", "t-FHS", "t-EVT", "st-FP", "st-FHS", "st-EVT", "opt")
nm=dim(VaRout1)[1]; n=dim(VaRout1)[2]

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


##------------------------------------------------------------------------------------------------##

## VaR

## alpha = 0.9
for(i in 1:nm)
{
  smatVaR1[i,] <- sfVaR(r=VaRout1[i,],x=y,a=avec[1],h=1) # 1-homogeneous scoring functions
  smatVaR0[i,] <- sfVaR(r=VaRout1[i,],x=y,a=avec[1],h=0) # 0-homogeneous scoring functions
}

# Calculate the difference of VaR-score functions between methods

egrel.mat.VaR.max <- egrel.mat.VaR <- matrix(nrow = nm, ncol = nm) # matrix to store maximum and final e-values

options(warn = 0)

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
          gamma.temp <- (1-avec[1])*(r-r.star) - (r<=r.star)*(min(r,-M) - min(r.star,-M)) - (r>r.star)*(min(r,M) - min(r.star,M))
        }else{
          gamma.temp <- (1-avec[1])*(log(r)-log(r.star)) - (r>r.star)*(log((min(max(r.star,M),r))/(r.star)))
        }
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        
        
        ## GREL method
        
        sdiff.grel <- sfVaR(r=rep(r,w), x=y[(l-w):(l-1)], a=avec[1], h=sf.code) - sfVaR(r=rep(r.star,w), x=y[(l-w):(l-1)], a=avec[1], h=sf.code)
        
        
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

(egrel.mat.VaR1 <- t(egrel.mat.VaR.max))




##------------------------------------------------------------------------------------------------##

## VaR

## alpha = 0.99
for(i in 1:nm)
{
  smatVaR1[i,] <- sfVaR(r=VaRout3[i,],x=y,a=avec[3],h=1)
  smatVaR0[i,] <- sfVaR(r=VaRout3[i,],x=y,a=avec[3],h=0)
}

# Calculate the difference of VaR-score functions between methods

egrel.mat.VaR.max <- egrel.mat.VaR <- matrix(nrow = nm, ncol = nm)

options(warn = 0)

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
      # M <- max(abs(y),abs(VaRout3[c(i,j),]))
      
      # E value calculation
      
      for (l in (1+w):n) {
        
        r = VaRout3[i,l]
        r.star = VaRout3[j,l]
        
        ## gamma: upper bound of lambda
        if (sf.code == 1) {
          gamma.temp <- (1-avec[3])*(r-r.star) - (r<=r.star)*(min(r,-M) - min(r.star,-M)) - (r>r.star)*(min(r,M) - min(r.star,M))
        }else{
          gamma.temp <- (1-avec[3])*(log(r)-log(r.star)) - (r>r.star)*(log((min(max(r.star,M),r))/(r.star)))
        }
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        # GREL method
        
          sdiff.grel <- sfVaR(r=rep(r,w), x=y[(l-w):(l-1)], a=avec[3], h=sf.code) - sfVaR(r=rep(r.star,w), x=y[(l-w):(l-1)], a=avec[3], h=sf.code)
        
        
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

(egrel.mat.VaR3 <- t(egrel.mat.VaR.max))


##------------------------------------------------------------------------------------------------##
# Plotting the results
# postscript(file="plot_sim_tlm_VaR.eps", width=19, height= 4, horizontal=FALSE)
# par(mfrow=c(1,2))
# plotTLM.ecomp(egrel.mat.VaR1,rm = 1,lev = avec[1])
# plotTLM.ecomp(egrel.mat.VaR3,rm = 1,lev = avec[3])
# dev.off()

postscript(file="plot_sim_tlm_VaR_rej2.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaR1,rm = 1,lev = avec[1],thre=2)
plotTLM.ecomp.thre(egrel.mat.VaR3,rm = 1,lev = avec[3],thre=2)
dev.off()

postscript(file="plot_sim_tlm_VaR_rej5.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaR1,rm = 1,lev = avec[1],thre=5)
plotTLM.ecomp.thre(egrel.mat.VaR3,rm = 1,lev = avec[3],thre=5)
dev.off()

postscript(file="plot_sim_tlm_VaR_rej10.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaR1,rm = 1,lev = avec[1],thre=10)
plotTLM.ecomp.thre(egrel.mat.VaR3,rm = 1,lev = avec[3],thre=10)
dev.off()


##------------------------------------------------------------------------------------------------##

## Expectile

## tau = 0.96561

for(i in 1:nm)
{
  smatExp2[i,] <- sfEXP(r=EXPout1[i,],x=y,a=tvec[1],h=1) # 2-homogeneous scoring function
  smatExp0[i,] <- sfEXP(r=EXPout1[i,],x=y,a=tvec[1],h=0) # 0-homogeneous scoring function
}

egrel.mat.EXP.max <- egrel.mat.EXP <- matrix(nrow = nm, ncol = nm)

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
          gamma.temp <- (((r.star <= r) & (r <= M)) | ((r < r.star) & (r.star < -M)))*tvec[1]*(r^2 - r.star^2 - 2*M*abs(r-r.star)) +
            ((r >= r.star) & (r > M))*( (1-tvec[1])*(r^2 - r.star^2 - 2*M*(r-r.star)) + (1-2*tvec[1])*(max(M,r.star) - r.star)^2 ) +
            ((r < r.star) & (r.star >= -M))*( (1-tvec[1])*(r^2 - r.star^2 + 2*M*(r-r.star)) - (1-2*tvec[1])*(max(-M,r)-r)^2 )
          
        } else {
          gamma.temp <- (r <= r.star) * (1-tvec[1]) * (log(r/r.star)-M/r+M/r.star)
                      - ((r.star < M) & (r > M)) * (1 - 2*tvec[1]) * (log(M/r.star)+1-M/r.star)
                      + (r > r.star) * abs(1-tvec[1]-(r <= M)) * (log(r/r.star)+M/r-M/r.star)
        }
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        ## GREL method
        sdiff.grel <- sfEXP(r=rep(r,w), x=y[(l-w):(l-1)], a=tvec[1], h=sf.code) - sfEXP(r=rep(r.star,w), x=y[(l-w):(l-1)], a=tvec[1], h=sf.code)

        
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

(egrel.mat.EXP1 <- t(egrel.mat.EXP.max))



##------------------------------------------------------------------------------------------------##

## tau = 0.99855

for(i in 1:nm)
{
  smatExp2[i,] <- sfEXP(r=EXPout3[i,],x=y,a=tvec[3],h=1)
  smatExp0[i,] <- sfEXP(r=EXPout3[i,],x=y,a=tvec[3],h=0)
}

egrel.mat.EXP.max <- egrel.mat.EXP <- matrix(nrow = nm, ncol = nm)

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
      
      #M <- max(abs(y),abs(EXPout3[c(i,j),]))
      M <- max(abs(y))
      
      
      # E value calculation
      for (l in (1+w):n) {
        
        r = EXPout3[i,l]
        r.star = EXPout3[j,l]
        
        if (sf.code == 1) {
          gamma.temp <- (((r.star <= r) & (r <= M)) | ((r < r.star) & (r.star < -M)))*tvec[3]*(r^2 - r.star^2 - 2*M*abs(r-r.star)) +
            ((r >= r.star) & (r > M))*( (1-tvec[3])*(r^2 - r.star^2 - 2*M*(r-r.star)) + (1-2*tvec[3])*(max(M,r.star) - r.star)^2 ) +
            ((r < r.star) & (r.star >= -M))*( (1-tvec[3])*(r^2 - r.star^2 + 2*M*(r-r.star)) - (1-2*tvec[3])*(max(-M,r)-r)^2 )
        
          # optimization
          # gamma.fun <- function(x){sfEXP(r=r,x=x,a=tvec[3],h=1) - sfEXP(r=r.star,x=x,a=tvec[3],h=1)}
          # gamma.temp <- optimize(gamma.fun, c(-M,M), tol = 1e-6)$objective
        
        } else {
          gamma.temp <- (r <= r.star) * (1-tvec[3]) * (log(r/r.star)-M/r+M/r.star)
                      - ((r.star < M) & (r > M)) * (1 - 2*tvec[3]) * (log(M/r.star)+1-M/r.star)
                      + (r > r.star) * abs(1-tvec[3]-(r <= M)) * (log(r/r.star)+M/r-M/r.star)
        }
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        ## GREL method
        sdiff.grel <- sfEXP(r=rep(r,w), x=y[(l-w):(l-1)], a=tvec[3], h=sf.code) - sfEXP(r=rep(r.star,w), x=y[(l-w):(l-1)], a=tvec[3], h=sf.code)
      
        
        if(sum(sdiff.grel) == 0) {lambda.grel.temp[l-w] <- 0} else {
          lambda.grel.temp[l-w] <- (sum(sdiff.grel))/(sum(sdiff.grel^2))}
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

(egrel.mat.EXP3 <- t(egrel.mat.EXP.max))

##------------------------------------------------------------------------------------------------##
## Plotting the results
# postscript(file="plot_sim_tlm_EXP.eps", width=19, height= 4, horizontal=FALSE)
# par(mfrow=c(1,2))
# plotTLM.ecomp(egrel.mat.EXP1,rm = 2,lev = tvec[1])
# plotTLM.ecomp(egrel.mat.EXP3,rm = 2,lev = tvec[3])
# dev.off()

postscript(file="plot_sim_tlm_EXP_rej2.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.EXP1,rm = 2,lev = tvec[1],thre=2)
plotTLM.ecomp.thre(egrel.mat.EXP3,rm = 2,lev = tvec[3],thre=2)
dev.off()

postscript(file="plot_sim_tlm_EXP_rej5.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.EXP1,rm = 2,lev = tvec[1],thre=5)
plotTLM.ecomp.thre(egrel.mat.EXP3,rm = 2,lev = tvec[3],thre=5)
dev.off()

postscript(file="plot_sim_tlm_EXP_rej10.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.EXP1,rm = 2,lev = tvec[1],thre=10)
plotTLM.ecomp.thre(egrel.mat.EXP3,rm = 2,lev = tvec[3],thre=10)
dev.off()


##------------------------------------------------------------------------------------------------##



##------------------------------------------------------------------------------------------------##

## (VaR,ES)

## niu = 0.754
for(i in 1:nm)
{
  smatVaRES1[i,] <- sfVaRES(r1=VaRout1b[i,], r2=ESout1[i,], x=y, a=nvec[1], h=1) # 1/2-homogeneous scoring functions
  smatVaRES0[i,] <- sfVaRES(r1=VaRout1b[i,], r2=ESout1[i,], x=y, a=nvec[1], h=0) # 0-homogeneous scoring functions
}

egrel.mat.VaRES.max <- egrel.mat.VaRES <- matrix(nrow = nm, ncol = nm)
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
          else{gamma.temp <- (1-nvec[1])*((r+z)/(2*sqrt(r))-(r.star+z.star)/(2*sqrt(r.star)))+
            (r<=r.star)*((max(z,-M)-z)/(2*sqrt(r))-max((max(min(M,z),-M)-z.star)/(2*sqrt(r.star)),0))+
            (r>r.star)*min((M-min(z,M))/(2*sqrt(r))-(M-min(z.star,M))/(2*sqrt(r.star)),(max(z,-M)-z)/(2*sqrt(r)))}
        } else {
          gamma.temp <- (1-nvec[1]) * (z/r - z.star/r.star +log(r/r.star)) + (r<=r.star) * ((max(z,-M)-z)/r - max((max(min(M,z),-M)-z.star)/r.star,0))
                      + (r>r.star) * (min((M-min(z,M))/r-(M-min(z.star,M))/r.star, (max(z,-M)-z)/r))
        }
        
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        
        ## GREL method
        sdiff.grel <- sfVaRES(r1=rep(z,w), r2=rep(r,w), x=y[(l-w):(l-1)], a=nvec[1], h=sf.code) - sfVaRES(r1=rep(z.star,w), r2=rep(r.star,w), x=y[(l-w):(l-1)], a=nvec[1], h=sf.code)
        
        
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

(egrel.mat.VaRES1 <- t(egrel.mat.VaRES.max))

#-------------------------#
## niu = 0.975
for(i in 1:nm)
{
  smatVaRES1[i,] <- sfVaRES(r1=VaRout3b[i,], r2=ESout3[i,], x=y, a=nvec[3], h=1)
  smatVaRES0[i,] <- sfVaRES(r1=VaRout3b[i,], r2=ESout3[i,], x=y, a=nvec[3], h=0)
}
egrel.mat.VaRES.max <- egrel.mat.VaRES <- matrix(nrow = nm, ncol = nm)
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
      
      # M <- max(abs(y),abs(ESout3[c(i,j),]),abs(VaRout3b[c(i,j),]))
      M <- max(abs(y))
      
      # E value calculation
      
      for (l in (1+w):n) {
        r = max(ESout3[i,l],0)
        r.star = max(ESout3[j,l],0)
        z = VaRout3b[i,l]
        z.star = VaRout3b[j,l]
        
        # gamma: upper bound of lambda
        if (sf.code == 1) {
          if(r==0 | r.star==0){gamma.temp<-Inf}
          else{gamma.temp <- (1-nvec[3])*((r+z)/(2*sqrt(r))-(r.star+z.star)/(2*sqrt(r.star)))+
            (r<=r.star)*((max(z,-M)-z)/(2*sqrt(r))-max((max(min(M,z),-M)-z.star)/(2*sqrt(r.star)),0))+
            (r>r.star)*min((M-min(z,M))/(2*sqrt(r))-(M-min(z.star,M))/(2*sqrt(r.star)),(max(z,-M)-z)/(2*sqrt(r)))}
        } else {
          gamma.temp <- (1-nvec[3]) * (z/r - z.star/r.star +log(r/r.star)) + (r<=r.star) * ((max(z,-M)-z)/r - max((max(min(M,z),-M)-z.star)/r.star,0))
                      + (r>r.star) * (min((M-min(z,M))/r-(M-min(z.star,M))/r.star, (max(z,-M)-z)/r)) 
        }
        
        
        
        ## New method
        c <- final_c
        if (gamma.temp>=0) {gamma.vec[l-w] <- Inf} else {gamma.vec[l-w] <- -c/gamma.temp} 
        gamma <- gamma.vec[l-w]
        
        
        
        ## GREL method
        sdiff.grel <- sfVaRES(r1=rep(z,w), r2=rep(r,w), x=y[(l-w):(l-1)], a=nvec[3], h=sf.code) - sfVaRES(r1=rep(z.star,w), r2=rep(r.star,w), x=y[(l-w):(l-1)], a=nvec[3], h=sf.code)
        
        
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

(egrel.mat.VaRES3 <- t(egrel.mat.VaRES.max))

##------------------------------------------------------------------------------------------------##
## Plotting the results
# postscript(file="plot_sim_tlm_VaRES01.eps", width=19, height= 4, horizontal=FALSE)
# par(mfrow=c(1,2))
# plotTLM.ecomp(egrel.mat.VaRES1,rm = 3,lev = nvec[1])
# plotTLM.ecomp(egrel.mat.VaRES3,rm = 3,lev = nvec[3])
# dev.off()

postscript(file="plot_sim_tlm_VaRES01_rej2.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaRES1,rm = 3,lev = nvec[1],thre=2)
plotTLM.ecomp.thre(egrel.mat.VaRES3,rm = 3,lev = nvec[3],thre=2)
dev.off()

postscript(file="plot_sim_tlm_VaRES01_rej5.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaRES1,rm = 3,lev = nvec[1],thre=5)
plotTLM.ecomp.thre(egrel.mat.VaRES3,rm = 3,lev = nvec[3],thre=5)
dev.off()

postscript(file="plot_sim_tlm_VaRES01_rej10.eps", width=19, height= 4, horizontal=FALSE)
par(mfrow=c(1,2))
plotTLM.ecomp.thre(egrel.mat.VaRES1,rm = 3,lev = nvec[1],thre=10)
plotTLM.ecomp.thre(egrel.mat.VaRES3,rm = 3,lev = nvec[3],thre=10)
dev.off()


# Plotting e-processes

# VaR_0.99, n-FHS v.s. t-FHS
i=2; j=5

pdf(paste0(c("Sim_Eprocess/VaR/VaR99_ll_",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)

par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
par(mgp = c(1.5, 0.5, 0))
plotmat.tmp <- get(paste0("elist",".VaR"))
plotmat <- cbind(log(t(plotmat.tmp[[i]])[1:(n-w),j]), log(t(plotmat.tmp[[j]])[1:(n-w),i]), log(t(plotmat.tmp[[i]])[1:(n-w),j])-log(t(plotmat.tmp[[j]])[1:(n-w),i]))
matplot(1:(n-w), plotmat, type = "l", lty = 1,col = c("red","blue", "green"), ylim = c(-1,10),xlab = "number of data",ylab = "e-process (log scale)",main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
grid(nx = NA, ny = NULL)

abline(h=log(2),lty=2)
abline(h=log(5),lty=2)
abline(h=log(10),lty=2)
legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+")),expression(log("M"^"-"~"/"~"M"^"+"))),lty = 1,col =  c("red","blue","green"),bty = "n", cex = 0.85)

dev.off()

plot(1:(n-w), lambda.VaR.list[[i]][j,1:(n-w)],type = "l", lty = 1)
plot(1:(n-w), galist.VaR[[i]][j,1:(n-w)],type = "l", lty = 1)
plot(1:n, smatVaR1[i,] - smatVaR1[j,],type = "l", lty = 1, ylim = c(-0.01,.02))


# EXP_0.99855, t-EVT v.s. st-EVT
i=6; j=9

pdf(paste0(c("Sim_Eprocess/EXP/EXP998_lm_",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)

par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
par(mgp = c(1.5, 0.5, 0))
plotmat.tmp <- get(paste0("elist",".EXP"))
plotmat <- cbind(log(t(plotmat.tmp[[i]])[1:(n-w),j]), log(t(plotmat.tmp[[j]])[1:(n-w),i]), log(t(plotmat.tmp[[i]])[1:(n-w),j])-log(t(plotmat.tmp[[j]])[1:(n-w),i]))
matplot(1:(n-w), plotmat, type = "l", lty = 1,col = c("red","blue", "green"), ylim = c(-1,10),xlab = "number of data",ylab = "e-process (log scale)",main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
grid(nx = NA, ny = NULL)

abline(h=log(2),lty=2)
abline(h=log(5),lty=2)
abline(h=log(10),lty=2)
legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+")),expression(log("M"^"-"~"/"~"M"^"+"))),lty = 1,col =  c("red","blue","green"),bty = "n", cex = 0.85)

dev.off()


plot(1:(n-w), lambda.EXP.list[[i]][j,1:(n-w)],type = "l", lty = 1)
plot(1:(n-w), galist.EXP[[i]][j,1:(n-w)],type = "l", lty = 1)
plot(1:n, smatExp2[i,] - smatExp2[j,],type = "l", lty = 1, ylim = c(-0.01,.02))


# VaRES_0.975, n-FP v.s. n-EVT
i=1; j=3

pdf(paste0(c("Sim_Eprocess/VaRES/ES975_lr_",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)

par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
par(mgp = c(1.5, 0.5, 0))
plotmat.tmp <- get(paste0("elist",".VaRES"))
plotmat <- cbind(log(t(plotmat.tmp[[i]])[1:(n-w),j]), log(t(plotmat.tmp[[j]])[1:(n-w),i]), log(t(plotmat.tmp[[i]])[1:(n-w),j])-log(t(plotmat.tmp[[j]])[1:(n-w),i]))
matplot(1:(n-w), plotmat, type = "l", lty = 1,col = c("red","blue", "green"), ylim = c(-1,10),xlab = "number of data",ylab = "e-process (log scale)",main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
grid(nx = NA, ny = NULL)

abline(h=log(2),lty=2)
abline(h=log(5),lty=2)
abline(h=log(10),lty=2)
legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+")),expression(log("M"^"-"~"/"~"M"^"+"))),lty = 1,col =  c("red","blue","green"),bty = "n", cex = 0.85)

dev.off()
