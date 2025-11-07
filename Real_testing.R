# Illustration of backtesting and method comparison methodologies 
# using the NASDAQ Composite index


library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library("plotrix") 
library("lubridate")

source("Ecomp-Rfns.R")

# =======================================================
# COMPARATIVE E-BACKTESTING
# Real data analysis with NASDAQ index
# =======================================================

# Data

dat = rev(read.csv("nasdaq0325.csv")$Price) # NASDAQ Composite index (Jan 3, 2003- May 30, 2025)
dates = mdy(rev(read.csv("nasdaq0325.csv")$Date))[-1] # dates for returns
x = - log(dat[-1]/dat[-length(dat)])*100 # negated percentage log returns
N=length(x) # sample size
w0=500 # moving window size for forecasting
w=500 # moving window size for betting process calculation
(n=N-w) 

# risk measure levels
avec=.99 # vector of alpha levels for VaR
tvec=0.99855 # vector of tau levels for expectile
nvec=.975 # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own and in pair with ES
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR


load("NASDAQnorm0325.RDATA")
VaRout = rbind(out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(out$ES[,1],out$ESfhs[,1],out$ESevt[,1])
sigt.norm = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

#load("NASDAQstd.RDATA")
#VaRout = rbind(VaRout,out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
#EXPout = rbind(EXPout,out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
#VaRout2 = rbind(VaRout2,out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
#ESout = rbind(ESout,out$ES[,1],out$ESfhs[,1],out$ESevt[,1])
#sigt.sstd = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)

load("NASDAQsstd0325.RDATA")
VaRout = rbind(VaRout,out$VaR[,1],out$VaRfhs[,1],out$VaRevt[,1])
EXPout = rbind(EXPout,out$EXP[,1],out$EXPfhs[,1],out$EXPevt[,1])
VaRout2 = rbind(VaRout2,out$VaR[,2],out$VaRfhs[,2],out$VaRevt[,2])
ESout = rbind(ESout,out$ES[,1],out$ESfhs[,1],out$ESevt[,1])
sigt.sstd = matrix(rep(out$sigt,3), nrow=3,byrow=TRUE)


# Average forecasts and method ranks
nm=dim(VaRout)[1]; n=dim(VaRout)[2]
y = tail(x,n)
plot(y, type="l")
method.name <- c("n-FP", "n-FHS", "n-EVT", "st-FP", "st-FHS", "st-EVT")



#---------------------------------------#
#
#      Comparative E-backtest
#
#---------------------------------------#


# Matrices to store score values
smatVaRES1 <- smatVaRES0 <- smatExp2 <- smatExp0 <- smatVaR1 <-smatVaR0 <- matrix(nrow=nm,ncol=n)
#Matrices to store score values' difference
smatdiff <- matrix(nrow = nm*(nm-1), ncol = n)
# List to save gamma (i.e. upper bound)
galist.VaRES <- galist.EXP <- galist.VaR <- list()
for (i in 1:nm) {galist.VaRES[[i]] <- galist.EXP[[i]] <- galist.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save lambda
mono.test.list <- lambda.VaRES.list <-  lambda.VaR.list <- lambda.EXP.list <- list()
for (i in 1:nm) {mono.test.list[[i]] <- lambda.VaRES.list[[i]] <- lambda.VaR.list[[i]] <- lambda.EXP.list[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save e-process
elist.VaRES <- elist.EXP <- elist.VaR <- list()
for (i in 1:nm) {elist.VaRES[[i]] <- elist.EXP[[i]] <- elist.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)}
# List to save e-process with stopping rule
estop10.VaRES <- estop5.VaRES <- estop2.VaRES <- estop10.EXP <- estop5.EXP <- estop2.EXP <- estop10.VaR <- estop5.VaR <- estop2.VaR <- list()
for (i in 1:nm) 
  estop10.VaRES[[i]] <- estop5.VaRES[[i]] <- estop2.VaRES[[i]] <- estop10.EXP[[i]] <- estop5.EXP[[i]] <- estop2.EXP[[i]] <- estop10.VaR[[i]] <- estop5.VaR[[i]] <- estop2.VaR[[i]] <- matrix(0,nrow = nm, ncol = n-w)


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

stop.num <- 5

for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatVaR1[i,] - smatVaR1[j,]
    
    ##Adding the stopping rule in to the calculation
    #stop at 2, 5 or 10
    l = 1
    estop.temp <- vector()
    for (k in 1:(n-w)) {
      if (prod(1 + lambda.VaR.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)]) < stop.num & 
          prod(1 + lambda.VaR.list[[j]][i,l:k] * (-sdiff.temp[(l+w):(k+w)])) < stop.num) {
        estop.temp[k] <- prod(1 + lambda.VaR.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
      }else{
        estop.temp[k] <- prod(1 + lambda.VaR.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
        l = k + 1
      }}
    
    # estop2.VaR[[i]][j,1:(n-w)] <- estop.temp
    estop5.VaR[[i]][j,1:(n-w)] <- estop.temp
    # estop10.VaR[[i]][j,1:(n-w)] <- estop.temp
  }}

#plot the result
# method.name <- c("n-FP","n-FHS", "n-EVT", "st-FP","st-FHS", "st-EVT")
# for (i in 1:5) {
#   for (j in (i+1):nm) {
#     if (i==j) {
#       plot.new()
#     } else {
# 
#       pdf(paste0(c("NASDAQ 0325/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
#       par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2))
#       par(mgp = c(1.5, 0.5, 0))
#       #plotmat.tmp <- get(paste0("elist",".VaRES"))
#       #plotmat.tmp <- get(paste0("elist",".EXP"))
#       plotmat.tmp <- get(paste0("elist",".VaR"))
#       plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i],t(plotmat.tmp[[i]])[,j]-t(plotmat.tmp[[j]])[,i])
#       #plotmat <- cbind(t(elist.VaR[[i]])[,j], t(elist.VaR[[j]])[,i],t(elist.VaR[[i]])[,j]-t(elist.VaR[[j]])[,i])
#       matplot(tail(dates,n-w),plotmat, type = "l", lty = 1,col =  c("red","blue","green"),xaxt = "n",xlab = "year",ylab = "e-value",ylim = c(-1,4),main = paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""))
#       axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
#       grid(nx = NA, ny = NULL)
#       #abline(h=stop.num,lty=2)
#       abline(h=2,lty=2)
#       legend("topleft",legend = c(expression("M"^"-"),expression("M"^"+"),expression("M"^"-"~" - "~"M"^"+")),lty = 1,col =  c("red","blue","green"),bty = "n", cex = 0.85)
#       dev.off()
#     }}
# }



for (i in 1:5) {
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("Real_VaR stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",stop.num,".VaR"))
      plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i])
      matplot(tail(dates,n-w),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xaxt = "n",xlab = "year",ylab = "e-process (log scale)",
              main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
      grid(nx = NA, ny = NULL)
      
      for (k in 1:(nrow(plotmat) - 1)) {
        if (plotmat[k,1] <= stop.num & plotmat[k + 1, 1] <= stop.num & plotmat[k, 2] <= stop.num & plotmat[k + 1, 2] <= stop.num) {
          # Both points are <= stop.num: solid line
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),2]), col = "blue", lty = 1, lwd = 2)
        } else {
          # At least one point > stop.num: dashed line
          abline(v = tail(dates,n-w)[k], col = "grey", lty = 2, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,2]), col = "blue", lty = 1, lwd = 2)
        }
      }
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      # abline(h=log(10),lty=2)
      legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+"))),lty = 1,col = c("red", "blue"),bty = "n", cex = 1)
      
      dev.off()
    }}}
##------------------------------------------------------------------------------------------------##



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

stop.num <- 5

for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatExp2[i,] - smatExp2[j,]
    
    # Adding the stopping rule in to the calculation
    # stop at 2, 5 or 10
    l = 1
    estop.temp <- vector()
    for (k in 1:(n-w)) {
      if (prod(1 + lambda.EXP.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)]) < stop.num & 
          prod(1 + lambda.EXP.list[[j]][i,l:k] * (-sdiff.temp[(l+w):(k+w)])) < stop.num) {
        estop.temp[k] <- prod(1 + lambda.EXP.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
      }else{
        estop.temp[k] <- prod(1 + lambda.EXP.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
        l = k + 1
      }}
    # estop2.EXP[[i]][j,1:(n-w)] <- estop.temp
    estop5.EXP[[i]][j,1:(n-w)] <- estop.temp
    # estop10.EXP[[i]][j,1:(n-w)] <- estop.temp
  }}


# plot the result

# for (i in 1:5) {
#   for (j in (i+1):nm) {
#     if (i==j) {
#       plot.new()
#     } else {
#       
#       #pdf(paste0(c("NASDAQ 0323/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
#       
#       par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
#       par(mgp = c(1.5, 0.5, 0))
#       #plotmat.tmp <- get(paste0("elist",".VaRES"))
#       plotmat.tmp <- get(paste0("elist",".EXP"))
#       #plotmat.tmp <- get(paste0("elist",".VaR"))
#       plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i],t(plotmat.tmp[[i]])[,j]-t(plotmat.tmp[[j]])[,i])
#       #plotmat <- cbind(t(elist.VaR[[i]])[,j], t(elist.VaR[[j]])[,i],t(elist.VaR[[i]])[,j]-t(elist.VaR[[j]])[,i])
#       matplot(tail(dates,n-w),plotmat, type = "l", lty = 1,col =  c("red","blue","green"),xaxt = "n",xlab = "year",ylab = "e-value",ylim = c(-1.5,3),main = paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""))
#       axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
#       grid(nx = NA, ny = NULL)
#       #abline(h=stop.num,lty=2)
#       abline(h=2,lty=2)
#       legend("topleft",legend = c(expression("M"^"-"),expression("M"^"+"),expression("M"^"-"~" - "~"M"^"+")),lty = 1,col =  c("red","blue","green"),bty = "n", cex = 0.85)
#       
#       
#       #dev.off()
#     }}
# }

for (i in 1:5) {
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("Real_EXP stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",stop.num,".EXP"))
      plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i])
      matplot(tail(dates,n-w),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xaxt = "n",xlab = "year",ylab = "e-process (log scale)",main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
      grid(nx = NA, ny = NULL)
      
      for (k in 1:(nrow(plotmat) - 1)) {
        if (plotmat[k,1] <= stop.num & plotmat[k + 1, 1] <= stop.num & plotmat[k, 2] <= stop.num & plotmat[k + 1, 2] <= stop.num) {
          # Both points are <= stop.num: solid line
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),2]), col = "blue", lty = 1, lwd = 2)
        } else {
          # At least one point > stop.num: dashed line
          abline(v = tail(dates,n-w)[k], col = "grey", lty = 2, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,2]), col = "blue", lty = 1, lwd = 2)
        }
      }
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      # abline(h=log(10),lty=2)
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


stop.num <- 5

for (i in 1:nm) {
  for (j in (1:nm)[-i]) {
    
    sdiff.temp <- smatVaRES1[i,] - smatVaRES1[j,]
    
    ##Adding the stopping rule in to the calculation
    #stop at 2
    l = 1
    estop.temp <- vector()
    for (k in 1:(n-w)) {
      if (prod(1 + lambda.VaRES.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)]) < stop.num & 
          prod(1 + lambda.VaRES.list[[j]][i,l:k] * (-sdiff.temp[(l+w):(k+w)])) < stop.num) {
        estop.temp[k] <- prod(1 + lambda.VaRES.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
      }else{
        estop.temp[k] <- prod(1 + lambda.VaRES.list[[i]][j,l:k] * sdiff.temp[(l+w):(k+w)])
        l = k + 1
      }}
    # estop2.VaRES[[i]][j,1:(n-w)] <- estop.temp
    estop5.VaRES[[i]][j,1:(n-w)] <- estop.temp
    # estop10.VaRES[[i]][j,1:(n-w)] <- estop.temp
  }}


# plot the result

# for (i in 1:nm) {
#   stop.num = 2
#   #pdf(paste0(c("NASDAQ 0623/",method.name[i],"_VaRES_0623_s",stop.num,".pdf"),collapse = ""), width = 16.56/2, height = 23.2*0.33)
#   plotmat.tmp <- get(paste0("estop",stop.num,".VaRES"))
#   plotmat <- t(plotmat.tmp[[i]])[,-i]
#   
#   if (max(plotmat) <= 1.2*stop.num) {
#     matplot(tail(dates,n-w),plotmat,type = "l", ylim = c(0,1.2*stop.num),lty = 1,col = c(1:nm)[-i],xaxt = "n",xlab = "year",ylab = "e-value",main = method.name[i])
#   }else{
#     matplot(tail(dates,n-w),plotmat,type = "l",lty = 1,col = c(1:nm)[-i],xaxt = "n",xlab = "year",ylab = "e-value",main = method.name[i])
#   }
#   
#   axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
#   grid(nx = NA, ny = NULL)
#   abline(h=stop.num,lty=2)
#   legend("topleft",legend = method.name[-i],lty = 1,col = c(1:nm)[-i],bty = "n",cex = 1)
#   #dev.off()
# }

for (i in 1:5) {
  for (j in (i+1):nm) {
    if (i==j) {
      plot.new()
    } else {
      
      pdf(paste0(c("Real_VaRES stop/",paste0(c(method.name[i], " vs ",method.name[j]), collapse = ""),".pdf"),collapse = ""), width = 16.56/3, height = 23.2*0.18)
      
      par(mar = c(2.5, 2.5, 1.5, 0), oma = c(2, 2, 1, 2) ) 
      par(mgp = c(1.5, 0.5, 0))
      plotmat.tmp <- get(paste0("estop",stop.num,".VaRES"))
      plotmat <- cbind(t(plotmat.tmp[[i]])[,j], t(plotmat.tmp[[j]])[,i])
      matplot(tail(dates,n-w),log(plotmat), type = "l", lty = 1,col = c("white","white"), ylim = c(-1,2), xaxt = "n",xlab = "year",ylab = "e-process (log scale)",main = paste0(c(method.name[i], " v.s. ",method.name[j]), collapse = ""))
      axis.Date(1, at = seq(min(tail(dates,n-w)), max(tail(dates,n-w)), by = "2 years"), format = "%Y")
      grid(nx = NA, ny = NULL)
      
      for (k in 1:(nrow(plotmat) - 1)) {
        if (plotmat[k,1] <= stop.num & plotmat[k + 1, 1] <= stop.num & plotmat[k, 2] <= stop.num & plotmat[k + 1, 2] <= stop.num) {
          # Both points are <= stop.num: solid line
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[k:(k + 1)], log(plotmat[k:(k + 1),2]), col = "blue", lty = 1, lwd = 2)
        } else {
          # At least one point > stop.num: dashed line
          abline(v = tail(dates,n-w)[k], col = "grey", lty = 2, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,1]), col = "red", lty = 1, lwd = 2)
          lines(tail(dates,n-w)[(k-1):k], log(plotmat[(k-1):k,2]), col = "blue", lty = 1, lwd = 2)
        }
      }
      
      abline(h=log(2),lty=2)
      abline(h=log(5),lty=2)
      # abline(h=log(10),lty=2)
      legend("topleft",legend = c(expression(log("M"^"-")),expression(log("M"^"+"))),lty = 1,col = c("red", "blue"),bty = "n", cex = 1)
      
      dev.off()
    }}}





