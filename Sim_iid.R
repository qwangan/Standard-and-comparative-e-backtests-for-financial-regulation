set.seed(20251007)

##  parameter
l <- 10 # training set
n <- 1000 # data set
Tlen <- l + n
p <- 0.95 # confidence levek
gamma <- 0.5
B <- 1000      # repetition
thresh <- 5    # threshold

## e value for ES 
e_ES <- function(x, r, z, p = 0.95) {
  num <- pmax(x - z, 0)
  den <- (1 - p) * (r - z)
  if (length(den) == 1L) den <- rep(den, length(num))
  ifelse(den > 0, num / den, 0)
}

## combination of underestimate VaR or/and ES
param_sets <- list(
  c(k = 0.95, l_mult = 1.00),
  c(k = 0.9, l_mult = 1.00),  # underestimate VaR
  c(k = 1.00, l_mult = 0.95),
  c(k = 1.00, l_mult = 0.90),  # underestimate ES
  c(k = 0.95, l_mult = 0.95), # underestimate both
  c(k = 1.00, l_mult = 1.00)  # baseline
)

## save e-process and rejection
avg_list <- vector("list", length(param_sets))
rej_list <- numeric(length(param_sets))
labels   <- character(length(param_sets))


for (idx in seq_along(param_sets)) {
  k <- param_sets[[idx]]["k"]
  l_mult <- param_sets[[idx]]["l_mult"]
  
  # save M_t
  M_mat <- matrix(NA_real_, nrow = Tlen, ncol = B)
  rejects <- numeric(B)
  
  for (b in 1:B) {
    # date generated process
    vals <- c(0, 0.1 * 1:5, -0.1 * 1:5)
    eps  <- sample(vals, size = Tlen, replace = TRUE) # uniform
    
    z <- 1.64 * k + eps
    r <- 2.06 * l_mult + eps
    L <- rnorm(Tlen, mean = 0, sd = 1)
    
    M <- 1.0 # initial e-process
    M_path <- rep(1.0, Tlen)
    
    for (t in (l + 1):Tlen) {
      # grel method
      e_vec <- e_ES(L[1:(t - 1)], r = r[t], z = z[t], p = p)
      num <- sum(e_vec) - (t - 1)
      den <- sum((e_vec - 1)^2)
      ratio <- if (den > 0) num / den else 0
      lambda_t <- min(gamma, max(0, ratio))
      
      X_t <- e_ES(L[t], r = r[t], z = z[t], p = p)
      
      # e-process
      M <- M * (1 - lambda_t + lambda_t * X_t)
      M_path[t] <- M
    }
    
    M_mat[, b] <- M_path
    rejects[b] <- as.integer(M_path[Tlen] > thresh)
  }
  
  # average e-process and rejection
  M_avg <- rowMeans(log(M_mat)) # log scale e-process
  rej_rate <- mean(rejects)
  
  avg_list[[idx]] <- M_avg
  rej_list[idx] <- rej_rate
  labels[idx] <- sprintf("k=%.2f, l=%.2f (rej=%.3f)", k, l_mult, rej_rate)
}




## ---------- plot ----------
cols <- c("black","#d62728","#2ca02c","#9467bd","#ff7f0e","#17becf")
plot(1:Tlen, avg_list[[6]], type = "l", lwd = 1.6, col = cols[1],
     xlab = "number of data (t)",
     ylab = "e-process (log-scale)"
     #main = sprintf("Average log e-process over %d simulations", B)
     )

abline(h = log(2),  lty = 2, col = "gray")
abline(h = log(5),  lty = 2, col = "gray")
abline(h = log(10), lty = 2, col = "gray")
abline(h = log(20), lty = 2, col = "gray")

lines(1:Tlen, avg_list[[1]], col = cols[2], lwd = 1.6)
lines(1:Tlen, avg_list[[2]], col = cols[3], lwd = 1.6)
lines(1:Tlen, avg_list[[3]], col = cols[4], lwd = 1.6)
lines(1:Tlen, avg_list[[4]], col = cols[5], lwd = 1.6)
lines(1:Tlen, avg_list[[5]], col = cols[6], lwd = 1.6)


legend("topleft",
       legend = c("Baseline", "-5% VaR", "-10% VaR",
                  "-5% ES","-10% ES", "-5% Both"),
       col = cols, lty = 1, lwd = 1.6, bty = "n",
       x.intersp = 0.6,
       y.intersp = 0.7)


## ---------- rejection rate ----------
for (idx in seq_along(param_sets)) {
  k <- param_sets[[idx]]["k"]; l_mult <- param_sets[[idx]]["l_mult"]
  cat(sprintf("k=%.2f, l=%.2f -> Rejection rate over %d sims: %.4f\n",
              k, l_mult, B, rej_list[idx]))
}



## -------------------------------------------------------------------------------- 
## loss and prediction

set.seed(20251007)
l <- 10; n <- 1000; Tlen <- l + n
p <- 0.95
k <- 0.95     
ell <- 1.00 

## data generated process
vals <- c(0, 0.1 * 1:5, -0.1 * 1:5)
eps  <- sample(vals, size = Tlen, replace = TRUE)
z <- 1.64 * k + eps
r <- 2.06 * ell + eps
L <- rnorm(Tlen, 0, 1)

ymin <- min(L, z, r, 1.64, 2.06) - 0.7
ymax <- max(L, z, r, 1.64, 2.06) + 0.5

plot(L, type = "l", col = "black", lwd = 1.2,
     xlab = "number of data (t)", ylab = "Loss and forecast",
     #main = sprintf("Loss and Forecasts (k=%.2f, l=%.2f)", k, ell),
     ylim = c(ymin, ymax))

# predicted VaR/ES
lines(z, col = rgb(1, 0, 0, 0.5),  lwd = 1)
lines(r, col = rgb(0, 0, 1, 0.5), lwd = 1)

# true VaR / ES
abline(h = 1.64, col = "red",  lty = 2, lwd = 1.4)
abline(h = 2.06, col = "blue", lty = 2, lwd = 1.4)
legend("bottomleft",
       legend = c("Loss", "VaR forecast", "ES forecast",
                  "True VaR", "True ES"),
       col = c("black", rgb(1,0,0,0.5), rgb(0,0,1,0.5), "red", "blue"),
       lty = c(1, 1, 1, 2, 2), lwd = c(1.2, 1, 1, 1.4, 1.4),
       bty = "n", cex = 0.85, y.intersp = 0.6, x.intersp = 0.6)

