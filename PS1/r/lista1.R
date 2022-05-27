options(scipen = 999)
library(magrittr)

# Exercise 1 --------------------------------------------------------------

## Calibration
rho = 0.95
sig = 0.007

N = 9        # Number of points (states);
m = 3        # Scaling parameter (I don't know why!); 

T = 10000

## Grid and Transition matrix function of Tauchen's Method

tauchen = function(N, sig, rho, m) {
  
  # Grid
  grid = rep(0, N)
  grid[1] = - m * sqrt(sig^2/(1-rho^2))
  grid[N] = + m * sqrt(sig^2/(1-rho^2))
  step = (grid[N] - grid[1]) / (N - 1)
  for (i in 2:(N-1)) {
    grid[i] = grid[i-1] + step
  }
  
  # Transition matrix
  if (N > 1) {
    step = grid[2] - grid[1]
    P = array(0, dim = c(N, N))
    for (j in 1:N) {
      for (k in 1:N) {
        if (k == 1) {
          P[j, k] = pnorm((grid[k] - rho * grid[j] + (step / 2)) / sig)
        }
        else if (k == N) {
          P[j, k] = 1 - pnorm((grid[k] - rho * grid[j] - (step / 2)) / sig)
        }
        else {
          P[j, k] = pnorm((grid[k] - rho * grid[j] + (step / 2)) / sig) - pnorm((grid[k] - rho * grid[j] - (step / 2)) / sig)
        }
      }
    }
  } else {
    P = 1
  }
  return(list(zgrid = round(grid, 4), P = round(P,4)))
}

## Grid and Transition matrix (Tauchen's Method)
Tauchen95 = tauchen(N, sig, rho, m)
Tauchen95

# Exercise 2 --------------------------------------------------------------

## Grid and Transition matrix function of Rouwenhorst's Method

rouwen <- function(rho, sigma, n){
  
  zgrid <- seq(
    from = - (sig / sqrt(1-rho^2)) * sqrt(N-1),
    to = + (sig / sqrt(1-rho^2)) * sqrt(N-1),
    length = n
  )

  p  <- (rho+1)/2
  nu <- ((n-1)/(1-rho^2))^(1/2) * sigma
  P  <- matrix(c(p, 1 - p, 1 - p, p), nrow = 2, ncol = 2)
  
  for (i in 3:n){
    zeros <- matrix(0, nrow = i-1, ncol = 1)
    zzeros <- matrix(0, nrow = 1, ncol = i-1)
    
    P        <- p * rbind(cbind(P, zeros), cbind(zzeros, 0)) + 
               (1-p) * rbind(cbind(zeros, P), cbind(0, zzeros)) + 
               (1-p) * rbind(cbind(zzeros, 0), cbind(P, zeros)) + 
                p * rbind(cbind(0, zzeros), cbind(zeros, P))
  }
  
  for (j in 1:N){
    P[j,] = P[j,]/sum(P[j,])
  }
  
  return(list(zgrid = round(zgrid, 4), P = round(P, 4)))
}

## Grid and Transition matrix (Rouwenhorst's Method)
Rouwen95 <- rouwen(rho, sig, N)
Rouwen95

# Exercise 3 --------------------------------------------------------------

## Function of a AR(1) continuous process 

ar1_simulation <- function(rho, sig, T){
  
  Z = rep(0, T); Z[1] = 0
  
  eps = rnorm(T, mean = 0, sd = sig)
  
  for(i in 2:T){
    Z[i] = rho*Z[i-1] + eps[i]
  }
  
  return(list(Z = round(Z, 4), eps = round(eps, 4)))
}

## Simulation of the process 
## (rho = 0.95, mu = 0, sig = 0.007, T = 10000, Z[1] = 0)
set.seed(1347)
ar1_95 <- ar1_simulation(rho, sig, T)

## Show the first 10 numbers
ar1_95 %>% purrr::map(.f = ~ head(.x, 10))

## Function for discretize the process
discret <- function(th0, sig, eps, P, T){ 
  
  idx = rep(1, T)
  idx[1] = th0
  cum = t(apply(P, 1, cumsum))
  
  for(i in 2:T){ 
    x = which(pnorm(eps[i], mean = 0, sd = sig) <= cum[idx[i-1],])
    idx[i] = x[1]
  }
  
  return(idx)
}

## Tauchen's Method, rho = 0.95 ####

## Defining the initial state
th0 <- which(Tauchen95$zgrid == median(Tauchen95$zgrid))

## Returning the indices of the grid
idx = discret(th0, sig, ar1_95$eps, Tauchen95$P, T)

## Simulation the discretized process
ztauchen95 <- Tauchen95$zgrid[idx]

## Plotting
par(mfrow=c(2,1))

plot(ar1_95$Z, type = 'S', col = 1, 
     main = "Realização do Processo Contínuo para Z(t), rho = 0.95", 
     xlab = "Período de Tempo", ylab = "Realização de Z(t)", 
     ylim = c(-0.08, 0.08))
lines(ztauchen95, col = 2)
abline(h = Tauchen95$zgrid, col = scales::alpha('black', 0.1), lty = 2)
legend("topright",
  legend = c("Realização do Processo Contínuo para Z(t)", 
             "Processo Z(t) discretizado via Método de Tauchen"),
  col = c("black", "red"),
  lty = c(1, 1)
)

## Rouwenhorst's Method, rho = 0.95 ####

## Defining the initial state
th0 <- which(Rouwen95$zgrid == median(Rouwen95$zgrid))

## Returning the indices of the grid
idx = discret(th0, sig, ar1_95$eps, Rouwen95$P, T)

## Simulation the discretized process
zrouwen95 <- Rouwen95$zgrid[idx]

## Plotting
plot(ar1_95$Z, type = 'S', col = 1, 
     main = "Realização do Processo Contínuo para Z(t), rho = 0.95", 
     xlab = "Período de Tempo", ylab = "Realização de Z(t)", 
     ylim = c(-0.08, 0.08))
lines(zrouwen95, col = 2)
abline(h = Rouwen95$zgrid, col = scales::alpha('black', 0.1), lty = 2)
legend("topright",
       legend = c("Realização do Processo Contínuo para Z(t)", 
                  "Processo Z(t) discretizado via Método de Rouwenhorst"),
       col = c("black", "red"),
       lty = c(1, 1), text.font = 1
)

# Exercise 4 --------------------------------------------------------------

## Compute the lag (Tauchen's method)
ztauchen95_lag1 <- dplyr::lag(ztauchen95, 1)

## Run the regression and show the results
lm_tauchen95 <- lm(ztauchen95 ~ 0 + ztauchen95_lag1)
broom::tidy(lm_tauchen95)

## Hyphotesis test (H0: rho = 0.95 vs H1: rho ≠ 0.95)
tauchen95_tstat <- (lm_tauchen95$coefficients - 0.95)/sqrt(diag(vcov(lm_tauchen95)))
tauchen95_tstat

## Confidence interval (if inside, we do not reject the null)
qt(c(.025, .975), df = lm_tauchen95$df.residual)

## Compute the lag (Tauchen's method)
zrouwen95_lag1 <- dplyr::lag(zrouwen95, 1)

## Run the regression and show the results
lm_rouwen95 <- lm(zrouwen95 ~ 0 + zrouwen95_lag1)
broom::tidy(lm_rouwen95)

## Hyphotesis test (H0: rho = 0.95 vs H1: rho ≠ 0.95)
rouwen95_tstat <- (lm_rouwen95$coefficients - 0.95)/sqrt(diag(vcov(lm_rouwen95)))
rouwen95_tstat

## Confidence interval (if inside, we do not reject the null)
qt(c(.025, .975), df = lm_rouwen95$df.residual)

# Exercise 5 --------------------------------------------------------------

## Recalibration
rho = 0.99
sig = 0.007

N = 9        # Number of points (states);
m = 3        # Scaling parameter (I don't know why!); 

T = 10000

## 5.1 Tauchen's Method ####
Tauchen99 <- tauchen(N, sig, rho, m)

## 5.2 Rouwenhorst's Method ####
Rouwen99 <- rouwen(rho, sig, N)

## 5.3 Process simulation and discretize #### 

## Simulation of the process 
## (rho = 0.99, mu = 0, sig = 0.007, T = 10000, Z[1] = 0)
set.seed(1347)
ar1_99 <- ar1_simulation(rho, sig, T)

## Show the first 10 numbers
ar1_99 %>% purrr::map(.f = ~ head(.x, 10))

## Tauchen's Method, rho = 0.99 ####

## Defining the initial state
th0 <- which(Tauchen99$zgrid == median(Tauchen99$zgrid))

## Returning the indices of the grid
idx = discret(th0, sig, ar1_99$eps, Tauchen99$P, T)

## Simulation the discretized process
ztauchen99 <- Tauchen99$zgrid[idx]

## Plotting
par(mfrow = c(2, 1))

plot(ar1_99$Z, type = 'S', col = 1, 
     main = "Realização do Processo Contínuo para Z(t), rho = 0.99", 
     xlab = "Período de Tempo", ylab = "Realização de Z(t)", 
     ylim = c(-0.15, 0.15))
lines(ztauchen99, col = 2)
abline(h = Tauchen99$zgrid, col = scales::alpha('black', 0.1), lty = 2)
legend("topright",
       legend = c("Realização do Processo Contínuo para Z(t)", 
                  "Processo Z(t) discretizado via Método de Tauchen"),
       col = c("black", "red"),
       lty = c(1, 1)
)

## Rouwenhorst's Method, rho = 0.99 ####

## Defining the initial state
th0 <- which(Rouwen99$zgrid == median(Rouwen99$zgrid))

## Returning the indices of the grid
idx = discret(th0, sig, ar1_99$eps, Rouwen99$P, T)

## Simulation the discretized process
zrouwen99 <- Rouwen99$zgrid[idx]

## Plotting
plot(ar1_99$Z, type = 'S', col = 1, 
     main = "Realização do Processo Contínuo para Z(t), rho = 0.99", 
     xlab = "Período de Tempo", ylab = "Realização de Z(t)", 
     ylim = c(-0.15, 0.15))
lines(zrouwen99, col = 2)
abline(h = Rouwen99$zgrid, col = scales::alpha('black', 0.1), lty = 2)
legend("topright",
       legend = c("Realização do Processo Contínuo para Z(t)", 
                  "Processo Z(t) discretizado via Método de Rouwenhorst"),
       col = c("black", "red"),
       lty = c(1, 1), text.font = 1)

## 5.4 Regressions #### 

## Compute the lag (Tauchen's method)
ztauchen99_lag1 <- dplyr::lag(ztauchen99, 1)

## Run the regression and show the results
lm_tauchen99 <- lm(ztauchen99 ~ 0 + ztauchen99_lag1)
broom::tidy(lm_tauchen99)

## Hyphotesis test (H0: rho = 0.99 vs H1: rho ≠ 0.99)
tauchen99_tstat <- (lm_tauchen99$coefficients - 0.99)/sqrt(diag(vcov(lm_tauchen99)))
tauchen99_tstat %>% as.vector()

## Confidence interval (if inside, we do not reject the null)
qt(c(.025, .975), df = lm_tauchen99$df.residual)

## Compute the lag (Tauchen's method)
zrouwen99_lag1 <- dplyr::lag(zrouwen99, 1)

## Run the regression and show the results
lm_rouwen99 <- lm(zrouwen99 ~ 0 + zrouwen99_lag1)
broom::tidy(lm_rouwen99)

## Hyphotesis test (H0: rho = 0.99 vs H1: rho ≠ 0.99)
rouwen99_tstat <- (lm_rouwen99$coefficients - 0.99)/sqrt(diag(vcov(lm_rouwen99)))
rouwen99_tstat %>% as.vector()

## Confidence interval (if inside, we do not reject the null)
qt(c(.025, .975), df = lm_rouwen99$df.residual)
