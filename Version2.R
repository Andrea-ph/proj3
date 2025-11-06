################################################################################
## Practical 3: Smooth Deconvolution - COVID-19 Deaths Analysis
## GitHub Repository: [https://github.com/Andrea-ph/proj3]
###########################################################
####proj3 - Group 33 - Extended Statistical Programming ###
#### Group members as below ################################
#### Shuo Li (s2795688), Zhe Zhu (s2841606), Antrea Filippou (s2766374)
#### Contributions as below ################################
#### Shuo Li:  (%) #########################
#### Zhe Zhu: (%) ######################
#### Antrea Filippou: (%) ###################
############################################################
## This code implements a smooth deconvolution model to infer daily new COVID-19
## infections from observed death data, using B-splines with a smoothing penalty.
################################################################################

library(splines)
dat <- read.table("engcov.txt", header = TRUE)
## Extract relevant variables
y <- dat$nhs ## Deaths in NHS hospitals
t <- dat$julian ## Day of the year 2020
n <- length(y) ## Number of observations

## Calculate probability distribution for infection-to-death duration

## Duration from infection to death follows log-normal distribution
## log(d) ~ N(3.151, 0.469^2) where d is duration in days
edur <- 3.151  ## Mean of log-duration
sdur <- 0.469  ## Standard deviation of log-duration

## The probability of each fatal disease duration from 1 day up tp the max 80 days
d <- 1:80
pd <- dlnorm(d, edur, sdur) ## Log-normal density for each duration
pd <- pd / sum(pd) ## Normalization: sum(pd)=1

## Function to create model matrices X, X_tilde, and penalty matrix S.
## X_tilde is the B-spline basis for f(t) evaluated at extended time points.
## X is the model matrix for deaths incorporating the convolution structure.
## S is the penalty matrix for smoothing.

setup_matrices <- function(t, pd, k = 80, backdays = 30) {
## t: a vector of observation days (julian days).
## pd: probability distribution of infection-to-death duration.
## k: number of B-spline basis functions.
## backdays: 30 to avoid looking too far back for the origins of the infections causing the first deaths.
  
  n <- length(t)  ## Number of observations
  
  ## Define the range for f(t): from t[1]-30 to t[n]
  ## We need to model infections that occurred up to 30 days before first death
  t_min <- t[1] - backdays
  t_max <- t[n]
  
  ## Create knot sequence for B-splines
  ## Need k+4 evenly spaced knots, with middle k-2 covering our time range
  ## This ensures the splines span the entire period of interest
  knots <- seq(t_min, t_max, length = k + 4)
  
  ## Create time points where f(t) will be evaluated.
  t_cover <- seq(t_min, t_max, by=1) ## From t[1]-30 to t[n], covering all possible infection times (integers).
  
  ## Generate B-spline basis matrix X_tilde.
  ## Each column is a B-spline basis function evaluated at t_cover.
  ## Unless outer.ok is true, the values in x must be between the “inner” knots knots[ord] and knots[ length(knots) - (ord-1)].
  X_tilde <- splineDesign(knots, t_cover, outer.ok = TRUE)
  
  ## Initialize model matrix X for deaths
  X <- matrix(0, n, k) ## X[i,j] gives contribution of beta_j to expected deaths on day t[i].
  
  ## Build X by convolving X_tilde with probability distribution pd
  for (i in 1:n) { ## For each observation day i
    ## Determine range of infection days that contribute to deaths on day t[i].
    j_max <- min(29 + i, 80) ## Upper limit: either 29 days before day i, or 80 days before, whichever is less.
    
    ## For each possible duration j from infection to death
    for (j in 1:j_max) {
      ## Index in extended time vector for infection day
      ## Day t[i] - j is when infection occurred that leads to death on day t[i]
      subscript <- 30 + i - j
      
      ## Add weighted contribution of B-splines at infection time
      ## Weighted by probability of duration j
      if (subscript > 0 && subscript <= nrow(X_tilde)) {
        X[i, ] <- X[i, ] + X_tilde[subscript, ] * pd[j]
      }
    }
  }
  
  ## Create penalty matrix S.
  ## Penalizes second differences: (beta_{k-1} - 2*beta_k + beta_{k+1})^2.
  ## This encourages smooth changes in the beta coefficients.
  D <- diff(diag(k), differences = 2) ## Second difference matrix
  S <- crossprod(D) ## S = D^T D
  
  return(list(X = X, X_tilde = X_tilde, S = S, t_cover = t_cover)) ## Return all matrices as a list
}

## Create the matrices
matrices <- setup_matrices(t, pd, k = 80, backdays = 30)
X <- matrices$X
X_tilde <- matrices$X_tilde
S <- matrices$S
t_cover <- matrices$t_cover

## Negative log-likelihood function for Poisson model
## Working with gamma = log(beta) to ensure beta > 0
## Returns penalized negative log-likelihood

nll <- function(gamma, X, y, S, lambda) {
## gamma: log-transformed parameters (ensures positivity)
## X: model matrix
## y: observed deaths
## S: penalty matrix
## lambda: smoothing parameter
  
  beta <- exp(gamma) ## Transform to ensure beta > 0
  mu <- X %*% beta  ## Expected deaths: mu = X * beta
  mu <- as.vector(mu) ## Convert to vector
  ## Poisson log-likelihood: sum(y*log(mu) - mu - log(y!))
  ll <- sum(y * log(mu) - mu) ## We drop log(y!) as it doesn't depend on parameters
  ## Add smoothing penalty: lambda * beta^T * S * beta / 2
  penalty <- (lambda / 2) * sum(beta * (S %*% beta))
  
  return(-ll + penalty) ## Return negative penalized log-likelihood
}

## Gradient of negative log-likelihood with respect to gamma
## This uses the chain rule: d(nll)/d(gamma) = d(nll)/d(beta) * d(beta)/d(gamma)
grad_nll <- function(gamma, X, y, S, lambda) {
  
  beta <- exp(gamma) ## Transform parameters
  mu <- X %*% beta ## Expected deaths
  mu <- as.vector(mu) ## Convert to vector
  
  ## Gradient of log-likelihood w.r.t. beta
  ## d(ll)/d(beta) = X^T * (y/mu - 1)
  grad_ll_beta <- t(X) %*% (y / mu - 1)
  
  ## Gradient of penalty w.r.t. beta
  ## d(penalty)/d(beta) = lambda * S * beta
  grad_penalty_beta <- lambda * (S %*% beta)
  
  ## Using Chain rule: d(beta)/d(gamma) = diag(beta)
  ## So d(nll)/d(gamma) = -diag(beta) * (grad_ll_beta - grad_penalty_beta)
  grad_gamma <- -beta * (grad_ll_beta - grad_penalty_beta)
  return(as.vector(grad_gamma))
}

## Function to test gradient accuracy
test_gradient <- function() {
  ## Use random starting values for testing
  gamma_test <- rnorm(ncol(X))
  lambda_test <- 5e-5
  
  ## Compute analytical gradient
  grad_analytical <- grad_nll(gamma_test, X, y, S, lambda_test)
  ## Compute numerical gradient by finite differences
  delta <- 1e-6  # Small step for finite difference
  grad_numerical <- numeric(length(gamma_test))
  
  for (i in 1:length(gamma_test)) { 
    gamma_plus <- gamma_test 
    gamma_plus[i] <- gamma_plus[i] + delta 
    nll_plus <- nll(gamma_plus, X, y, S, lambda_test)
    gamma_minus <- gamma_test
    gamma_minus[i] <- gamma_minus[i] - delta
    nll_minus <- nll(gamma_minus, X, y, S, lambda_test)
    
    grad_numerical[i] <- (nll_plus - nll_minus) / (2 * delta) ## Using centered difference
  }
  
  ## Compare gradients
  max_diff <- max(abs(grad_analytical - grad_numerical))
  cat("Maximum difference between analytical and numerical gradient:", max_diff, "\n")
  
  if (max_diff < 1e-5) { ## Test passes if difference is very small (< 1e-5)
    cat("Gradient test PASSED!\n")
  } else {
    cat("Gradient test FAILED!\n")
    cat("Sample differences (first 10 elements):\n")
    print(head(cbind(Analytical = grad_analytical, Numerical = grad_numerical, 
                     Difference = grad_analytical - grad_numerical), 10))
  }
}

test_gradient() ## Run gradient test

cat("\n=== Sanity Check: Fitting model with lambda = 5e-5 ===\n")
lambda_initial <- 5e-5 ## Set smoothing parameter for initial fit

gamma_init <- rnorm(ncol(X), mean = 0, sd = 0.1) 
## Starting values: initialize gamma to small random values. This corresponds to beta ≈ 1

fit_initial <- optim( ## Optimize using BFGS method from optim.
  par = gamma_init, fn = nll, gr = grad_nll, X = X, y = y, S = S,
  lambda = lambda_initial, method = "BFGS", control = list(maxit = 1000)
) ## BFGS is a quasi-Newton method that builds up Hessian approximation

if (fit_initial$convergence != 0) { ## Check convergence
  warning("Optimization did not converge!")
} else {
  cat("Optimization converged successfully\n")
}

gamma_hat <- fit_initial$par ## Extract fitted parameters
beta_hat <- exp(gamma_hat)
mu_fitted <- as.vector(X %*% beta_hat) ## Calculate fitted values
f_fitted <- as.vector(matrices$X_tilde %*% beta_hat) ## f(t), the infection curve evaluated at extended time points.

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)) ## Plot results

## Plot 1: Actual vs Fitted Deaths
plot(t, y, type = "p", pch = 19, col = "black",
     xlab = "Day of Year 2020", ylab = "Deaths", main = "Actual vs Fitted Deaths")
lines(t, mu_fitted, col = "red", lwd = 2)
legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"),
       pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2))

## Plot 2: Inferred Daily Infections
## Please note: t_cover includes days before first observation (from t[1]-30)
plot(t_cover, f_fitted, type = "l", col = "blue", lwd = 2,
     xlab = "Day of Year 2020", ylab = "Daily New Infections", main = "Inferred Infection Curve f(t)")
abline(v = min(t), lty = 2, col = "gray") ## Mark first observation day
text(min(t) + 5, max(f_fitted) * 0.9, "First death", pos = 4, col = "gray")
cat("\nSanity check plots created\n")







