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

## Load required library for B-spline basis functions
library(splines)

## Read COVID-19 death data from English hospitals
## Data contains daily deaths (nhs) and corresponding julian day (julian) in 2020
dat <- read.table("engcov.txt", header = TRUE)

## Extract relevant variables
y <- dat$nhs ## Deaths in NHS hospitals
t <- dat$julian ## Day of the year 2020
n <- length(y) ## Number of observations

## Calculate probability distribution for infection-to-death duration

## Duration from infection to death follows log-normal distribution
## Based on epidemiological data: log(d) ~ N(3.151, 0.469^2)
## where d is the duration in days from infection to death
edur <- 3.151  ## Mean of log-duration
sdur <- 0.469  ## Standard deviation of log-duration

## The probability of each fatal disease duration from 1 day up to 80 days max.
## We cap at 80 days as longer durations have negligible probability
d <- 1:80
pd <- dlnorm(d, edur, sdur) ## Log-normal probability density for each duration
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
  
  n <- length(t)  ## Number of death observations

  ## Define time range for modeling infections from t[1]-30 to t[n]
  ## We need to model infections that occurred up to 30 days before first death
  t_min <- t[1] - backdays ## Start: 30 days before first death
  t_max <- t[n]            ## End: last observation day

  
  ## Create knot sequence for B-splines
  ## B-splines require k+4 evenly spaced knots, with middle k-2 covering our time range
  ##   - The first and last 2 knots are boundary knots
  ##   - The middle k knots define the basis functions
  ## This ensures the splines span the entire period of interest
  knots <- seq(t_min, t_max, length = k + 4)
  
  ## Create time points where f(t) will be evaluated.
  ## Using by=1 ensures we have integer days only (not fractional days)
  ## This gives us exactly t_max - t_min + 1 points
  t_cover <- seq(t_min, t_max, by=1) ## From t[1]-30 to t[n], covering all possible infection times.
  
  ## Generate B-spline basis matrix X_tilde.
  ## Each column is a B-spline basis function evaluated at t_cover.
  ## Unless outer.ok is true, the values in x must be between the “inner” knots knots[ord] and knots[ length(knots) - (ord-1)].
  X_tilde <- splineDesign(knots, t_cover, outer.ok = TRUE)
  
  ## Initialize model matrix X for deaths
  X <- matrix(0, n, k) ## X[i,j] gives contribution of beta_j to expected deaths on day t[i].
                       ## This accounts for the convolution of infections with the duration distribution
  
  ## Build X by convolving X_tilde with probability distribution pd
  for (i in 1:n) { ## For each observation day i
    ## Determine range of infection days that contribute to deaths on day t[i].
    j_max <- min(29 + i, 80) ## Upper limit: either 29 days before day i, or 80 days before, whichever is less.
    
    ## Loop over possible durations from infection to death
    for (j in 1:j_max) {
      ## Calculate index in t_cover corresponding to infection time
      ## If death occurs on day t[i] and duration is j days, then infection occurred on day t[i] - j
      ## Since t_cover starts at t[1]-backdays, the index is:
      subscript <- 30 + i - j
      
      ## Add weighted contribution of B-splines at infection time
      ## Weighted by probability of duration j
      if (subscript > 0 && subscript <= nrow(X_tilde)) {
        X[i, ] <- X[i, ] + X_tilde[subscript, ] * pd[j]
      }
    }
  }
  
  ## Create penalty matrix S for controlling smoothness.
  ## Penalizes second differences: (beta_{k-1} - 2*beta_k + beta_{k+1})^2.
  ## This encourages smooth changes in the beta coefficients.
  D <- diff(diag(k), differences = 2) ## Second difference matrix
  S <- crossprod(D) ## S = D^T D      ## S = D^T * D, dimension k x k
  
  return(list(X = X, X_tilde = X_tilde, S = S, t_cover = t_cover)) ## Return all matrices as a list
}

## Create the matrices using constructed function.
matrices <- setup_matrices(t, pd, k = 80, backdays = 30)
X <- matrices$X               ## Model matrix for deaths (150 x 80)
X_tilde <- matrices$X_tilde   ## B-spline basis matrix (180 x 80)
S <- matrices$S               ## Penalty matrix (80 x 80)
t_cover <- matrices$t_cover   ## Time points for f(t) (length 180)


nll <- function(gamma, X, y, S, lambda) {
## Compute the negative log-likelihood function for Poisson model
## Working with gamma = log(beta) to ensure beta > 0 (required for rates).
## The model is: y_i ~ Poi(mu_i) where mu_i = sum_j X[i,j] * beta_j
## The log-likelihood is: l = sum_i (y_i * log(mu_i) - mu_i - log(y_i!))
## The penalty is: P = (lambda/2) * beta^T * S * beta

## Arguments:
 ##   gamma: log-transformed parameters (length k vector)
 ##   X: model matrix (n x k)
 ##   y: observed deaths (length n vector)
 ##   S: penalty matrix (k x k)
 ##   lambda: smoothing parameter (scalar > 0)

## Returns 
 ## Scalar value of penalized negative log-likelihood

  
  beta <- exp(gamma) ## Transform to ensure beta > 0
  mu <- X %*% beta  ## Expected deaths: mu = X * beta
  mu <- as.vector(mu) ## Convert from matrix to vector
  ## Poisson log-likelihood: sum(y*log(mu) - mu - log(y!))
  ll <- sum(y * log(mu) - mu) ## We drop log(y!) as it doesn't depend on parameters
  ## Add smoothing penalty: lambda * beta^T * S * beta / 2
  penalty <- (lambda / 2) * sum(beta * (S %*% beta))
  
  return(-ll + penalty) ## Return negative penalized log-likelihood,We minimize this, equivalent to maximizing penalized likelihood.
}


grad_nll <- function(gamma, X, y, S, lambda) {
## Gradient of negative log-likelihood with respect to gamma
## This uses the chain rule: d(nll)/d(gamma) = d(nll)/d(beta) * d(beta)/d(gamma)
  
  ## Arguments:
     ##   gamma: log-transformed parameters (length k vector)
     ##   X: model matrix (n x k)
     ##   y: observed deaths (length n vector)
     ##   S: penalty matrix (k x k)
     ##   lambda: smoothing parameter (scalar > 0)
  
## Returns:
##   Gradient vector (length k)
  
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
## This is essential to verify correctness of derivative calculations.
## Uses central difference approximation: f'(x) ≈ (f(x+h) - f(x)) / h
test_gradient <- function() {
  ## Use random starting values for testing
  gamma_test <- rnorm(ncol(X))
  lambda_test <- 5e-5
  
  ## Compute analytical gradient using the formula
  grad_analytical <- grad_nll(gamma_test, X, y, S, lambda_test)
  
  ## Compute numerical gradient by finite differences
  eps <- 1e-7  ## Small step size for finite difference
  grad_numerical <- numeric(length(gamma_test))

   ## Get function value at test point
  nll_0 <- nll(gamma_test, X, y, S, lambda_test)

   ## Loop over each parameter, computing partial derivative
  for (i in 1:length(gamma_test)) {
    gamma_temp <- gamma_test
    gamma_temp[i] <- gamma_temp[i] + eps  ## Perturb i-th parameter
    nll_1 <- nll(gamma_temp, X, y, S, lambda_test)
    grad_numerical[i] <- (nll_1 - nll_0) / eps ## Finite difference
  }
  
  ## Compare gradients
  cat("Sample gradient comparison (first 5 elements):\n")
  print(head(cbind(Analytical = grad_analytical, 
                   Numerical = grad_numerical,
                   Difference = grad_analytical - grad_numerical), 5))
}

## Run gradient test
cat("\n=== Testing Gradient Function ===\n")
test_gradient()

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
f_fitted <- as.vector(matrices$X_tilde %*% beta_hat) ## f(t), the infection curve evaluated at extended time points
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

compute_bic <- function(gamma, lambda, X, y, S) {
## Function to compute effective degrees of freedom and BIC
## EDF = trace(H_lambda^{-1} * H_0) where:
## H_lambda = X^T W X + lambda*S (Hessian with penalty)
## H_0 = X^T W X (Hessian without penalty)
## W = diag(y_i / mu_i^2)
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  nll_value <- -sum(y * log(mu) - mu) ## Negative log-likelihood (without penalty)
  W <- y / (mu^2) ## Weight matrix W for Hessian computation, W_ii = y_i / mu_i^2
  
  ## Compute Hessian matrices, w.r.t. beta, not gamma.
  XtWX <- t(X) %*% (W * X)  ## X^T W X (using vectorized multiplication)
  H_lambda <- XtWX + lambda * S
  H_0 <- XtWX
  L <- chol(H_lambda) ## strictly positive definite
  H_lambda_inv <- chol2inv(L) ## get H_lambda^{-1}
  edf <- sum(diag(H_lambda_inv %*% H_0)) ## compute trace(H_lambda^{-1} * H_0)
  
  ## BIC = -2*log-likelihood + log(n)*EDF
  bic <- 2 * nll_value + log(length(y)) * edf
  return(list(bic = bic, edf = edf, nll = nll_value))
}

## Grid search over log(lambda)
log_lambda_seq <- seq(-13, -7, length = 50)
lambda_seq <- exp(log_lambda_seq)

## Storage for results
bic_values <- numeric(length(lambda_seq))
edf_values <- numeric(length(lambda_seq))
nll_values <- numeric(length(lambda_seq))
gamma_fits <- matrix(0, length(gamma_init), length(lambda_seq))

for (i in 1:length(lambda_seq)) { ## Loop over lambda values
  lambda_current <- lambda_seq[i]
  
  fit_current <- optim( ## Optimize for current lambda
    par = gamma_init,  ## Using same starting values
    fn = nll, gr = grad_nll, X = X, y = y,
    S = S, lambda = lambda_current, method = "BFGS",control = list(maxit = 1000))

  ## Store results
  gamma_fits[, i] <- fit_current$par 
  bic_result <- compute_bic(fit_current$par, lambda_current, X, y, S)
  bic_values[i] <- bic_result$bic
  edf_values[i] <- bic_result$edf
  nll_values[i] <- bic_result$nll
  
  if (i %% 10 == 0) { ## Progress update every 10 iterations
    cat("  Completed", i, "of", length(lambda_seq), "\n")
  }
}

## Find optimal lambda
optimal_index <- which.min(bic_values)
lambda_optimal <- lambda_seq[optimal_index]
log_lambda_optimal <- log_lambda_seq[optimal_index]
gamma_optimal <- gamma_fits[, optimal_index]

cat("\nOptimal log(lambda):", log_lambda_optimal, "\n")
cat("Optimal lambda:", lambda_optimal, "\n")
cat("Minimum BIC:", bic_values[optimal_index], "\n")
cat("Effective degrees of freedom at optimum:", edf_values[optimal_index], "\n")

## Plot BIC curve
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1)) 
plot(log_lambda_seq, bic_values, type = "l", lwd = 2,
     xlab = "log(lambda)", ylab = "BIC", main = "BIC vs log(lambda)")
points(log_lambda_optimal, bic_values[optimal_index], 
       col = "red", pch = 19, cex = 1.5)
abline(v = log_lambda_optimal, lty = 2, col = "red")

plot(edf_values, bic_values, type = "l", lwd = 2,
     xlab = "Effective Degrees of Freedom", ylab = "BIC", main = "BIC vs EDF")
points(edf_values[optimal_index], bic_values[optimal_index],
       col = "red", pch = 19, cex = 1.5)

## obtain 200 bootstrap replicates of f_hat
n_bootstrap <- 200

f_bootstrap <- matrix(0, length(t_cover), n_bootstrap) ## Storage for bootstrap results
                                                       ## Each column is f(t) for one bootstrap replicate


for (b in 1:n_bootstrap) { ## Bootstrap loop
  wb <- tabulate(sample(n, replace = TRUE), n) ## Generate bootstrap weights, equivalent to resampling data with replacement
  nll_bootstrap <- function(gamma, X, y, S, lambda, weights) {
    ## Modified negative log-likelihood with bootstrap weights
    beta <- exp(gamma)
    mu <- X %*% beta
    mu <- as.vector(mu)
    ll <- sum(weights * (y * log(mu) - mu)) ## Weighted log-likelihood
    penalty <- (lambda / 2) * sum(beta * (S %*% beta))
    return(-ll + penalty)
  }
  
  grad_nll_bootstrap <- function(gamma, X, y, S, lambda, weights) {
    ## Modified gradient with bootstrap weights
    beta <- exp(gamma)
    mu <- X %*% beta
    mu <- as.vector(mu)
    grad_ll_beta <- t(X) %*% (weights * (y / mu - 1))
    grad_penalty_beta <- lambda * (S %*% beta)
    grad <- -beta * (grad_ll_beta - grad_penalty_beta)
    return(as.vector(grad))
  }
  
  fit_boot <- optim( ## Optimize for bootstrap sample
    par = gamma_optimal,  ## Start from optimal parameters
    fn = nll_bootstrap, gr = grad_nll_bootstrap,
    X = X, y = y, S = S, lambda = lambda_optimal,
    weights = wb, method = "BFGS", control = list(maxit = 1000)
  )
  
  beta_boot <- exp(fit_boot$par) ## Calculate f(t) for this bootstrap replicate
  f_bootstrap[, b] <- as.vector(matrices$X_tilde %*% beta_boot)
  
  if (b %% 20 == 0) { ## Progress update
    cat("  Completed", b, "of", n_bootstrap, "bootstrap replicates\n")
  }
}

## Calculate 95% confidence intervals
f_lower <- apply(f_bootstrap, 1, quantile, probs = 0.025) ## 2.5th percentiles for each time point
f_upper <- apply(f_bootstrap, 1, quantile, probs = 0.975) ## 97.5th percentiles for each time point

f_mean <- apply(f_bootstrap, 1, mean) ## calculate mean bootstrap estimate

## Calculate final fitted values with optimal lambda
beta_optimal <- exp(gamma_optimal)
mu_optimal <- as.vector(X %*% beta_optimal)
f_optimal <- as.vector(matrices$X_tilde %*% beta_optimal)

par(mfrow = c(2, 1), mar = c(4, 4, 3, 1)) ## Create comprehensive plot with two panels

## Panel 1: Deaths - Actual vs Fitted
plot(t, y, type = "p", pch = 19, cex = 0.8, col = "black",
     xlab = "Day of Year 2020", ylab = "Number of Deaths",
     main = "COVID-19 Deaths: Actual vs Model Fit", ylim = c(0, max(y) * 1.1))
lines(t, mu_optimal, col = "red", lwd = 2)
legend("topright", legend = c("Observed Deaths", "Model Fit"),
       col = c("black", "red"), pch = c(19, NA),
       lty = c(NA, 1), lwd = c(NA, 2), bty = "n")

## Add text with model information
text(min(t) + 5, max(y) * 1.05,
     paste0("Optimal log(λ) = ", round(log_lambda_optimal, 2)),
     pos = 4, cex = 0.8)

## Panel 2: Inferred Daily Infections with 95% CI
plot(t_cover, f_optimal, type = "l", col = "blue", lwd = 2,
     xlab = "Day of Year 2020", ylab = "Daily New Infections",
     main = "Inferred Daily Infection Rate with 95% Confidence Intervals",
     ylim = c(0, max(f_upper) * 1.1))

## Add confidence interval as shaded region
polygon(c(t_cover, rev(t_cover)), c(f_lower, rev(f_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)

## Add confidence interval boundary lines
lines(t_cover, f_lower, col = "blue", lty = 2, lwd = 1)
lines(t_cover, f_upper, col = "blue", lty = 2, lwd = 1)
lines(t_cover, f_optimal, col = "blue", lwd = 2)
abline(v = min(t), lty = 2, col = "gray") ## Mark first observation day
legend("topright",
       legend = c("Estimated f(t)", "95% CI", "First Death"),
       col = c("blue", "blue", "gray"), lty = c(1, 2, 2), 
       lwd = c(2, 1, 1), bty = "n")












