################################################################################
# Practical 3: Smooth Deconvolution - COVID-19 Deaths Analysis
# GitHub Repository: [https://github.com/Andrea-ph/proj3]
#
###########################################################
####proj3 - Group 33 - Extended Statistical Programming ###
#### Group members as below ################################
#### Shuo Li (s2795688), Zhe Zhu (s2841606), Antrea Filippou (s2766374)
#### Contributions as below ################################
#### Shuo Li:  (%) #########################
#### Zhe Zhu: (%) ######################
#### Antrea Filippou: (%) ###################
############################################################
# This code implements a smooth deconvolution model to infer daily new COVID-19
# infections from observed death data, using B-splines with a smoothing penalty.
################################################################################

# Load required library
library(splines)

#working directory 
 setwd("/Users/andreaphilippou/Desktop/Msc Statistics with Data Science/extended statistical programming")

################################################################################
# STEP 1: Load and prepare data
################################################################################

# Read the COVID-19 death data
# Data contains deaths in English hospitals by day of year 2020
dat <- read.table("engcov (1).txt", header = TRUE)

# Extract relevant variables
y <- dat$nhs      # Deaths in NHS hospitals
t <- dat$julian   # Day of year
n <- length(y)    # Number of observations

################################################################################
# STEP 2: Calculate probability distribution for infection-to-death duration
################################################################################

# Duration from infection to death follows log-normal distribution
# log(d) ~ N(3.151, 0.469^2) where d is duration in days
edur <- 3.151  # Mean of log-duration
sdur <- 0.469  # SD of log-duration

# Calculate probability for each duration from 1 to 80 days
d <- 1:80
pd <- dlnorm(d, edur, sdur)  # Log-normal density for each duration
pd <- pd / sum(pd)           # Normalize to sum to 1 (probability distribution)

################################################################################
# STEP 3: Setup model matrices and penalty matrix
################################################################################

# Function to create model matrices X, X_tilde, and penalty matrix S
# X_tilde is the B-spline basis for f(t) evaluated at extended time points
# X is the model matrix for deaths incorporating the convolution structure
# S is the penalty matrix for smoothing
setup_matrices <- function(t, pd, k = 80) {
  # Parameters:
  # t: vector of observation days (julian days)
  # pd: probability distribution of infection-to-death duration
  # k: number of B-spline basis functions
  
  n <- length(t)  # Number of observations
  
  # Define the range for f(t): from t[1]-30 to t[n]
  # We need to model infections that occurred up to 30 days before first death
  t_min <- t[1] - 30
  t_max <- t[n]
  
  # Create knot sequence for B-splines
  # Need k+4 evenly spaced knots, with middle k-2 covering our time range
  # This ensures the splines span the entire period of interest
  knots <- seq(t_min, t_max, length = k + 4)
  
  # Create time points where f(t) will be evaluated
  # From t[1]-30 to t[n], covering all possible infection times
  t_ext <- seq(t_min, t_max, length = k + 30)
  
  # Generate B-spline basis matrix X_tilde
  # Each column is a B-spline basis function evaluated at t_ext
  X_tilde <- splineDesign(knots, t_ext, outer.ok = TRUE)
  
  # Initialize model matrix X for deaths
  # X[i,j] gives contribution of beta_j to expected deaths on day t[i]
  X <- matrix(0, n, k)
  
  # Build X by convolving X_tilde with probability distribution pd
  # For each observation day i
  for (i in 1:n) {
    # Determine range of infection days that contribute to deaths on day t[i]
    # Upper limit: either 29 days before day i, or 80 days before, whichever is less
    j_max <- min(29 + i, 80)
    
    # For each possible duration j from infection to death
    for (j in 1:j_max) {
      # Index in extended time vector for infection day
      # Day t[i] - j is when infection occurred that leads to death on day t[i]
      idx <- 30 + i - j
      
      # Add weighted contribution of B-splines at infection time
      # Weighted by probability of duration j
      if (idx > 0 && idx <= nrow(X_tilde)) {
        X[i, ] <- X[i, ] + X_tilde[idx, ] * pd[j]
      }
    }
  }
  
  # Create penalty matrix S
  # Penalizes second differences: (beta_{k-1} - 2*beta_k + beta_{k+1})^2
  # This encourages smooth changes in the beta coefficients
  D <- diff(diag(k), differences = 2)  # Second difference matrix
  S <- crossprod(D)                    # S = D^T D
  
  # Return all matrices as a list
  return(list(X = X, X_tilde = X_tilde, S = S, t_ext = t_ext))
}

# Create the matrices
matrices <- setup_matrices(t, pd, k = 80)
X <- matrices$X
S <- matrices$S
t_ext <- matrices$t_ext

################################################################################
# STEP 4: Negative log-likelihood and gradient functions
################################################################################

# Negative log-likelihood function for Poisson model
# Working with gamma = log(beta) to ensure beta > 0
# Returns penalized negative log-likelihood
nll <- function(gamma, X, y, S, lambda) {
  # Parameters:
  # gamma: log-transformed parameters (ensures positivity)
  # X: model matrix
  # y: observed deaths
  # S: penalty matrix
  # lambda: smoothing parameter
  
  beta <- exp(gamma)           # Transform to ensure beta > 0
  mu <- X %*% beta             # Expected deaths: mu = X * beta
  mu <- as.vector(mu)          # Convert to vector
  
  # Poisson log-likelihood: sum(y*log(mu) - mu - log(y!))
  # We drop log(y!) as it doesn't depend on parameters
  ll <- sum(y * log(mu) - mu)
  
  # Add smoothing penalty: lambda * beta^T * S * beta / 2
  penalty <- (lambda / 2) * sum(beta * (S %*% beta))
  
  # Return negative penalized log-likelihood
  return(-ll + penalty)
}

# Gradient of negative log-likelihood with respect to gamma
# This uses the chain rule: d(nll)/d(gamma) = d(nll)/d(beta) * d(beta)/d(gamma)
grad_nll <- function(gamma, X, y, S, lambda) {
  # Parameters: same as nll()
  
  beta <- exp(gamma)           # Transform parameters
  mu <- X %*% beta             # Expected deaths
  mu <- as.vector(mu)          # Convert to vector
  
  # Gradient of log-likelihood w.r.t. beta
  # d(ll)/d(beta) = X^T * (y/mu - 1)
  grad_ll_beta <- t(X) %*% (y / mu - 1)
  
  # Gradient of penalty w.r.t. beta
  # d(penalty)/d(beta) = lambda * S * beta
  grad_penalty_beta <- lambda * (S %*% beta)
  
  # Chain rule: d(beta)/d(gamma) = diag(beta)
  # So d(nll)/d(gamma) = -diag(beta) * (grad_ll_beta - grad_penalty_beta)
  grad <- -beta * (grad_ll_beta - grad_penalty_beta)
  
  return(as.vector(grad))
}

################################################################################
# STEP 5: Test gradient function by finite differencing
################################################################################

# Function to test gradient accuracy
test_gradient <- function() {
  # Use random starting values for testing
  set.seed(123)
  gamma_test <- rnorm(ncol(X))
  lambda_test <- 5e-5
  
  # Compute analytical gradient
  grad_analytical <- grad_nll(gamma_test, X, y, S, lambda_test)
  
  # Compute numerical gradient by finite differences
  eps <- 1e-7  # Small step for finite difference
  grad_numerical <- numeric(length(gamma_test))
  
  nll_0 <- nll(gamma_test, X, y, S, lambda_test)
  
  for (i in 1:length(gamma_test)) {
    gamma_temp <- gamma_test
    gamma_temp[i] <- gamma_temp[i] + eps
    nll_1 <- nll(gamma_temp, X, y, S, lambda_test)
    grad_numerical[i] <- (nll_1 - nll_0) / eps
  }
  
  # Compare gradients
  max_diff <- max(abs(grad_analytical - grad_numerical))
  cat("Maximum difference between analytical and numerical gradient:", max_diff, "\n")
  
  # Test passes if difference is very small (< 1e-5)
  if (max_diff < 1e-5) {
    cat("Gradient test PASSED!\n")
  } else {
    cat("Gradient test FAILED!\n")
    cat("Sample differences (first 10 elements):\n")
    print(head(cbind(Analytical = grad_analytical, 
                     Numerical = grad_numerical,
                     Difference = grad_analytical - grad_numerical), 10))
  }
}

# Run gradient test
cat("\n=== Testing Gradient Function ===\n")
test_gradient()

################################################################################
# STEP 6: Sanity check with lambda = 5e-5
################################################################################

cat("\n=== Sanity Check: Fitting model with lambda = 5e-5 ===\n")

# Set smoothing parameter for initial fit
lambda_initial <- 5e-5

# Starting values: initialize gamma to small random values
# This corresponds to beta ≈ 1
set.seed(123)
gamma_init <- rnorm(ncol(X), mean = 0, sd = 0.1)

# Optimize using BFGS method from optim
# BFGS is a quasi-Newton method that builds up Hessian approximation
fit_initial <- optim(
  par = gamma_init,
  fn = nll,
  gr = grad_nll,
  X = X,
  y = y,
  S = S,
  lambda = lambda_initial,
  method = "BFGS",
  control = list(maxit = 1000)
)

# Check convergence
if (fit_initial$convergence != 0) {
  warning("Optimization did not converge!")
} else {
  cat("Optimization converged successfully\n")
}

# Extract fitted parameters
gamma_hat <- fit_initial$par
beta_hat <- exp(gamma_hat)

# Calculate fitted values
mu_fitted <- as.vector(X %*% beta_hat)

# Calculate f(t) - the infection curve
# f(t) is evaluated at extended time points
f_fitted <- as.vector(matrices$X_tilde %*% beta_hat)

# Plot results
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

# Plot 1: Actual vs Fitted Deaths
plot(t, y, 
     type = "p", 
     pch = 19, 
     col = "black",
     xlab = "Day of Year 2020", 
     ylab = "Deaths",
     main = "Actual vs Fitted Deaths")
lines(t, mu_fitted, col = "red", lwd = 2)
legend("topright", 
       legend = c("Actual", "Fitted"), 
       col = c("black", "red"),
       pch = c(19, NA),
       lty = c(NA, 1),
       lwd = c(NA, 2))

# Plot 2: Inferred Daily Infections
# Note: t_ext includes days before first observation (from t[1]-30)
plot(t_ext, f_fitted,
     type = "l",
     col = "blue",
     lwd = 2,
     xlab = "Day of Year 2020",
     ylab = "Daily New Infections",
     main = "Inferred Infection Curve f(t)")
abline(v = min(t), lty = 2, col = "gray")  # Mark first observation day
text(min(t) + 5, max(f_fitted) * 0.9, 
     "First death", pos = 4, col = "gray")

cat("\nSanity check plots created\n")

################################################################################
# STEP 7: Optimize lambda using BIC
################################################################################

cat("\n=== Optimizing smoothing parameter lambda using BIC ===\n")

# Function to compute effective degrees of freedom and BIC
# EDF = trace(H_lambda^{-1} * H_0) where:
# H_lambda = X^T W X + lambda*S (Hessian with penalty)
# H_0 = X^T W X (Hessian without penalty)
# W = diag(y_i / mu_i^2)
compute_bic <- function(gamma, lambda, X, y, S) {
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Negative log-likelihood (without penalty)
  nll_value <- -sum(y * log(mu) - mu)
  
  # Weight matrix W for Hessian computation
  # W_ii = y_i / mu_i^2
  W <- y / (mu^2)
  
  # Compute Hessian matrices
  # These are w.r.t. beta, not gamma
  XtWX <- t(X) %*% (W * X)  # X^T W X (using vectorized multiplication)
  H_lambda <- XtWX + lambda * S
  H_0 <- XtWX
  
  # Compute effective degrees of freedom
  # EDF = trace(H_lambda^{-1} * H_0)
  # Use solve to compute H_lambda^{-1} * H_0
  H_inv_H0 <- solve(H_lambda, H_0)
  edf <- sum(diag(H_inv_H0))
  
  # BIC = -2*log-likelihood + log(n)*EDF
  bic <- 2 * nll_value + log(length(y)) * edf
  
  return(list(bic = bic, edf = edf, nll = nll_value))
}

# Grid search over log(lambda)
log_lambda_seq <- seq(-13, -7, length = 50)
lambda_seq <- exp(log_lambda_seq)

# Storage for results
bic_values <- numeric(length(lambda_seq))
edf_values <- numeric(length(lambda_seq))
nll_values <- numeric(length(lambda_seq))
gamma_fits <- matrix(0, length(gamma_init), length(lambda_seq))

# Progress indicator
cat("Optimizing over", length(lambda_seq), "lambda values...\n")

# Loop over lambda values
for (i in 1:length(lambda_seq)) {
  lambda_current <- lambda_seq[i]
  
  # Optimize for current lambda
  fit_current <- optim(
    par = gamma_init,  # Use same starting values
    fn = nll,
    gr = grad_nll,
    X = X,
    y = y,
    S = S,
    lambda = lambda_current,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Store results
  gamma_fits[, i] <- fit_current$par
  
  # Compute BIC
  bic_result <- compute_bic(fit_current$par, lambda_current, X, y, S)
  bic_values[i] <- bic_result$bic
  edf_values[i] <- bic_result$edf
  nll_values[i] <- bic_result$nll
  
  # Progress update every 10 iterations
  if (i %% 10 == 0) {
    cat("  Completed", i, "of", length(lambda_seq), "\n")
  }
}

# Find optimal lambda
optimal_idx <- which.min(bic_values)
lambda_optimal <- lambda_seq[optimal_idx]
log_lambda_optimal <- log_lambda_seq[optimal_idx]
gamma_optimal <- gamma_fits[, optimal_idx]

cat("\nOptimal log(lambda):", log_lambda_optimal, "\n")
cat("Optimal lambda:", lambda_optimal, "\n")
cat("Minimum BIC:", bic_values[optimal_idx], "\n")
cat("Effective degrees of freedom at optimum:", edf_values[optimal_idx], "\n")

# Plot BIC curve
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

plot(log_lambda_seq, bic_values,
     type = "l",
     lwd = 2,
     xlab = "log(lambda)",
     ylab = "BIC",
     main = "BIC vs log(lambda)")
points(log_lambda_optimal, bic_values[optimal_idx], 
       col = "red", pch = 19, cex = 1.5)
abline(v = log_lambda_optimal, lty = 2, col = "red")

plot(edf_values, bic_values,
     type = "l",
     lwd = 2,
     xlab = "Effective Degrees of Freedom",
     ylab = "BIC",
     main = "BIC vs EDF")
points(edf_values[optimal_idx], bic_values[optimal_idx],
       col = "red", pch = 19, cex = 1.5)

################################################################################
# STEP 8: Bootstrap to assess uncertainty
################################################################################

cat("\n=== Bootstrap Analysis (200 replicates) ===\n")

n_bootstrap <- 200

# Storage for bootstrap results
# Each column is f(t) for one bootstrap replicate
f_bootstrap <- matrix(0, length(t_ext), n_bootstrap)

# Progress tracking
cat("Running bootstrap...\n")

# Bootstrap loop
for (b in 1:n_bootstrap) {
  # Generate bootstrap weights
  # Equivalent to resampling data with replacement
  wb <- tabulate(sample(n, replace = TRUE), n)
  
  # Modified negative log-likelihood with bootstrap weights
  nll_bootstrap <- function(gamma, X, y, S, lambda, weights) {
    beta <- exp(gamma)
    mu <- X %*% beta
    mu <- as.vector(mu)
    
    # Weighted log-likelihood
    ll <- sum(weights * (y * log(mu) - mu))
    penalty <- (lambda / 2) * sum(beta * (S %*% beta))
    
    return(-ll + penalty)
  }
  
  # Modified gradient with bootstrap weights
  grad_nll_bootstrap <- function(gamma, X, y, S, lambda, weights) {
    beta <- exp(gamma)
    mu <- X %*% beta
    mu <- as.vector(mu)
    
    grad_ll_beta <- t(X) %*% (weights * (y / mu - 1))
    grad_penalty_beta <- lambda * (S %*% beta)
    grad <- -beta * (grad_ll_beta - grad_penalty_beta)
    
    return(as.vector(grad))
  }
  
  # Optimize for bootstrap sample
  fit_boot <- optim(
    par = gamma_optimal,  # Start from optimal parameters
    fn = nll_bootstrap,
    gr = grad_nll_bootstrap,
    X = X,
    y = y,
    S = S,
    lambda = lambda_optimal,
    weights = wb,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Calculate f(t) for this bootstrap replicate
  beta_boot <- exp(fit_boot$par)
  f_bootstrap[, b] <- as.vector(matrices$X_tilde %*% beta_boot)
  
  # Progress update
  if (b %% 20 == 0) {
    cat("  Completed", b, "of", n_bootstrap, "bootstrap replicates\n")
  }
}

# Calculate 95% confidence intervals
# 2.5th and 97.5th percentiles for each time point
f_lower <- apply(f_bootstrap, 1, quantile, probs = 0.025)
f_upper <- apply(f_bootstrap, 1, quantile, probs = 0.975)

# Also calculate mean bootstrap estimate
f_mean <- apply(f_bootstrap, 1, mean)

cat("\nBootstrap analysis complete\n")

################################################################################
# STEP 9: Final comprehensive plot
################################################################################

cat("\n=== Creating final plot ===\n")

# Calculate final fitted values with optimal lambda
beta_optimal <- exp(gamma_optimal)
mu_optimal <- as.vector(X %*% beta_optimal)
f_optimal <- as.vector(matrices$X_tilde %*% beta_optimal)

# Create comprehensive plot with two panels
par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

# Panel 1: Deaths - Actual vs Fitted
plot(t, y,
     type = "p",
     pch = 19,
     cex = 0.8,
     col = "black",
     xlab = "Day of Year 2020",
     ylab = "Number of Deaths",
     main = "COVID-19 Deaths: Actual vs Model Fit",
     ylim = c(0, max(y) * 1.1))

lines(t, mu_optimal, col = "red", lwd = 2)

legend("topright",
       legend = c("Observed Deaths", "Model Fit"),
       col = c("black", "red"),
       pch = c(19, NA),
       lty = c(NA, 1),
       lwd = c(NA, 2),
       bty = "n")

# Add text with model information
text(min(t) + 5, max(y) * 1.05,
     paste0("Optimal log(λ) = ", round(log_lambda_optimal, 2)),
     pos = 4, cex = 0.8)

# Panel 2: Inferred Daily Infections with 95% CI
plot(t_ext, f_optimal,
     type = "l",
     col = "blue",
     lwd = 2,
     xlab = "Day of Year 2020",
     ylab = "Daily New Infections",
     main = "Inferred Daily Infection Rate with 95% Confidence Intervals",
     ylim = c(0, max(f_upper) * 1.1))

# Add confidence interval as shaded region
polygon(c(t_ext, rev(t_ext)),
        c(f_lower, rev(f_upper)),
        col = rgb(0, 0, 1, 0.2),  # Transparent blue
        border = NA)

# Add confidence interval boundary lines
lines(t_ext, f_lower, col = "blue", lty = 2, lwd = 1)
lines(t_ext, f_upper, col = "blue", lty = 2, lwd = 1)

# Re-draw main curve on top
lines(t_ext, f_optimal, col = "blue", lwd = 2)

# Mark first observation day
abline(v = min(t), lty = 2, col = "gray")

legend("topright",
       legend = c("Estimated f(t)", "95% CI", "First Death"),
       col = c("blue", "blue", "gray"),
       lty = c(1, 2, 2),
       lwd = c(2, 1, 1),
       bty = "n")

cat("\nFinal plot created successfully\n")

################################################################################
# STEP 10: Summary statistics and results
################################################################################

cat("\n=== FINAL RESULTS SUMMARY ===\n")
cat("==============================\n")
cat("Data: English NHS COVID-19 deaths, year 2020\n")
cat("Number of observations:", n, "\n")
cat("Period:", min(t), "to", max(t), "(days of year)\n")
cat("\nModel Specifications:\n")
cat("- Number of B-spline basis functions:", ncol(X), "\n")
cat("- Infection-to-death distribution: log-normal(3.151, 0.469^2)\n")
cat("\nOptimization Results:\n")
cat("- Optimal log(lambda):", round(log_lambda_optimal, 4), "\n")
cat("- Optimal lambda:", format(lambda_optimal, scientific = TRUE), "\n")
cat("- Minimum BIC:", round(bic_values[optimal_idx], 2), "\n")
cat("- Effective degrees of freedom:", round(edf_values[optimal_idx], 2), "\n")
cat("\nModel Fit:\n")
cat("- Deviance:", round(2 * nll_values[optimal_idx], 2), "\n")
cat("- Mean fitted deaths:", round(mean(mu_optimal), 2), "\n")
cat("- Mean observed deaths:", round(mean(y), 2), "\n")
cat("\nPeak Infection Estimates:\n")
cat("- Peak infection rate:", round(max(f_optimal), 1), "per day\n")
cat("- Day of peak infection:", t_ext[which.max(f_optimal)], "\n")
cat("- 95% CI at peak: [", round(f_lower[which.max(f_optimal)], 1), ",",
    round(f_upper[which.max(f_optimal)], 1), "]\n")
cat("\nBootstrap Analysis:\n")
cat("- Number of bootstrap replicates:", n_bootstrap, "\n")
cat("- Average CI width:", round(mean(f_upper - f_lower), 1), "\n")
cat("==============================\n")

################################################################################
# END OF ANALYSIS
################################################################################

cat("\n=== Analysis Complete ===\n")
cat("All plots have been generated in the graphics window.\n")
cat("Code executed successfully!\n")
