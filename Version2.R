## proj3
library(splines)
dat <- read.table("engcov (1).txt", header = TRUE)
## Extract relevant variables
y <- dat$nhs # Deaths in NHS hospitals
t <- dat$julian # Day of the year 2020
n <- length(y) # Number of observations

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
  
  ## Create time points where f(t) will be evaluated
  ## From t[1]-30 to t[n], covering all possible infection times (integers).
  t_cover <- seq(t_min, t_max, by=1)
  cat(length(t_cover))
  
  ## Generate B-spline basis matrix X_tilde
  ## Each column is a B-spline basis function evaluated at t_cover
  ## Unless outer.ok is true, the values in x must be between the “inner” knots knots[ord] and knots[ length(knots) - (ord-1)].
  X_tilde <- splineDesign(knots, t_cover, outer.ok = TRUE)
  
  ## Initialize model matrix X for deaths
  ## X[i,j] gives contribution of beta_j to expected deaths on day t[i]
  X <- matrix(0, n, k)
  
  ## Build X by convolving X_tilde with probability distribution pd
  ## For each observation day i
  for (i in 1:n) {
    ## Determine range of infection days that contribute to deaths on day t[i]
    ## Upper limit: either 29 days before day i, or 80 days before, whichever is less
    j_max <- min(29 + i, 80)
    
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
  
  ## Create penalty matrix S
  ## Penalizes second differences: (beta_{k-1} - 2*beta_k + beta_{k+1})^2
  ## This encourages smooth changes in the beta coefficients
  D <- diff(diag(k), differences = 2) ## Second difference matrix
  S <- crossprod(D) ## S = D^T D
  
  ## Return all matrices as a list
  return(list(X = X, X_tilde = X_tilde, S = S, t_cover = t_cover))
}

## Create the matrices
matrices <- setup_matrices(t, pd, k = 80, backdays = 30)
X <- matrices$X
X_tilde <- matrices$X_tilde
S <- matrices$S

t_cover <- matrices$t_cover

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
  delta <- 1e-5  # Small step for finite difference
  grad_numerical <- numeric(length(gamma_test))
  
  for (i in 1:length(gamma_test)) {
    gamma_plus <- gamma_test
    gamma_plus[i] <- gamma_plus[i] + delta
    nll_plus <- nll(gamma_plus, X, y, S, lambda_test)
    
    gamma_minus <- gamma_test
    gamma_minus[i] <- gamma_minus[i] - delta
    nll_minus <- nll(gamma_minus, X, y, S, lambda_test)
    
    grad_numerical[i] <- (nll_plus - nll_minus) / (2 * delta)
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



