# loading data
load("dataex2.rdata")
load("dataex4.rdata")
load("dataex5.rdata")



### --- Q2 -------------------------
# Determine the maximum likelihood estimate of mu based on the data available 
# in the file dataex2.Rdata. Consider sigma2 known and equal to 1.5^2
mean.all <- mean(dataex2$X)
sigma2 <- 1.5

# Data extraction
X <- dataex2$X
R <- dataex2$R

# Define the negative log-likelihood function for the model
n.log.ll <- function(mu, sigma2, x, r) {
  log.ll <- 0
  for (i in 1:length(x)){
    # Update the log-likelihood value by summing over all data points
    log.ll <- log.ll - (r[i] * dnorm(x[i], mu, sqrt(sigma2)) + 
                         (1 - r[i]) * pnorm(x[i], mu, sqrt(sigma2)))
  }
  log.ll
}

# Optimize the negative log-likelihood function to find the maximum likelihood 
# estimate of mu
mle <- suppressWarnings(optim(mean.all, n.log.ll, 
                                 sigma2 = sigma2, x = X, r = R))

print(paste('The maximum likelihood estimate of mu based on the data is',
            mle$par))




### --- Q4 --------------------------
x <- dataex4$X
y <- dataex4$Y

# Extract observed and missing data
y.obs <- y[!is.na(y)]
x.obs <- x[!is.na(y)]
x.miss <- x[is.na(y)]

# Define the negative log-likelihood function
nll.q4 <- function(beta, y, x) {
  p <- exp(beta[1] + x * beta[2]) / (1 + exp(beta[1] + x * beta[2]))
  -sum(y * log(p) + (1 - y) * log(1 - p))
}

# EM algorithm
# Initialize the parameter estimates and convergence settings
beta <- c(0, 0)
tolerance <- 1e-4
converged <- FALSE

while (converged==FALSE) {
  # E-step: estimate the missing response values based on the current beta 
  # values
  miss <- exp(beta[1] + x.miss * beta[2]) / 
    (1 + exp(beta[1] + x.miss * beta[2]))
  
  # M-step: update the beta values using the complete data 
  # (observed + estimated missing values)
  y.comp <- c(y.obs, miss)
  x.comp <- c(x.obs, x.miss)
  
  # Optimize the negative log-likelihood function
  mle.beta <- optim(beta, function(beta) nll.q4(beta, y.comp, x.comp))
  beta.new <- mle.beta$par
  
  # Check convergence: if the change in beta values is less than 
  # the tolerance, break
  if (sum(abs(beta.new - beta)) < tolerance) {
    converged <- TRUE
    beta <- beta.new
  } else {
    beta <- beta.new
  }
}

print(paste('The maximum likelihood estimate of beta0 based on the data is',
            beta[1]))
print(paste('The maximum likelihood estimate of beta1 based on the data is',
            beta[2]))



### --- Q5 -------------------------------
y <- dataex5

# Define the density functions for the mixture model components
fx <- function(y, lambda) {
  lambda * y^(-lambda - 1)
}

fy <- function(y, mu) {
  mu * y^(-mu - 1)
}

# EM algorithm
# Initialize the parameter estimates and convergence settings
theta <- c(0.3, 0.3, 0.4)
tolerance <- 1e-4
converged <- FALSE

while (converged == FALSE) {
  # E-step: Compute the posterior probability of each observation 
  # belonging to component 1
  p <- (theta[1] * fx(y, theta[2])) / 
    (theta[1] * fx(y, theta[2]) + (1 - theta[1]) * fy(y, theta[3]))
  
  # M-step: Update the parameter estimates
  theta.new <- theta
  # Update the proportion of component 1
  theta.new[1] <- mean(p)
  mle.lambda <- suppressWarnings(
    optim(theta[2], 
          function(lambda) -sum(p * (log(lambda) - (lambda + 1) * log(y)))))
  # Update lambda using the MLE
  theta.new[2] <- mle.lambda$par
  mle.mu <- suppressWarnings(
    optim(theta[3], 
          function(mu) -sum((1 - p) * (log(mu) - (mu + 1) * log(y)))))
  # Update mu using the MLE
  theta.new[3] <- mle.mu$par
  
  # Check convergence: If the change in parameter estimates is less 
  # than the tolerance, break
  if (sum(abs(theta.new - theta)) < tolerance) {
    converged <- TRUE
    theta <- theta.new
  } else {
    theta <- theta.new
  }
}

print(paste('The maximum likelihood estimate of p based on the data is',
            theta[1]))
print(paste('The maximum likelihood estimate of lambda based on the data is',
            theta[2]))
print(paste('The maximum likelihood estimate of mu based on the data is',
            theta[3]))

# Calculate density
density <- theta[1] * fx(sort(y), theta[2]) + 
  (1 - theta[1]) * fy(sort(y), theta[3])

# Plot histogram and density
hist(y, freq = FALSE, breaks = 'FD', border = "black", 
     xlab = "Y", ylab = "Density", xlim = c(1,8),
     main = "Histogram of Y with Estimated Density Superimposed")
lines(sort(y), density, col = "red", lwd = 2)