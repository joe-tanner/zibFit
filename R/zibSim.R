zibSim <- function(n, t, a, b, sigma1, sigma2, phi, alpha, beta, X, Z, seed = NULL) {

  # n = number of samples to simulate
  # t = number of time points to simulate
  # a = random effect mean (probability of being 0)
  # b = random effect mean (average of beta distribution)
  # sigma1 = random effect variance (probability of being 0)
  # sigma2 = random effect variance (average of beta distribution)
  # phi = beta dispersion parameter
  # alpha = regression coefficients (same length as the number number of columns in X)
  # beta = regression coefficients (same length as the number number of columns in Z)
  # X = covariate matrix (logistic) (must include 'Subject' and 'Time' columns)
  # Z = covariate matrix (beta) (must include 'Subject' and 'Time' columns)
  # seed = simulation seed

  # check inputs
  if(length(n) != 1 | length(t) != 1) stop("n, t must be integers of length 1")
  if(!is.numeric(n)) stop("n must be an integer of length 1")
  if(!is.numeric(t)) stop("t must be an integer of length 1")

  if(length(a) != 1 | length(b) != 1) stop("a, b must be integers of length 1")
  if(!is.numeric(a)) stop("a must be an integer of length 1")
  if(!is.numeric(b)) stop("b must be an integer of length 1")

  if(length(sigma1) != 1 | length(sigma2) != 1) stop("sigma 1, 2 must be integers of length 1")
  if(!is.numeric(sigma1)) stop("sigma 1 must be an integer of length 1")
  if(!is.numeric(sigma2)) stop("sigma 2 must be an integer of length 1")
  if(sigma1 < 0) stop("sigma 1 must be strictly greater than 0")
  if(sigma2 < 0) stop("sigma 2 must be strictly greater than 0")

  if(length(phi) != 1) stop("phi must be an integer of length 1")
  if(!is.numeric(phi)) stop("phi must be an integer of length 1")

  if(length(alpha) != ncol(X)) stop("alpha must have the same number of integers as X has columns")
  if(length(beta) != ncol(Z)) stop("beta must have the same number of integers as Z has columns")

  if(n*t != nrow(X)) stop("covariate X must be the same length as the total number of observations (n*t)")
  if(n*t != nrow(Z)) stop("covariate Z must be the same length as the total number of observations (n*t)")

  if(is.null(seed)) seed = round(a*b*n*t*phi, 0)

  set.seed(seed)
  results = rep(0, times = n)

  # sample random effects from normal distributions
  sample_a <- rnorm(n = n, mean = a, sd = sqrt(sigma1))
  sample_b <- rnorm(n = n, mean = b, sd = sqrt(sigma2))

  for(i in 1:n) {
    for(j in 1:t) {
      m = (i-1)*t + j # works provided it's ordered correctly - look at??

      inv_logit <- function(x) exp(x)/(1 + exp(x))
      p = inv_logit(sample_a[i] + (X[m, ]%*%alpha))

      if(runif(1) < p) {
        u = inv_logit(sample_b[i] + (Z[m,]%*%beta))
        y = rbeta(n = 1, shape1 = u*phi, shape2 = (1-u)*phi)
      } else {
        y = 0
      }

      # catch all just in case numbers are going to be rounded to 1
      if (y > 1 - 10^(-6)) {
        y <- 1 - 10^(-6)
      }

      results[m] = y

    }
  }

  subject_ind = rep(1:n, each = t)
  time_ind = rep(1:t, n)
  sample_ind = rep(1:(n*t))
  parameters = list(a = a, b = b, sigma1 = sigma1, sigma2 = sigma2, phi = phi,
                    alpha = alpha, beta = beta)
  log_covariates = data.frame(X)
  beta_covariates = data.frame(Z)

  l <- list(rel_abundance = results,
            log_covariates = log_covariates,
            beta_covariates = beta_covariates,
            subject_ind = subject_ind,
            time_ind = time_ind,
            sample_ind = sample_ind,
            theta = parameters)

  class(l) <- c("zibSim", "listof")

  return(l)
}
