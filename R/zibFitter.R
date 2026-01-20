#' Fit Zero-Inflated Beta Regression Models to Longitudinal Microbiome
#' Composition Data
#'
#' Given an object of \code{"zibData"}, fit a zero-inflated beta regression
#' model using the \code{ZIBR} package. The zero-inflated beta
#' regression model, as proposed by Eric Z. Chen and Hongzhe Li in 2016,
#' proposes a two part zero-inflated beta regression model with random effects
#' for testing the association between microbial abundance and clinical
#' covariates for longitudinal microbiome data. The model includes a logistic
#' component to model presence/absence of the microbe in samples, and a Beta
#' component to model non-zero microbial abundance. Both components include
#' random effects to take into account the correlation among repeated
#' measurements on the same subject. Basic details on the model are as follows;
#' \deqn{
#' y_{it} \sim
#' \begin{cases}
#' 0 & \text{with probability } 1 - p_{it} \\
#' \mathrm{Beta}(\mu_{it}\phi,\ (1-\mu_{it})\phi) & \text{with probability } p_{it}
#' \end{cases}
#' }
#' This function models the \code{"zibData"} object, returning
#' a list of each microbe with statistically significant (p-value < 0.05)
#' estimates for each covariate included in the model. See
#' \code{\link{zibClean}} for more information on \code{"zibData"} objects.
#' This function requires a minimum of two different microbes to use it.
#'
#' @param zibData An object of class \code{"zibData"}, as returned from the
#' \code{\link{zibClean}} function.
#'
#' @returns A list including the following;
#'  \item{\code{p.species}: }{List of the p-values associated with each parameter
#'  estimate for each microbe included in the model.}
#'  \item{\code{sig_stats}: }{List detailing each microbe with statistically
#'  significant (p-value < 0.05) estimates for each covariate.}
#' If p-values are returned as NA, this means that the data wasn't suitable for
#' modelling with the \code{zibr} function, usually meaning that there are
#' either too many/few zeroes or the relative abundance of non-zero samples are
#' too small to fit a model to; please see \code{ZIBR::zibr} for more
#' information.
#'
#' @export
#'
#' @importFrom ZIBR "zibr"
#'
#' @seealso \code{zibClean}, \code{zibSim}, \code{ZIBR::zibr}
#'
#' @examples
#' library(zibFit)
#' # simulate zibr data and plot
#'
#' # initialise covariate matrix, simulating a clinical study where half
#' # recieve the treatment (1) and half do not (0)
#' n = 1000
#' t = 3
#' X <- as.matrix(c(rep(0, n/2 * t), rep(1, n/2 * t)))
#'
#' y1 <- zibSim(n = n, t = t, a = 0.5, b = 0.8, sigma1 = 2.5, sigma2 = 0.5,
#' phi = 7.2, alpha = 0.5, beta = -0.5, X = X, Z = X, seed = 6874)
#'
#' # fit a zibr model to longitudinal microbiome composition data
#' # simulate five different microbiomes
#' y2 <- zibSim(n = n, t = t, a = 0.5, b = -0.5, sigma1 = 0.8, sigma2 = 1.2,
#' phi = 5, alpha = -0.5, beta = -0.5, X = X, Z = X, seed = 6874)
#' y3 <- zibSim(n = n, t = t, a = 0.2, b = 0.3, sigma1 = 3, sigma2 = 2,
#' phi = 4, alpha = -0.5, beta = -0.5, X = X, Z = X, seed = 6874)
#' y4 <- zibSim(n = n, t = t, a = 0.1, b = 0.9, sigma1 = 0.4, sigma2 = 0.5,
#' phi = 8.1, alpha = 0.5, beta = 0.5, X = X, Z = X, seed = 6874)
#' y5 <- zibSim(n = n, t = t, a = 0.2, b = 1, sigma1 = 2, sigma2 = 0.5,
#' phi = 3.4, alpha = -0.5, beta = 0.5, X = X, Z = X, seed = 6874)
#'
#' ra <- cbind(y1$rel_abundance, y2$rel_abundance, y3$rel_abundance,
#' y4$rel_abundance, y5$rel_abundance)
#' rownames(ra) <- rep(1:3000) # rownames and column names are required
#' colnames(ra) <- rep(1:5)
#'
#' # create covariate data frame
#' sample_id <- rep(1:3000)
#' covariates <- data.frame(subject = y1$subject_ind, time = y1$time_ind,
#' treat = y1$log_covariates$X, sample = sample_id)
#'
#' zibData = zibClean(data = ra,
#' metadata = covariates,
#' log_covs = c("subject", "time", "treat"),
#' id_column = "sample",
#' subject_column = "subject",
#' time_column = "time")
#'
#' results <- zibFitter(zibData)
#'
zibFitter <- function(zibData) {

  if(!inherits(zibData, "zibData")) stop("This function only works on objects of class \"zibData\"")

  data <- zibData
  spe.all <- data$spec_names
  p.species.list.zibr <- list()

  for(spe in spe.all) {

    # extract covariates
    X <- data$log_covariates
    Z <- data$beta_covariates

    # include baseline abundance
    baseline_abundance <- data$rel_abundance[X$baseline_sample, spe]
    X <- cbind(X, baseline_abundance)[ , -which(names(X) == "baseline_sample")]
    Z <- cbind(Z, baseline_abundance)[ , -which(names(Z) == "baseline_sample")]

    subject.ind <- data$subject_ind
    time.ind <- data$time_ind
    sample.ind <- data$sample_ind

    Y <- data$rel_abundance[sample.ind, spe] # check that log and beta must have same indices

    print(spe)

    if(sum(Y > 0) < 10 || sum(Y == 0) < 10 || max(Y) < 0.01) {
      p.species.list.zibr[[spe]] <- NA
    } else {
      est <- ZIBR::zibr(
        logistic_cov = X, beta_cov = Z, Y = Y,
        subject_ind = subject.ind,
        time_ind = time.ind,
        quad_n = 30, verbose = TRUE
      )
      p.species.list.zibr[[spe]] <- est$joint_p
    }
  }

  p.species.zibr <- t(as.data.frame(p.species.list.zibr))
  vars <- colnames(p.species.zibr)

  # will need to change this code to base R
  p.species.zibr.adj <- data.frame(
    Species = rownames(p.species.zibr),
    p.species.zibr
  )

  sig <- list()
  for(name in vars) {
    sig[[name]] <- stats::na.omit(p.species.zibr.adj[p.species.zibr.adj[[name]] < 0.05, "Species"])
  }

  return <- list(p.species = p.species.list.zibr,
                 sig_stats = sig)

  return(return)

}
