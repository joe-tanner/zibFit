#' Prepare longitudinal microbiome composition data for modelling
#'
#' Given two data sets, one detailing the relative abundances of different
#' microbes for different samples (\code{data}) and one providing information
#' on the covariates for each individual sample (\code{metadata}), this function
#' tidies the data up into a \code{"zibData"} object, suitable for fitting
#' zero-inflated beta regression models to the data using the
#' \code{\link{zibFitter}} package.
#'
#' @param data Response variable (relative abundance) for the regression model;
#' row names must equal the \code{id_column} vector and column names must
#' equal the names of the different microbes to analyse.
#' @param metadata Metadata including information for each individual sample on
#' all covariates to be included in the model.
#' @param log_covs Names of the covariates to be included in the log part of the
#' model - should equal names of columns in the \code{metadata} data set.
#' @param beta_covs Names of the covariates to be included in the beta part of
#' the model - should equal names of columns in the \code{metadata} data set.
#' Equal to the \code{log_covs} by default.
#' @param id_column Name of the column in \code{metadata} corresponding to the
#' ID of each sample; should correspond to the rownames of \code{data}; should
#' be inputted as a character string.
#' @param subject_column Name of the column in \code{metadata} corresponding to
#' the subject each sample is associated with; should be inputted as a
#' character string.
#' @param time_column Name of the column in \code{metadata} corresponding to the
#' time point each sample was taken at; should be inputted as a character
#' string.
#'
#' @returns An object of class \code{"zibData"}, which is a list containing the
#' following items;
#'  \item{\code{rel_abundance}: }{Relative abundance of each species for each sample,
#'  of class \code{data.frame}.}.
#'  \item{\code{log_covariates}: }{Covariate matrix X used to simulate the data
#'  (logistic part), of class \code{data.frame}.}.
#'  \item{\code{beta_covariates}: }{Covariate matrix Z used to simulate the data
#'  (beta part), of class \code{data.frame}.}
#'  \item{\code{spec_names}: }{Vector of the unique species names to be modelled}.
#'  \item{\code{subject_ind}: }{Ordered index of subject ID's.}
#'  \item{\code{time_ind}: }{Ordered index of time points for each data point.}
#'  \item{\code{sample_ind}: }{Ordered index of sample ID's.}
#'
#' @export
#' @seealso \code{zibFitter}, \code{zibSim}
#'
#' @examples
#' library(zibFit)
#' # simulate zibr data and plot
#'
#' # initialise covariate matrix, simulating a clinical study where half recieve
#' # the treatment (1) and half do not (0)
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
zibClean <- function(data, metadata, log_covs, beta_covs = log_covs, id_column,
                     subject_column, time_column) {

  # id, subject and time are the names of the metadata columns in ""
  # log,beta are string vectors with the names of the metadata columns holding the covs
  # rownames of the data should equal the id values in the 'id_column' in metadata

  # data doesn't need to be the same size as the metadata, function will sort that
  # size of beta and log covs shouldn't matter either
  # data must be between 0 and 1

  if(!is.data.frame(metadata)) stop("metadata must be a dataframe")
  if(length(colnames(data)) == 0) stop("data requires column names of the different species in each column")
  if(length(rownames(data)) == 0) stop("data requires row names of the sample id's for each row")

  covs <- unique(c(log_covs, beta_covs))
  options <- colnames(metadata)
  if(!time_column %in% options) stop("time_column not found in metadata")
  if(!subject_column %in% options) stop("subject_column not found in metadata")
  if(!id_column %in% options) stop("id_column not found in metadata")

  missing <- covs[!covs %in% c(options)]
  if(length(missing > 0)) stop("some selected covariates cannot be found in the metadata dataset")

  spe_names <- colnames(data)
  sample_id <- rownames(data)

  sample_id <- sample_id[which(sample_id %in% metadata[[id_column]])]
  if(length(sample_id) == 0) stop("zero sample ids from data match the sample ids in metadata")
  data <- data[sample_id, ]

  # convert time if not numeric
  if(!is.numeric(metadata[[time_column]])) {

    time_levels <- unique(metadata[[time_column]])
    time_levels <- time_levels[order(as.numeric(sub("\\D+", "", time_levels)))]
    metadata[[time_column]] <- as.numeric(as.character(factor(metadata[[time_column]],
                                                              levels = time_levels,
                                                              labels = rep(0:(length(time_levels) - 1)))))
  }

  colnames(metadata)[which(colnames(metadata) == time_column)] <- "Time"

  # figure out which covs are wanted
  covs_new <- covs[!covs %in% c(time_column, id_column)]

  # convert other variables
  for(i in 1:length(covs_new)) {

    var = covs_new[i]
    temp = unique(metadata[[var]])
    metadata[[var]] <- as.numeric(as.character(factor(metadata[[var]],
                                                      levels = temp,
                                                      labels = rep(1:(length(temp))))))

  }

  # create covariates
  reg.cov <- data.frame(sample_id = sample_id, stringsAsFactors = FALSE)
  reg.cov <- merge(reg.cov, metadata, by.x = "sample_id", by.y = id_column,
                   all.y = FALSE, all.x = FALSE)

  # remove duplicates
  reg.cov <- reg.cov[!duplicated(reg.cov), ]
  rownames(reg.cov) <- reg.cov$sample_id

  # take out the baseline sample
  t0 <- min(metadata$Time)
  reg.cov.t1 <- subset(reg.cov, Time == t0)
  reg.cov.t1 <- reg.cov.t1[c("sample_id", subject_column)]
  reg.cov.t <- subset(reg.cov, Time != t0)
  sample_id <- reg.cov.t$sample_id

  all_cov <- merge(reg.cov.t, reg.cov.t1, by = subject_column)
  all_cov <- all_cov[match(sample_id, all_cov$sample_id.x), ]
  colnames(all_cov)[which(colnames(all_cov) == "sample_id.y")] <- "baseline_sample"
  colnames(all_cov)[which(colnames(all_cov) == "sample_id.x")] <- id_column
  rownames(all_cov) <- sample_id
  all_cov <- all_cov[ , -which(names(all_cov) == id_column)]

  colnames(all_cov)[which(colnames(all_cov) == "Time")] <- time_column

  # separate into log and beta covariates
  log_cov <- all_cov[c(log_covs, "baseline_sample")]
  beta_cov <- all_cov[c(beta_covs, "baseline_sample")]

  # update indices without baseline abundances
  time_ind <- as.numeric(all_cov[[time_column]])
  sample_ind <- all_cov[[id_column]]
  subject_ind <- as.numeric(all_cov[[subject_column]])

  # return
  list <- list(rel_abundance = data,
               log_covariates = log_cov,
               beta_covariates = beta_cov,
               spec_names = spe_names,
               subject_ind = subject_ind,
               time_ind = time_ind,
               sample_ind = sample_id)

  class(list) <- c("zibData", "listof")
  return(list)

}

