zibClean <- function(data, metadata, log_covs, beta_covs = log_covs, id_column,
                     subject_column, time_column) {

  # id, subject and time are the names of the metadata columns in ""
  # log,beta are string vectors with the names of the metadata columns holding the covs
  # rownames of the data should equal the id values in the 'id_column' in metadata

  # data doesn't need to be the same size as the metadata, function will sort that
  # size of beta and log covs shoudln't matter either

  if(!is.data.frame(metadata)) stop("metadata must be a dataframe")

  if(!time_column %in% options) stop("time_column not found in metadata")
  if(!subject_column %in% options) stop("subject_column not found in metadata")
  if(!id_column %in% options) stop("id_column not found in metadata")

  covs <- unique(c(log_covs, beta_covs))
  options <- colnames(metadata)
  missing <- covs[!covs %in% c(options)]
  if(length(missing > 0)) stop("some selected covariates cannot be found in the metadata dataset")

  spe_names <- colnames(data)
  sample_id <- rownames(data)

  sample_id <- sample_id[which(sample_id %in% metadata[[id_column]])]
  if(length(sample_id) == 0) stop("no sample ids from data match sample ids in metadata")
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
  reg.cov.t1 <- subset(reg.cov, Time == 0)
  reg.cov.t1 <- reg.cov.t1[c("sample_id", subject_column)]
  reg.cov.t <- subset(reg.cov, Time != 0)
  sample_id <- reg.cov.t$sample_id

  all_cov <- merge(reg.cov.t, reg.cov.t1, by = subject_column)
  all_cov <- all_cov[match(sample_id, all_cov$sample_id.x), ]
  colnames(all_cov)[which(colnames(all_cov) == "sample_id.y")] <- "baseline_sample"
  colnames(all_cov)[which(colnames(all_cov) == "sample_id.x")] <- id_column
  rownames(all_cov) <- sample_id
  all_cov <- all_cov[ , -which(names(all_cov) == id_column)]

  # separate into log and beta covariates
  log_cov <- all_cov[c(log_covs, "baseline_sample")]
  beta_cov <- all_cov[c(beta_covs, "baseline_sample")]

  # update indices without baseline abundances
  time_ind <- as.numeric(all_cov$Time)
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

