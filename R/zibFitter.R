zibrFitter <- function(zibData) {

  if(!inherits(zibData, "zibData")) stop("This function only works on objects of class \"zibData\"")

  # library(ZIBR)

  data <- zibData
  spe.all <- data$spec_names
  p.species.list.zibr <- list()

  for(spe in spe.all) {

    # extract covariates
    X <- data$log_covariates
    Z <- data$beta_covariates

    # include baseline abundance
    baseline_abundance <- data$rel_abundance[X$baseline_sample, spe]/100
    X <- cbind(X, baseline_abundance)[ , -which(names(X) == "baseline_sample")]
    Z <- cbind(Z, baseline_abundance)[ , -which(names(Z) == "baseline_sample")]

    subject.ind <- data$subject_ind
    time.ind <- data$time_ind
    sample.ind <- data$sample_ind


    Y <- data$rel_abundance[sample.ind, spe]/100 # check that log and beta must have same indices

    print(spe)

    if(sum(Y > 0) < 10 || sum(Y == 0) < 10 || max(Y) < 0.01) {
      p.species.list.zibr[[spe]] <- NA
    } else {
      est <- zibr(
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
  p.species.zibr.adj <-
    add_rownames(as.data.frame(p.species.zibr), var = "Species") %>% mutate_each(funs(p.adjust(., "fdr")), -Species)

  sig <- list()
  for(name in vars) {
    sig[[name]] <- na.omit(p.species.zibr.adj[p.species.zibr.adj[[name]] < 0.05, "Species"])
  }

  return <- list(p.species = p.species.list.zibr,
                 sig_stats = sig)

  return(return)

}
