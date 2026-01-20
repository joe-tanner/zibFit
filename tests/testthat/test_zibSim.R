library(zibFit)

test_that("zibSim produces valid data", {
  expect_s3_class(zibSim(n = 1000, t = 3, a = 0.1, b = 0.9, sigma1 = 3.2, sigma2 = 2.1,
                         phi = 3.2, alpha = -0.5, beta = -0.5,
                         X = as.matrix(c(rep(0, 1500), rep(1, 1500))),
                         Z = as.matrix(c(rep(0, 1500), rep(1, 1500))),
                         seed = 6874), "zibSim")
  expect_error(zibSim(n = 1000, t = 3, a = 0.1, b = 0.9, sigma1 = -3.2, sigma2 = 2.1,
                      phi = 3.2, alpha = -0.5, beta = -0.5,
                      X = as.matrix(c(rep(0, 1500), rep(1, 1500))),
                      Z = as.matrix(c(rep(0, 1500), rep(1, 1500))),
                      seed = 6874))
  expect_error(zibSim(n = 1000, t = 3, a = 0.1, b = 0.9, sigma1 = 3.2, sigma2 = 2.1,
                      phi = 3.2, alpha = -0.5, beta = -0.5,
                      X = as.matrix(c(rep(0, 750), rep(1, 1500))),
                      Z = as.matrix(c(rep(0, 1500), rep(1, 1500))),
                      seed = 6874))
})

