test_that("joint_summarize input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are
  #'   behaving as expected.

  # run traditional model to do tests with
  out <- traditional_model(data = green_crab_data, family = "negbin",
                           multicore = FALSE,
                           n_chain = 1, n_iter = 1000)

  #1. make sure model fit is of class stanfit
  data <- data.frame(y = c(1, 2, 3), x = c(1, 2, 3))
  lm_out <- lm(y ~ x, data = data)
  expect_error(joint_summarize(lm_out$model),
               "model_fit must be of class 'stanfit'.")

  #2. make sure probs is a numeric vector
  expect_error(joint_summarize(out$model, probs = "95%"),
               "probs must be a numeric vector.")

  #3. make sure all values of probs are between 0 and 1
  expect_error(joint_summarize(out$model, probs = c(5, 95)),
               "probs must be between 0 and 1.")

  #4. make sure par is a character vector
  expect_error(joint_summarize(out$model, par = c(1, 2, 3)),
               "par must be a character vector.")

  #5. make sure model fit contains all par input
  expect_error(joint_summarize(out$model, par = "alpha"),
               "model_fit must contain all selected parameters: alpha")
})

test_that("joint_summarize outputs work - q, phi, beta, p10", {
  testthat::skip_on_cran()

  ## 1.
  # model includes "p10", "q", "phi", "beta"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  phi <- 1.2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i, j] <- rnbinom(n = 1, mu = mu[i] * q, size = phi)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
  )

  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    phi = phi,
    q = q
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "phi", "q")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "negbin",
                initial_values = inits, n_warmup = 25, n_iter = 100,
                n_chain = 1, multicore = FALSE, seed = 10)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "phi", "alpha[1]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 4)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2", "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3]),
                  is.numeric(out[, 4])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1, mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - p10, beta, q", {
  testthat::skip_on_cran()

  ## 2.
  # model includes "p10", "beta", "q"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rpois(1, mu[i])
      } else {
        count[i, j] <- rpois(1, mu[i] * q)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    q = q
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "q")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, n_chain = 1, multicore = FALSE,
                seed = 10, initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 4)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2", "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3]),
                  is.numeric(out[, 4])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1, mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - gamma, p10, beta, q", {
  testthat::skip_on_cran()

  ## 3.
  # model includes "p10", "beta", "q", "alpha_gamma", "beta_gamma"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i], rate = beta_gamma)
      } else {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i] * q, rate = beta_gamma)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = beta,
    q = q
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha", "q")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "gamma",
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]") %in% output_params))

})

test_that("joint_summarize outputs work - p10, phi, beta", {
  testthat::skip_on_cran()

  ## 4.
  # model includes "p10", "phi", "beta"

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  phi <- 1.2
  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "phi")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin", initial_values = inits,
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "phi", "alpha[1]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional",
                                                 "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - p10, beta", {
  testthat::skip_on_cran()

  ## 5.
  # model includes "p10", "beta"

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, initial_values = inits,
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional",
                                                 "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1, mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - gamma, p10, beta", {
  testthat::skip_on_cran()

  ## 6.
  # model includes "p10", "beta", "alpha_gamma", "alpha_beta"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta <- 0.5
  log_p10 <- -4.5
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rgamma(nobs_count, shape = alpha_gamma[i], rate = beta_gamma)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = beta
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, initial_values = inits, family = "gamma",
                n_chain = 1, multicore = FALSE, seed = 10,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]") %in% output_params))

})

test_that("joint_summarize outputs work - p10, q, phi, alpha", {
  testthat::skip_on_cran()

  ## 7.
  # model includes "p10", "q", "phi", "alpha"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  phi <- 10
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i, j] <- rnbinom(n = 1, mu = mu[i] * q, size = phi)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    q = q,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha", "q", "phi")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin", q = TRUE,
                cov = c("var_a", "var_b"), n_chain = 1, multicore = FALSE,
                seed = 10, initial_values = inits, adapt_delta = 0.99,
                n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "phi", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                             cov_val = c(0, 0))

  # test dimensions
  expect_true(all(dim(out) == c(10, 4)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2", "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3]),
                  is.numeric(out[, 4])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1, cov_val = c(0, 0))

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - p10, alpha, q", {
  testthat::skip_on_cran()


  ## 8.
  # model includes "p10", "alpha", "q"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rpois(n = 1, mu[i])
      } else {
        count[i, j] <- rpois(n = 1, mu[i] * q)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]", "alpha[2]",
                    "alpha[3]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                             cov_val = c(0, 0))

  # test dimensions
  expect_true(all(dim(out) == c(10, 4)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2", "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3]),
                  is.numeric(out[, 4])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1, cov_val = c(0, 0))

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - gamma, p10, alpha, q", {
  testthat::skip_on_cran()


  ## 9.
  # model includes "p10", "alpha", "q", "alpha_gamma", "alpha_beta"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i], rate = beta_gamma)
      } else {
        count[i, j] <- rgamma(1, shape = alpha_gamma[i] * q, rate = beta_gamma)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    count_type = count_type,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = alpha,
    q = q
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha", "q")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, q = TRUE, family = "gamma",
                cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "q[1]", "alpha[1]", "alpha[2]",
                    "alpha[3]") %in% output_params))

})

test_that("joint_summarize outputs work - p10, phi, alpha", {
  testthat::skip_on_cran()


  ## 10.
  # model includes "p10", "phi", "alpha"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  phi <- 10

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha",
                         "phi")
  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "negbin",
                cov = c("var_a", "var_b"), n_warmup = 25,
                n_iter = 100, n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "phi", "alpha[1]",
                    "alpha[2]", "alpha[3]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                             cov_val = c(0, 0))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional",
                                                 "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1, cov_val = c(0, 0))

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - p10, alpha", {
  testthat::skip_on_cran()


  ## 11.
  # model includes "p10", "alpha"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rpois(nobs_count, mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("mu", "p10", "alpha")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]", "alpha[2]",
                    "alpha[3]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                             cov_val = c(0, 0))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional",
                                                 "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1, cov_val = c(0, 0))

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - gamma, p10, alpha", {
  testthat::skip_on_cran()


  ## 12.
  # model includes "p10", "alpha", "alpha_gamma", "beta_gamma"

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rgamma(nobs_count, shape = alpha_gamma[i], rate = beta_gamma)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[, 1] <- 1 # intercept
  for (i in 2:length(alpha)) {
    mat_site[, i] <- rnorm(nsite, 0, 1)
  }
  colnames(mat_site) <- c("int", "var_a", "var_b")
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA, nsite)
  p <- rep(NA, nsite)
  for (i in 1:nsite) {
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i, ] * alpha)))
    p[i] <- min(p11[i] + exp(log_p10), 1)
  }
  # pcr_n (# qPCR observations)
  pcr_n <- matrix(nobs_pcr, nrow = nsite, ncol = nobs_pcr)
  # pcr_k (# positive qPCR detections)
  pcr_k <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite) {
    pcr_k[i, ] <- rbinom(nobs_pcr, pcr_n[i, ], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    pcr_n = pcr_n,
    pcr_k = pcr_k,
    count = count,
    site_cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1, length(mu)),
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c("alpha_gamma", "beta_gamma", "p10", "alpha")

  # run model
  fit <- suppressWarnings({
    joint_model(data = data, family = "gamma", cov = c("var_a", "var_b"),
                n_chain = 1, multicore = FALSE, seed = 10,
                initial_values = inits, n_warmup = 25,
                n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("p10", "alpha[1]", "alpha[2]",
                    "alpha[3]") %in% output_params))

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                             cov_val = c(0, 0))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional",
                                                 "n_eDNA")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1, cov_val = c(0, 0))

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - q, phi", {
  testthat::skip_on_cran()


  ## 13.
  # model includes "q", "phi" (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  phi <- 1.2
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i, j] <- rnbinom(n = 1, mu = mu[i] * q, size = phi)
      }
    }
  }

  # collect data
  data <- list(
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "phi")
  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, q = TRUE, family = "negbin",
                      n_chain = 1, multicore = FALSE, seed = 10,
                      initial_values = inits,
                      n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("q[1]", "phi") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - q", {
  testthat::skip_on_cran()

  ## 14.
  # model includes "q" (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    for (j in 1:nobs_count) {
      if (count_type[i, j] == 1) {
        count[i, j] <- rpois(n = 1, mu[i])
      } else {
        count[i, j] <- rpois(n = 1, mu[i] * q)
      }
    }
  }

  # collect data
  data <- list(
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu
  )
  names(inits[[1]]) <- c("mu")
  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, q = TRUE, n_chain = 1, multicore = FALSE,
                      seed = 10, initial_values = inits, n_warmup = 25,
                      n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("q[1]") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - gamma, q", {
  testthat::skip_on_cran()

  ## 15.
  # model includes "q", "alpha_gamma", "beta_gamma" (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count / 2),
                      matrix(2, nrow = nsite, ncol = nobs_count / 2))

  # collect data
  data <- list(
    count = count,
    count_type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha = mu,
    beta = rep(1, length(mu)),
    q = q
  )
  names(inits[[1]]) <- c("alpha", "beta", "q")
  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, q = TRUE, family = "gamma",
                      n_chain = 1, multicore = FALSE, seed = 10,
                      initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("q[1]") %in% output_params))

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 3)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional_1",
                                                 "n_traditional_2")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2]),
                  is.numeric(out[, 3])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - phi", {
  testthat::skip_on_cran()


  ## 16.
  # model includes "phi" (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  phi <- 1.2

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)

  }

  # collect data
  data <- list(
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    phi = phi
  )
  names(inits[[1]]) <- c("mu", "phi")
  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, family = "negbin",
                      n_chain = 1, multicore = FALSE, seed = 10,
                      initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(c("phi") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 1)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 2)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - trad pois", {
  testthat::skip_on_cran()


  ## 17.
  # model, pois (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rpois(n = nobs_count, mu[i])

  }

  # collect data
  data <- list(
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu
  )
  names(inits[[1]]) <- c("mu")
  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, n_chain = 1, multicore = FALSE, seed = 10,
                      initial_values = inits, n_warmup = 25, n_iter = 100)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(!c("p10", "alpha[1]", "q", "phi") %in% output_params))

  # test expectation
  expect_true(fit$model@par_dims$phi == 0)

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 2)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

})

test_that("joint_summarize outputs work - trad, gamma", {
  testthat::skip_on_cran()

  ## 18.
  # model, gamma (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for (i in 1:nsite) {
    count[i, ] <- rgamma(nobs_count, shape = alpha_gamma[i], rate = beta_gamma)
  }

  # collect data
  data <- list(
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha = mu,
    beta = rep(1, length(mu))
  )
  names(inits[[1]]) <- c("alpha", "beta")

  # run model
  fit <- suppressWarnings({
    traditional_model(data = data, n_chain = 1, family = "gamma",
                      multicore = FALSE, seed = 10, n_warmup = 25,
                      n_iter = 100, initial_values = inits)
  })

  # get output params
  output_params <- rownames(as.data.frame(joint_summarize(fit$model)))

  # test expectation
  expect_true(all(!c("p10", "q", "phi") %in% output_params))

  # detection_calculate and detection_plot
  out <- detection_calculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10, 2)))

  # test names
  expect_true(all(names(as.data.frame(out)) == c("mu", "n_traditional")))

  # test numeric
  expect_true(all(is.numeric(out[, 1]),
                  is.numeric(out[, 2])), TRUE)

  # test plot
  out_plot <- detection_plot(fit$model, mu_min = 0.1,
                             mu_max = 1)

  # test plot type
  expect_equal(mode(out_plot), "list")

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


})
