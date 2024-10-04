test_that("jointSummarize input checks work", {
  testthat::skip_on_cran()
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are
  #'   behaving as expected.

  # run traditional model to do tests with
  data <- data("greencrabData")

  out <- traditionalModel(data = greencrabData, family = 'negbin',
                          multicore = FALSE,
                          n.chain = 1, n.iter.sample = 1000)

  #1. make sure model fit is of class stanfit
  data <- data.frame(y = c(1,2,3),x = c(1,2,3))
  lm_out <- lm(y ~ x, data = data)
  expect_error(jointSummarize(lm_out$model),
               "modelfit must be of class 'stanfit'.")

  #2. make sure probs is a numeric vector
  expect_error(jointSummarize(out$model,probs = c('95%')),
               "probs must be a numeric vector.")

  #3. make sure all values of probs are between 0 and 1
  expect_error(jointSummarize(out$model,probs = c(5,95)),
               "probs must be between 0 and 1.")

  #4. make sure par is a character vector
  expect_error(jointSummarize(out$model,par = c(1,2,3)),
               "par must be a character vector.")

  #5. make sure model fit contains all par input
  expect_error(jointSummarize(out$model,par = 'alpha'),
               "modelfit must contain all selected parameters: alpha")
})

test_that("jointSummarize outputs work", {
  testthat::skip_on_cran()

  ## 1.
  # model includes 'p10','q','phi','beta'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i,j] <- rnbinom(n = 1, mu = mu[i]*q, size = phi)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )

  # initial values
  inits <- list()
  inits[[1]] <- list(
  mu = mu,
  p10 = exp(log_p10),
  beta = beta,
  phi = phi,
  q = q
  )
  names(inits[[1]]) <- c('mu','p10','beta','phi','q')

  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'negbin',
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75,
                    n.chain = 1, multicore = FALSE, seed = 10
                    )})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','phi','beta') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,4)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3]),
                  is.numeric(out[,4])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 2.
  # model includes 'p10','beta','q'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(1, mu[i])
      } else {
        count[i,j] <- rpois(1, mu[i]*q)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i], nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    beta = beta,
    q = q
  )
  names(inits[[1]]) <- c('mu','p10','beta','q')
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE,
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','beta') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,4)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3]),
                  is.numeric(out[,4])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

  ## 3.
  # model includes 'p10','beta','q', 'alpha_gamma','beta_gamma'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rgamma(1,shape = alpha_gamma[i],rate = beta_gamma)
      } else {
        count[i,j] <- rgamma(1,shape = alpha_gamma[i]*q,rate = beta_gamma)
      }
    }
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    beta = beta,
    q = q
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','beta','q')
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'gamma',
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','beta') %in% output_params))

  ## 4.
  # model includes 'p10','phi','beta'

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
  for(i in 1:nsite){
    count[i,] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    beta = beta,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','beta','phi')
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin',
                                      initial_values = inits,
                    n.chain = 1, multicore = FALSE, seed = 10,
                    n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','phi','beta') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 5.
  # model includes 'p10','beta'

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
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count, mu[i])
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    beta = beta
  )
  names(inits[[1]]) <- c('mu','p10','beta')
  # run model
  fit <- suppressWarnings({jointModel(data = data, initial_values = inits,
                    n.chain = 1, multicore = FALSE, seed = 10,
                    n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','beta') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 6.
  # model includes 'p10','beta','alpha_gamma','alpha_beta'

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
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape = alpha_gamma[i],rate = beta_gamma)
  }
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(beta))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    beta = beta
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','beta')
  # run model
  fit <- suppressWarnings({jointModel(data = data, initial_values = inits,
                                      family = 'gamma',
                    n.chain = 1, multicore = FALSE, seed = 10,
                    n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','beta') %in% output_params))


  ## 7.
  # model includes 'p10','q','phi','alpha'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i,j] <- rnbinom(n = 1, mu = mu[i]*q, size = phi)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
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
  names(inits[[1]]) <- c('mu','p10','alpha','q','phi'
                         )
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin', q = TRUE,
                    cov = c('var_a','var_b'),
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, adapt_delta = 0.99,
                    n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','phi','alpha[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                            cov.val = c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,4)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3]),
                  is.numeric(out[,4])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 8.
  # model includes 'p10','alpha','q'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(n = 1, mu[i])
      } else {
        count[i,j] <- rpois(n = 1, mu[i]*q)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
                         )
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE,
                    cov = c('var_a','var_b'),
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                            cov.val = c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,4)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3]),
                  is.numeric(out[,4])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 9.
  # model includes 'p10','alpha','q', 'alpha_gamma', 'alpha_beta'

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
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rgamma(1,shape = alpha_gamma[i],rate = beta_gamma)
      } else {
        count[i,j] <- rgamma(1,shape = alpha_gamma[i]*q,rate = beta_gamma)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    count.type = count_type,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    alpha = alpha,
    q = q
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha','q'
                         )
  # run model
  fit <- suppressWarnings({jointModel(data = data, q = TRUE, family = 'gamma',
                    cov = c('var_a','var_b'),
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75
                    )})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','q[1]','alpha[1]') %in% output_params))


  ## 10.
  # model includes 'p10','phi','alpha'

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
  for(i in 1:nsite){
    count[i,] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha',
                         'phi')
  # run model
  fit <- suppressWarnings({jointModel(data = data, family = 'negbin',
                    cov = c('var_a','var_b'),n.iter.burn = 25,
                    n.iter.sample = 75,
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','phi','alpha[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                            cov.val = c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 11.
  # model includes 'p10','alpha'

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
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count, mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
                         )
  # run model
  fit <- suppressWarnings({jointModel(data = data,
                    cov = c('var_a','var_b'),
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                            cov.val = c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 12.
  # model includes 'p10','alpha','alpha_gamma','beta_gamma'

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
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape = alpha_gamma[i],rate = beta_gamma)
  }
  # site-level covariates
  mat_site <- matrix(NA, nrow = nsite, ncol = length(alpha))
  mat_site[,1] <- 1 # intercept
  for(i in 2:length(alpha)){
    mat_site[,i] <- rnorm(nsite,0,1)
  }
  colnames(mat_site) <- c('int','var_a','var_b')
  # p11 (probability of true positive eDNA detection) and p (probability
  # of eDNA detection)
  p11 <- rep(NA,nsite)
  p <- rep(NA,nsite)
  for (i in 1:nsite){
    p11[i] <- mu[i] / (mu[i] + exp(sum(mat_site[i,]*alpha)))
    p[i] <- min(p11[i] + exp(log_p10),1)
  }
  # qPCR.N (# qPCR observations)
  qPCR.N <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA, nrow = nsite, ncol = nobs_pcr)
  for (i in 1:nsite){
    qPCR.K[i,] <- rbinom(nobs_pcr, qPCR.N[i,], rep(p[i],nobs_pcr))
  }
  # collect data
  data <- list(
    qPCR.N = qPCR.N,
    qPCR.K = qPCR.K,
    count = count,
    site.cov = mat_site
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha_gamma = mu,
    beta_gamma = rep(1,length(mu)),
    p10 = exp(log_p10),
    alpha = alpha
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha'
                         )
  # run model
  fit <- suppressWarnings({jointModel(data = data,family = 'gamma',
                    cov = c('var_a','var_b'),
                    n.chain = 1, multicore = FALSE, seed = 10,
                    initial_values = inits, n.iter.burn = 25,
                    n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('p10','alpha[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1),
                            cov.val = c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 13.
  # model includes 'q','phi' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  phi <- 1.2
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n = 1, mu = mu[i], size = phi)
      } else {
        count[i,j] <- rnbinom(n = 1, mu = mu[i]*q, size = phi)
      }
    }
  }

  # collect data
  data <- list(
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','phi')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data, q = TRUE,
                                            family = 'negbin',
                                            n.chain = 1, multicore = FALSE,
                                            seed = 10,
                                            initial_values = inits,
                                            n.iter.burn = 25,
                                            n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('q[1]','phi') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 14.
  # model includes 'q' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(n = 1, mu[i])
      } else {
        count[i,j] <- rpois(n = 1, mu[i]*q)
      }
    }
  }

  # collect data
  data <- list(
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    mu = mu
  )
  names(inits[[1]]) <- c('mu')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data, q = TRUE,
                          n.chain = 1, multicore = FALSE, seed = 10,
                          initial_values = inits, n.iter.burn = 25,
                          n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('q[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

  ## 15.
  # model includes 'q','alpha_gamma','beta_gamma' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1, nrow = nsite, ncol = nobs_count/2),
                      matrix(2, nrow = nsite, ncol = nobs_count/2))

  # collect data
  data <- list(
    count = count,
    count.type = count_type
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha = mu,
    beta = rep(1,length(mu)),
    q = q
  )
  names(inits[[1]]) <- c('alpha','beta','q')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data, q = TRUE,
                                            family = 'gamma',
                                            n.chain = 1, multicore = FALSE,
                                            seed = 10,
                                            initial_values = inits,
                                            n.iter.burn = 25,
                                            n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('q[1]') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional_1',
                                               'n_traditional_2')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 16.
  # model includes 'phi' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)
  phi <- 1.2

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n = nobs_count, mu = mu[i], size = phi)

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
  names(inits[[1]]) <- c('mu','phi')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data,family = 'negbin',
                                            n.chain = 1, multicore = FALSE,
                                            seed = 10,
                                            initial_values = inits,
                                            n.iter.burn = 25,
                                            n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(c('phi') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 17.
  # model, pois (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite, meanlog = log(1), sdlog = 1)

  # count
  count <- matrix(NA, nrow = nsite, ncol = nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(n = nobs_count, mu[i])

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
  names(inits[[1]]) <- c('mu')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data,n.chain = 1,
                                            multicore = FALSE, seed = 10,
                                            initial_values = inits,
                                            n.iter.burn = 25,
                                            n.iter.sample = 75)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(!c('p10','beta','q','phi') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

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
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape = alpha_gamma[i],rate = beta_gamma)
  }

  # collect data
  data <- list(
    count = count
  )
  # initial values
  inits <- list()
  inits[[1]] <- list(
    alpha = mu,
    beta = rep(1,length(mu))
  )
  names(inits[[1]]) <- c('alpha','beta')
  # run model
  fit <- suppressWarnings({traditionalModel(data = data, n.chain = 1,
                                            family = 'gamma',
                                            multicore = FALSE, seed = 10,
                                            n.iter.burn = 25,
                                            n.iter.sample = 75,
                                            initial_values = inits)})

  # get output params
  output_params <- rownames(as.data.frame(jointSummarize(fit$model)))

  # test expectation
  expect_true(all(!c('p10','q','phi') %in% output_params))

  # detectionCalculate and detectionPlot
  out <- detectionCalculate(fit$model, mu = seq(from = 0.1, to = 1, by = 0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min = 0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


})


