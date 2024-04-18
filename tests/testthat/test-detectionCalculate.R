test_that("detectionCalculate input checks work", {
  #' @srrstats {G5.2,G5.2b} Tests the assure function input checks are behaving
  #'   as expected.
  # run joint model to do tests with
  model1 <- jointModel(data=gobyData, cov=c('Filter_time','Salinity'),
                       multicore=FALSE)

  model2 <- jointModel(data=greencrabData,family='negbin',multicore=FALSE)

  #1. make sure model fit is of class stanfit
  expect_error(detectionCalculate(as.matrix(model1$model), mu = c(0.1, 0.5),
                                  cov.val = c(0,0)),
               "modelfit must be of class 'stanfit'.")

  #2. make sure ci is valid
  expect_error(detectionCalculate(model1$model, mu = c('0', 0.5),
                                  cov.val = c(0,0)),
               "mu must be a numeric vector of positive values")

  #3. make sure mu is a numeric vector of positive values
  expect_error(detectionCalculate(model1$model, mu = c(0, 0.5),
                                  cov.val = c(0,0)),
               "mu must be a numeric vector of positive values")

  #4. make sure probability is a numeric value between 0 and 1
  expect_error(detectionCalculate(model1$model, mu = c(0.1, 0.5),
                                  cov.val = c(0,0), probability = 1.05),
               "probability must be a numeric value between 0 and 1")

  #5. cov.val is numeric, if provided
  expect_error(detectionCalculate(model1$model, mu = c(0.1, 0.5),
                                  cov.val = c('0',0)),
               "cov.val must be a numeric vector")

  #6. Only include input cov.val if covariates are included in model
  expect_error(detectionCalculate(model2$model, mu = c(0.1, 0.5),
                                  cov.val = c(0,0)),
               paste0("cov.val must be 'None' if the model does not ",
                      "contain site-level covariates."))

  #7. Input cov.val is the same length as the number of estimated covariates.
  expect_error(detectionCalculate(model1$model, mu = c(0.1, 0.5),
                                  cov.val = c(0,0,0)),
               paste0("cov.val must be of the same length as the number of ",
                      "non-intercept site-level coefficients in the model."))

  #8. If covariates are in model, cov.val must be provided
  expect_error(detectionCalculate(model1$model, mu = c(0.1, 0.5)),
               paste0("cov.val must be provided if the model contains ",
                      "site-level covariates."))

  #9. qPCR.N must be an integer
  expect_error(detectionCalculate(model1$model, mu = c(0.1, 0.5),
                                  cov.val = c(0,0), qPCR.N = 6.8),
               "qPCR.N should be an integer.")
})


test_that("detectionCalculate and detectionPlot output checks", {

  ## 1.
  # model includes 'p10','q','phi','beta'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  phi <- 1.2
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n=1,mu=mu[i],size=phi)
      } else {
        count[i,j] <- rnbinom(n=1,mu=mu[i]*q,size=phi)
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    beta = beta,
    phi = phi,
    q = q
  )
  names(inits[[1]]) <- c('mu','p10','beta','phi','q')

  # run model
  fit <- jointModel(data=data, q=TRUE, family = 'negbin',
                    initial_values = inits,
                    n.chain=1, multicore=FALSE, seed = 10#,
                    #adapt_delta = 0.99
  )

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
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
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  beta <- 0.5
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(1,mu[i])
      } else {
        count[i,j] <- rpois(1,mu[i]*q)
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    beta = beta,
    q = q
  )
  names(inits[[1]]) <- c('mu','p10','beta','q')
  # run model
  fit <- jointModel(data=data, q=TRUE,
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values = inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)



  ## 3.
  # model includes 'p10','phi','beta'

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  beta <- 0.5
  log_p10 <- -4.5
  phi <- 1.2
  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n=nobs_count,mu=mu[i],size=phi)
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    beta = beta,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','beta','phi')
  # run model
  fit <- jointModel(data=data, family = 'negbin', initial_values=inits,
                    n.chain=1, multicore=FALSE, seed = 10,
  )

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)



  ## 4.
  # model includes 'p10','beta'

  # constants
  nsite <- 50
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  beta <- 0.5
  log_p10 <- -4.5

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    beta = beta
  )
  names(inits[[1]]) <- c('mu','p10','beta')

  # run model
  fit <- jointModel(data=data, initial_values = inits,
                    n.chain=1, multicore=FALSE, seed = 10)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 5.
  # model includes 'p10','q','phi','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  phi <- 10
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n=1,mu=mu[i],size=phi)
      } else {
        count[i,j] <- rnbinom(n=1,mu=mu[i]*q,size=phi)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    alpha = alpha,
    q = q,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha','q','phi'
  )
  # run model
  fit <- jointModel(data=data, family = 'negbin', q=TRUE,
                    cov=c('var_a','var_b'),#n.iter.burn=5000,
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values=inits, adapt_delta = 0.99,
  )

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1),
                            cov.val=c(0,0))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 6.
  # model includes 'p10','alpha','q'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(n=1,mu[i])
      } else {
        count[i,j] <- rpois(n=1,mu[i]*q)
      }
    }
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
  )
  # run model
  fit <- jointModel(data=data, q=TRUE,
                    cov=c('var_a','var_b'),
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1),
                            cov.val=c(0,0))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)



  ## 7.
  # model includes 'p10','phi','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  phi <- 10

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n=nobs_count,mu=mu[i],size=phi)
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    alpha = alpha,
    phi = phi
  )
  names(inits[[1]]) <- c('mu','p10','alpha',
                         'phi')
  # run model
  fit <- jointModel(data=data, family = 'negbin',
                    cov=c('var_a','var_b'),
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1),
                            cov.val=c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 8.
  # model includes 'p10','alpha'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(nobs_count,mu[i])
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    alpha = alpha
  )
  names(inits[[1]]) <- c('mu','p10','alpha'
  )
  # run model
  fit <- jointModel(data=data,
                    cov=c('var_a','var_b'),
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1),
                            cov.val=c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 9.
  # model includes 'q','phi' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  phi <- 1.2
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rnbinom(n=1,mu=mu[i],size=phi)
      } else {
        count[i,j] <- rnbinom(n=1,mu=mu[i]*q,size=phi)
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
  fit <- traditionalModel(data=data, q=TRUE, family='negbin',
                          n.chain=1, multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 10.
  # model includes 'q' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  q <- 2
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    for(j in 1:nobs_count){
      if(count_type[i,j]==1){
        count[i,j] <- rpois(n=1,mu[i])
      } else {
        count[i,j] <- rpois(n=1,mu[i]*q)
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
  fit <- traditionalModel(data=data, q=TRUE,
                          n.chain=1, multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

  ## 11.
  # model includes 'phi' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  phi <- 1.2

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rnbinom(n=nobs_count,mu=mu[i],size=phi)

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
  fit <- traditionalModel(data=data,family='negbin',
                          n.chain=1, multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 12.
  # model, pois (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rpois(n=nobs_count,mu[i])

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
  fit <- traditionalModel(data=data,n.chain=1, multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 13.
  # model, gamma (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape=alpha_gamma[i],rate=beta_gamma)
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
  fit <- traditionalModel(data=data,n.chain=1, family='gamma',
                          multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

  # test dimensions
  expect_true(all(dim(out) == c(10,2)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)

  ## 14.
  # model includes 'q','alpha_gamma','beta_gamma' (traditional model)

  # constants
  nsite <- 20
  nobs_count <- 100
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  q <- 2
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma
  # traditional type
  count_type <- cbind(matrix(1,nrow=nsite,ncol=nobs_count/2),
                      matrix(2,nrow=nsite,ncol=nobs_count/2))

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
    q=q
  )
  names(inits[[1]]) <- c('alpha','beta','q')
  # run model
  fit <- traditionalModel(data=data, q=TRUE,family='gamma',
                          n.chain=1, multicore=FALSE, seed = 10,
                          initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1))

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
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1)

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)


  ## 15.
  # model includes 'p10','alpha','alpha_gamma','beta_gamma'

  # constants
  nsite <- 20
  nobs_count <- 100
  nobs_pcr <- 8
  # params
  mu <- rlnorm(nsite,meanlog=log(1),sdlog=1)
  alpha <- c(0.5, 0.1, -0.4)
  log_p10 <- -4.5
  beta_gamma <- 1
  alpha_gamma <- mu * beta_gamma

  # count
  count <- matrix(NA,nrow=nsite,ncol=nobs_count)
  for(i in 1:nsite){
    count[i,] <- rgamma(nobs_count,shape=alpha_gamma[i],rate=beta_gamma)
  }
  # site-level covariates
  mat_site <- matrix(NA,nrow=nsite,ncol=length(alpha))
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
  qPCR.N <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
  for(i in 1:nsite){
    qPCR.N[i,] <- rep(3,nobs_pcr)
  }
  # qPCR.K (# positive qPCR detections)
  qPCR.K <- matrix(NA,nrow=nsite,ncol=nobs_pcr)
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
    p10 = log_p10,
    alpha = alpha
  )
  names(inits[[1]]) <- c('alpha_gamma','beta_gamma','p10','alpha'
  )
  # run model
  fit <- jointModel(data=data,family='gamma',
                    cov=c('var_a','var_b'),
                    n.chain=1, multicore=FALSE, seed = 10,
                    initial_values=inits)

  out <- detectionCalculate(fit$model, mu=seq(from=0.1, to=1, by=0.1),
                            cov.val=c(0,0))

  # test dimensions
  expect_true(all(dim(out) == c(10,3)))

  # test names
  expect_true(all(names(as.data.frame(out))==c('mu','n_traditional','n_eDNA')))

  # test numeric
  expect_true(all(is.numeric(out[,1]),
                  is.numeric(out[,2]),
                  is.numeric(out[,3])),TRUE)

  # test plot
  out_plot <- detectionPlot(fit$model, mu.min=0.1,
                            mu.max = 1, cov.val = c(0,0))

  # test plot type
  expect_equal(mode(out_plot),'list')

  # test data in plot
  expect_equal(all(is.numeric(out_plot$data$mu),
                   is.character(out_plot$data$survey_type),
                   is.numeric(out_plot$data$value)), TRUE)



})

