#' Calculate the survey effort necessary to detect species presence, given the species expected catch rate.
#'
#' This function calculates the median number of survey effort units to necessary detect species presence. Detecting species presence is defined as producing at least one true positive eDNA detection or catching at least one individual.
#'
#' @export
#' @param modelfit An object of class `stanfit`.
#' @param mu A numeric vector of species densities/capture rates. If multiple traditional gear types are represented in the model, mu is the catch rate of gear type 1.
#' @param cov.val A numeric vector indicating the values of site-level covariates to use for prediction. Default is 'None'.
#' @param probability A numeric value indicating the probability of detecting presence. The default is 0.9.
#' @param qPCR.N An integer indicating the number of qPCR replicates per eDNA sample. The default is 3.
#' @return A summary table of survey efforts necessary to detect species presence, given mu, for each survey type.
#'
#' @note  Before fitting the model, this function checks to ensure that the function is possible given the inputs. These checks include:
#' \itemize{
#' \item  Input model fit is an object of class 'stanfit'.
#' \item  Input mu is a numeric vector.
#' \item  Input probability is a numeric value.
#' \item  If model fit contains alpha, cov.val must be provided.
#' \item  Input cov.val is numeric.
#' \item  Input cov.val is the same length as the number of estimated covariates.
#' \item  Input model fit has converged (i.e. no divergent transitions after warm-up).
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' @examples
#' \donttest{
#' # Ex. 1: Calculating necessary effort for detection with site-level covariates
#'
#' # Load data
#' data(gobyData)
#'
#' # Fit a model including 'Filter_time' and 'Hab_size' site-level covariates
#' fit.cov = jointModel(data=gobyData, cov=c('Filter_time','Hab_size'),
#'                      family='poisson', p10priors=c(1,20), q=FALSE)
#'
#' # Calculate at the mean covariate values (covariates are standardized, so mean=0)
#' detection.calculate(fit.cov, mu = seq(from=0.1,to=1,by=0.1), cov.val = c(0,0), qPCR.N = 3)
#'
#' # Calculate mu_critical at habitat size 0.5 z-scores greater than the mean
#' detection.calculate(fit.cov, mu = seq(from=0.1,to=1,by=0.1), cov.val = c(0,0.5), qPCR.N = 3)
#'
#' # Ex. 2: Calculating necessary effort for detection with multiple traditional gear types
#'
#' # Load data
#' data(greencrabData)
#'
#' # Fit a model with no site-level covariates
#' fit.q = jointModel(data=greencrabData, cov='None', family='negbin',
#'                    p10priors=c(1,20), q=TRUE)
#'
#' # Calculate
#' detection.calculate(fit.q, mu = seq(from=0.1,to=1,by=0.1), cov.val = 'None', qPCR.N = 3)
#'
#' # Change probability of detecting presence to 0.95
#' detection.calculate(fit.q, mu = 0.1, cov.val = 'None', probability = 0.95, qPCR.N = 3)
#' }
#'

detection.calculate <- function(modelfit, mu, cov.val = 'None', probability=0.9, qPCR.N = 3){

  ## #1. make sure model fit is of class stanfit
  if(!is(modelfit,'stanfit')) {
    errMsg = paste("modelfit must be of class 'stanfit'.")
    stop(errMsg)
  }

  ## #2. make sure mu is a numeric vector
  if(!is.numeric(mu)) {
    errMsg = paste("mu must be a numeric vector")
    stop(errMsg)
  }

  ## #3. make sure probability is a numeric value between 0 and 1
  if(!is.numeric(probability) | length(probability)>1 | any(probability < 0) | any(probability > 1)) {
    errMsg = paste("probability must be a numeric value between 0 and 1")
    stop(errMsg)
  }

  ## #4. cov.val is numeric, if provided
  if(all(cov.val != 'None') && !is.numeric(cov.val)) {
    errMsg = paste("cov.val must be a numeric vector")
    stop(errMsg)
  }

  ## #5. Only include input cov.val if covariates are included in model
  if(all(cov.val != 'None') && !c('alpha') %in% modelfit@model_pars) {
    errMsg = paste("cov.val must be 'None' if the model does not contain site-level covariates.")
    stop(errMsg)
  }

  ## #6. Input cov.val is the same length as the number of estimated covariates.
  if(all(cov.val != 'None') && length(cov.val)!=(modelfit@par_dims$alpha-1)) {
    errMsg = paste("cov.val must be of the same length as the number of non-intercept site-level coefficients in the model.")
    stop(errMsg)
  }

  if (!requireNamespace("rstan", quietly = TRUE)){
    stop ("The 'rstan' package is not installed.", call. = FALSE)
  }

  ## #7. check to see if there are any divergent transitions
  if(sum(rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[1]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[2]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[3]][,'divergent__'],
         rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[4]][,'divergent__']) > 0 ){

    sum <- sum(rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[1]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[2]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[3]][,'divergent__'],
               rstan::get_sampler_params(modelfit,inc_warmup=FALSE)[[4]][,'divergent__'])

    warning <- paste0('Warning: There are ',sum,' divergent transitions in your model fit. ')
    print(warning)

  }


  ############################################
  #joint, no catchability, negbin, no sitecov#
  ############################################

  if('phi' %in% modelfit@model_pars && all(!c('q','alpha') %in% modelfit@model_pars)){
    #get median parameter estimates
    phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))
    beta <- stats::median(unlist(rstan::extract(modelfit,pars='beta')))

    ##
    #number of traditional survey effort
    ntrad_out <- vector(length = length(mu))

    for(i in 1:length(mu)){
      #number of traditional survey replicates
      ntrad <- seq(from = 0, to = 50000, by = 1)

      pr <- stats::pnbinom(q = 0, mu = mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
      #dnbinom: x = failed events; size = successful events
      prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

      #find value
      value <- match(min(prob[prob >= probability]), prob)
      ntrad_out[i] <- ntrad[value]
      }

    ##
    #number of qPCR reactions to get at least one true detection (among N replicates)
    ndna_out <- vector(length = length(mu))

    for(i in 1:length(mu)){
      #number of eDNA samples
      ndna <- seq(from = 0, to = 50000, by = 1)

      p11 <- mu[i]/(mu[i]+exp(beta))
      #P(at least one detection|mu) out of total bottles/PCR replicates
      prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

      #find value
      value <- match(min(prob[prob >= probability]), prob)
      ndna_out[i] <- ndna[value]
    }

    out <- cbind(mu,ntrad_out,ndna_out)
    colnames(out) <- c('mu','n_traditional','n_eDNA')

    ##########################################
    #joint, no catchability, pois, no sitecov#
    ##########################################
    } else if(all(!c('q','phi','alpha') %in% modelfit@model_pars)){
      #get median parameter estimates
      beta <- stats::median(unlist(rstan::extract(modelfit,pars='beta')))

      ##
      #number of traditional survey effort
      ntrad_out <- vector(length = length(mu))

      for(i in 1:length(mu)){
        #number of traditional survey replicates
        ntrad <- seq(from = 0, to = 50000, by = 1)

        pr <- stats::ppois(q = 0, lambda = mu[i]) #P(X = 0 | mu) in one traditional survey trial
        #dnbinom: x = failed events; size = successful events
        prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

        #find value
        value <- match(min(prob[prob >= probability]), prob)
        ntrad_out[i] <- ntrad[value]
      }

      ##
      #number of qPCR reactions to get at least one true detection (among N replicates)
      ndna_out <- vector(length = length(mu))

      for(i in 1:length(mu)){
        #number of eDNA samples
        ndna <- seq(from = 0, to = 50000, by = 1)

        p11 <- mu[i]/(mu[i]+exp(beta))
        #P(at least one detection|mu) out of total bottles/PCR replicates
        prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

        #find value
        value <- match(min(prob[prob >= probability]), prob)
        ndna_out[i] <- ndna[value]
      }

      out <- cbind(mu,ntrad_out,ndna_out)
      colnames(out) <- c('mu','n_traditional','n_eDNA')

      #########################################
      #joint, catchability, negbin, no sitecov#
      #########################################
    } else if(all(c('phi','q') %in% modelfit@model_pars) && !c('alpha') %in% modelfit@model_pars){
      #get median parameter estimates
      phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))
      beta <- stats::median(unlist(rstan::extract(modelfit,pars='beta')))

      #create empty q list
      q_list <- list()

      ##fill in q list
      for(i in 1:modelfit@par_dims$q){
       name <- paste('q',i,sep='')
        q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
      }
      q_list <- c(1,unlist(q_list))

      ##
      #number of traditional survey effort
      ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

      for(j in 1:length(q_list)){
        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::pnbinom(q = 0, mu = q_list[j]*mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i,j] <- ntrad[value]
        }
      }

      ##
      #number of qPCR reactions to get at least one true detection (among N replicates)
      ndna_out <- vector(length = length(mu))

      for(i in 1:length(mu)){
        #number of eDNA samples
        ndna <- seq(from = 0, to = 50000, by = 1)

        p11 <- mu[i]/(mu[i]+exp(beta))
        #P(at least one detection|mu) out of total bottles/PCR replicates
        prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

        #find value
        value <- match(min(prob[prob >= probability]), prob)
        ndna_out[i] <- ndna[value]
      }

      out <- cbind(mu,ntrad_out,ndna_out)
      #rename columns
      for(i in 1:modelfit@par_dims$q){
        trad_names <- paste('n_traditional_',i+1,sep='')
      }
      colnames(out) <- c('mu','n_traditional_1',trad_names,'n_eDNA')

      #######################################
      #joint, catchability, pois, no sitecov#
      #######################################
    } else if('q' %in% modelfit@model_pars && all(!c('phi','alpha') %in% modelfit@model_pars)){
      #get median parameter estimates
      beta <- stats::median(unlist(rstan::extract(modelfit,pars='beta')))

      #create empty q list
      q_list <- list()

      ##fill in q list
      for(i in 1:modelfit@par_dims$q){
        name <- paste('q',i,sep='')
        q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
      }
      q_list <- c(1,unlist(q_list))

      ##
      #number of traditional survey effort
      ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

      for(j in 1:length(q_list)){
        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::ppois(q = 0, lambda = q_list[j]*mu[i]) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i,j] <- ntrad[value]
        }
      }

      ##
      #number of qPCR reactions to get at least one true detection (among N replicates)
      ndna_out <- vector(length = length(mu))

      for(i in 1:length(mu)){
        #number of eDNA samples
        ndna <- seq(from = 0, to = 50000, by = 1)

        p11 <- mu[i]/(mu[i]+exp(beta))
        #P(at least one detection|mu) out of total bottles/PCR replicates
        prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

        #find value
        value <- match(min(prob[prob >= probability]), prob)
        ndna_out[i] <- ndna[value]
      }

      out <- cbind(mu,ntrad_out,ndna_out)
      #rename columns
      for(i in 1:modelfit@par_dims$q){
        trad_names <- paste('n_traditional_',i+1,sep='')
      }
      colnames(out) <- c('mu','n_traditional_1',trad_names,'n_eDNA')

      #########################################
      #joint, no catchability, negbin, sitecov#
      #########################################

      } else if(all(c('phi','alpha') %in% modelfit@model_pars) && !c('q') %in% modelfit@model_pars){
        #get median parameter estimates
        phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))
        alpha <- apply(rstan::extract(modelfit,pars='alpha')$alpha,2,'median')
        beta <- alpha %*% c(1,cov.val)

        ##
        #number of traditional survey effort
        ntrad_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::pnbinom(q = 0, mu = mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i] <- ntrad[value]
        }

        ##
        #number of qPCR reactions to get at least one true detection (among N replicates)
        ndna_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of eDNA samples
          ndna <- seq(from = 0, to = 50000, by = 1)

          p11 <- mu[i]/(mu[i]+exp(beta))
          #P(at least one detection|mu) out of total bottles/PCR replicates
          prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ndna_out[i] <- ndna[value]
        }

        out <- cbind(mu,ntrad_out,ndna_out)
        colnames(out) <- c('mu','n_traditional','n_eDNA')

        #######################################
        #joint, no catchability, pois, sitecov#
        #######################################
      } else if('alpha' %in% modelfit@model_pars && all(!c('q','phi') %in% modelfit@model_pars)){
        #get median parameter estimates
        alpha <- apply(rstan::extract(modelfit,pars='alpha')$alpha,2,'median')
        beta <- alpha %*% c(1,cov.val)


        ##
        #number of traditional survey effort
        ntrad_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::ppois(q = 0, lambda = mu[i]) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i] <- ntrad[value]
        }

        ##
        #number of qPCR reactions to get at least one true detection (among N replicates)
        ndna_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of eDNA samples
          ndna <- seq(from = 0, to = 50000, by = 1)

          p11 <- mu[i]/(mu[i]+exp(beta))
          #P(at least one detection|mu) out of total bottles/PCR replicates
          prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ndna_out[i] <- ndna[value]
        }

        out <- cbind(mu,ntrad_out,ndna_out)
        colnames(out) <- c('mu','n_traditional','n_eDNA')

        ######################################
        #joint, catchability, negbin, sitecov#
        ######################################
      } else if(all(c('phi','q','alpha') %in% modelfit@model_pars)){
        #get median parameter estimates
        phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))
        alpha <- apply(rstan::extract(modelfit,pars='alpha')$alpha,2,'median')
        beta <- alpha %*% c(1,cov.val)

        #create empty q list
        q_list <- list()

        ##fill in q list
        for(i in 1:modelfit@par_dims$q){
          name <- paste('q',i,sep='')
          q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
        }
        q_list <- c(1,unlist(q_list))

        ##
        #number of traditional survey effort
        ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

        for(j in 1:length(q_list)){
          for(i in 1:length(mu)){
            #number of traditional survey replicates
            ntrad <- seq(from = 0, to = 50000, by = 1)

            pr <- stats::pnbinom(q = 0, mu = q_list[j]*mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
            #dnbinom: x = failed events; size = successful events
            prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

            #find value
            value <- match(min(prob[prob >= probability]), prob)
            ntrad_out[i,j] <- ntrad[value]
          }
        }

        ##
        #number of qPCR reactions to get at least one true detection (among N replicates)
        ndna_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of eDNA samples
          ndna <- seq(from = 0, to = 50000, by = 1)

          p11 <- mu[i]/(mu[i]+exp(beta))
          #P(at least one detection|mu) out of total bottles/PCR replicates
          prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ndna_out[i] <- ndna[value]
        }

        out <- cbind(mu,ntrad_out,ndna_out)
        #rename columns
        for(i in 1:modelfit@par_dims$q){
          trad_names <- paste('n_traditional_',i+1,sep='')
        }
        colnames(out) <- c('mu','n_traditional_1',trad_names,'n_eDNA')

        #######################################
        #joint, catchability, pois, no sitecov#
        #######################################
      } else if(all(c('alpha','q') %in% modelfit@model_pars) && !c('phi') %in% modelfit@model_pars){
        #get median parameter estimates
        alpha <- apply(rstan::extract(modelfit,pars='alpha')$alpha,2,'median')
        beta <- alpha %*% c(1,cov.val)

        #create empty q list
        q_list <- list()

        ##fill in q list
        for(i in 1:modelfit@par_dims$q){
          name <- paste('q',i,sep='')
          q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
        }
        q_list <- c(1,unlist(q_list))

        ##
        #number of traditional survey effort
        ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

        for(j in 1:length(q_list)){
          for(i in 1:length(mu)){
            #number of traditional survey replicates
            ntrad <- seq(from = 0, to = 50000, by = 1)

            pr <- stats::ppois(q = 0, lambda = q_list[j]*mu[i]) #P(X = 0 | mu) in one traditional survey trial
            #dnbinom: x = failed events; size = successful events
            prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

            #find value
            value <- match(min(prob[prob >= probability]), prob)
            ntrad_out[i,j] <- ntrad[value]
          }
        }

        ##
        #number of qPCR reactions to get at least one true detection (among N replicates)
        ndna_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of eDNA samples
          ndna <- seq(from = 0, to = 50000, by = 1)

          p11 <- mu[i]/(mu[i]+exp(beta))
          #P(at least one detection|mu) out of total bottles/PCR replicates
          prob <- 1 - stats::dbinom(x = 0, size = ndna*qPCR.N, prob = p11)

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ndna_out[i] <- ndna[value]
        }

        out <- cbind(mu,ntrad_out,ndna_out)
        #rename columns
        for(i in 1:modelfit@par_dims$q){
          trad_names <- paste('n_traditional_',i+1,sep='')
        }
        colnames(out) <- c('mu','n_traditional_1',trad_names,'n_eDNA')

        ###############################
        #trad, no catchability, negbin#
        ###############################

      } else if(c('phi') %in% modelfit@model_pars && all(!c('q','p10') %in% modelfit@model_pars)){
        #get median parameter estimates
        phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))

        ##
        #number of traditional survey effort
        ntrad_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::pnbinom(q = 0, mu = mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i] <- ntrad[value]
        }

        out <- cbind(mu,ntrad_out)
        colnames(out) <- c('mu','n_traditional')

        #############################
        #trad, no catchability, pois#
        #############################
      } else if(all(!c('q','p10','phi') %in% modelfit@model_pars)){

        ##
        #number of traditional survey effort
        ntrad_out <- vector(length = length(mu))

        for(i in 1:length(mu)){
          #number of traditional survey replicates
          ntrad <- seq(from = 0, to = 50000, by = 1)

          pr <- stats::ppois(q = 0, lambda = mu[i]) #P(X = 0 | mu) in one traditional survey trial
          #dnbinom: x = failed events; size = successful events
          prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

          #find value
          value <- match(min(prob[prob >= probability]), prob)
          ntrad_out[i] <- ntrad[value]
        }

        out <- cbind(mu,ntrad_out)
        colnames(out) <- c('mu','n_traditional')

        ############################
        #trad, catchability, negbin#
        ############################
      } else if(all(c('phi','q') %in% modelfit@model_pars) && !c('p10') %in% modelfit@model_pars){
        #get median parameter estimates
        phi <- stats::median(unlist(rstan::extract(modelfit,pars='phi')))

        #create empty q list
        q_list <- list()

        ##fill in q list
        for(i in 1:modelfit@par_dims$q){
          name <- paste('q',i,sep='')
          q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
        }
        q_list <- c(1,unlist(q_list))

        ##
        #number of traditional survey effort
        ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

        for(j in 1:length(q_list)){
          for(i in 1:length(mu)){
            #number of traditional survey replicates
            ntrad <- seq(from = 0, to = 50000, by = 1)

            pr <- stats::pnbinom(q = 0, mu = q_list[j]*mu[i], size = phi) #P(X = 0 | mu) in one traditional survey trial
            #dnbinom: x = failed events; size = successful events
            prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

            #find value
            value <- match(min(prob[prob >= probability]), prob)
            ntrad_out[i,j] <- ntrad[value]
          }
        }

        out <- cbind(mu,ntrad_out)
        #rename columns
        for(i in 1:modelfit@par_dims$q){
          trad_names <- paste('n_traditional_',i+1,sep='')
        }
        colnames(out) <- c('mu','n_traditional_1',trad_names)

        ##########################
        #trad, catchability, pois#
        ##########################
      } else if(c('q') %in% modelfit@model_pars && all(!c('phi','p10') %in% modelfit@model_pars)){

        #create empty q list
        q_list <- list()

        ##fill in q list
        for(i in 1:modelfit@par_dims$q){
          name <- paste('q',i,sep='')
          q_list[[name]] <- stats::median(rstan::extract(modelfit,pars='q')$q[,i])
        }
        q_list <- c(1,unlist(q_list))

        ##
        #number of traditional survey effort
        ntrad_out <- matrix(NA,nrow = length(mu),ncol = length(q_list))

        for(j in 1:length(q_list)){
          for(i in 1:length(mu)){
            #number of traditional survey replicates
            ntrad <- seq(from = 0, to = 50000, by = 1)

            pr <- stats::ppois(q = 0, lambda = q_list[j]*mu[i]) #P(X = 0 | mu) in one traditional survey trial
            #dnbinom: x = failed events; size = successful events
            prob <- 1 - stats::dnbinom(x = 0, size = ntrad, prob = pr) #probability of catching >0 animals based on ntrad

            #find value
            value <- match(min(prob[prob >= probability]), prob)
            ntrad_out[i,j] <- ntrad[value]
          }
        }

        out <- cbind(mu,ntrad_out)
        #rename columns
        for(i in 1:modelfit@par_dims$q){
          trad_names <- paste('n_traditional_',i+1,sep='')
        }
        colnames(out) <- c('mu','n_traditional_1',trad_names)
      }

    return(out)
  }

