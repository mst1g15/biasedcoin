#non-myopic logistic regression for weighted L-optimal design



#' Assuming the currrent response, future covariate value and future treatment,
#' Calculate optimality if horizon is 1. If horizon >1, iterate back to cr.exp.loss function.
#'
#' @param z.next vector of covariate values for future unit
#' @param t.next treatment of future unit
#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param z.probs vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param k number of covariates
#' @param beta vector of current estimates for regression coefficients
#' @param y vectir of responses generated up until current unit
#' @param cr  matrix of contrasts
#' @param wr matrix of weights
#' @param prior.scale prior scale parameter
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed
#' @param dyn set to NULL of there are no dynamic covariates, set to T if there are dynamic covariates
#' @param Nprobs a counter to be used if there are dynamic covariates
#'
#' @return optimality assuming current response, future covariate value and future treatment
#'
#'
#' @export
#'

cr.future.loss <- function(z.next, t.next, z.probs, N, design, k,  beta, y, cr, wr, prior.scale, bayes, dyn=NULL, Nprobs=NULL, ...){

    design.next <- rbind(design, c(1, as.numeric(z.next), t.next, as.numeric(z.next*t.next)))

  #if horizon is 1, calculate loss. else iterate.
  if (N==1){

   var <- rep(0, r)

    Imat <-  solve( Imat.beta(design.next,  beta) + diag(1/prior.scale, ncol(design))) #Optimality criterion

    for (j in 1:r){
      var[j] <-   as.matrix(wr[j]%*% t(cr[j,]) %*% Imat  %*% cr[j,])
    }

    loss <-  sum(var)

  } else{  if (!is.null(dyn)){

    dyn <- dyn+1
    z.probs <- Nprobs[dyn, ]
  }

    loss <- cr.exp.loss( z.now=z.next, t.now=t.next, z.probs, N-1, design, k, beta, y, cr, wr, dyn, bayes, prior.scale, Nprobs, ...)



  }

  return(loss)

}



#' Assuming a current response, break down the expected future optimality by cases for every combination of:
#' 1) future possible covariate
#' 2) future possible treatment
#' Find a weighted average across all these cases

#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param z.probs vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param k number of covariates
#' @param beta vector of current estimates for regression coefficients
#' @param y vectir of responses generated up until current unit
#' @param cr  matrix of contrasts
#' @param wr matrix of weights
#' @param prior.scale prior scale parameter
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed
#' @param dyn set to NULL of there are no dynamic covariates, set to T if there are dynamic covariates
#' @param Nprobs a counter to be used if there are dynamic covariates
#'
#'
#' @return expected optimality assuming a current response
#'
#'
#' @export

cr.future.y <- function(y.now, z.now, t.now, z.probs, N, design, k, beta, y, cr, wr,  prior.scale, bayes, dyn=NULL, Nprobs=NULL, ...){


  d.now <- c(1, z.now, t.now, z.now*t.now)
  row.names(design) <-NULL

  design <- rbind(design, d.now)

  #append y.now to seen data
  sim.y <- c(y, y.now)

  #given y.now, calculate new beta
  if(bayes==T){
    sim.beta <- coef(bayesglm(sim.y~design[,-1], family=binomial(link="logit")))
  }else{
    sim.beta <- coef(glm(sim.y~design[,-1], family=binomial(link="logit")))
  }

  all.z <- expand.grid(rep(list(c(-1,1)),k))  #grid for all possible covariates
  names(all.z) <- NULL

  #for each possible covariate value, calculate loss when t.next=1
  loss.p <- apply(all.z, 1, cr.future.loss , t.next=1, z.probs, N, design, k, sim.beta, sim.y, cr, wr, prior.scale, bayes, dyn, Nprobs, ...)#loss for each covariate combination when the future patient has tmt=1

  #for each possible covariate value, calculate loss when t.next=1
  loss.m <- apply(all.z, 1, cr.future.loss , t.next=-1, z.probs, N, design, k,  sim.beta, sim.y, cr, wr, prior.scale, bayes,  dyn, Nprobs, ...) #loss for each covariate combination when the future patient has tmt=-1

  #find loss for optimal treatment (t*) for each possible covariate value
  loss <- ifelse(loss.p < loss.m, loss.p, loss.m)

  #expected loss: weighed according to probabilities of each covariate value
  exploss <- sum(loss*z.probs)


  return(exploss)

}

#' Break down expected future optimality into two components:
#' 1) assuming that the current response is 0
#' 2) assuming that the current response is 1
#' Find the weighted average of the two cases

#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param z.probs vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param k number of covariates
#' @param beta vector of current estimates for regression coefficients
#' @param y vectir of responses generated up until current unit
#' @param cr  matrix of contrasts
#' @param wr matrix of weights
#' @param prior.scale prior scale parameter
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed
#' @param dyn set to NULL of there are no dynamic covariates, set to T if there are dynamic covariates
#' @param Nprobs a counter to be used if there are dynamic covariates

#'
#' @return expected optimality one step ahead in the future
#'
#'
#' @export
#'
cr.exp.loss <- function(z.now, t.now, z.probs, N, design, k, beta, y, cr, wr,  prior.scale, bayes, dyn=NULL, Nprobs=NULL,...){

  #probability that y.now=1

    pi <- probi(as.numeric(c(1, z.now, t.now, z.now*t.now)), beta)

  #expected loss, given y.now=1
  loss1 <- cr.future.y(y.now=1, z.now, t.now, z.probs, N, design, k, beta, y,  cr, wr,  prior.scale, bayes, dyn, Nprobs,  ...)

  #expected loss, given y.now=1
  loss0 <- cr.future.y(y.now=0, z.now, t.now, z.probs, N, design, k,  beta, y, cr, wr,  prior.scale, bayes, dyn, Nprobs,  ...)

  #weighted by probability that P(y=1)
  exp.loss<- pi*loss1 + (1-pi)*loss0


  return(exp.loss)
}






#' Allocate treatments according to a weighted L-optimal criterion allowing for a non-myopic approach.
#' We assume a logistic model for the response and simulate responses sequentially.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param z.probs vector of probabilities for each level of covariate z
#' @param N natural number greater than 0 for horizon
#' @param prior.scale the prior scale parameter
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.

#'
#' @return Design matrix D, all estimates of beta, final estimate of beta, responses y
#'
#'
#' @export
logit.Lbnon <- function(covar,  true.beta, threshold, kappa, init, z.probs, N,  prior.scale=100, stoc=T, bayes=T, ...){

  n <- nrow(covar)
  k <- ncol(covar) #covar must be a dataframe



  cr1 <- matrix(0, 2^k, k+1)
  cr2 <- expand.grid(rep(list(c(0,1)),k))
  cr <- as.matrix(cbind(cr1, 1, cr2))


  r <- nrow(cr)


  # randomly select treatment for first unit
  design <- as.matrix(cbind(1, covar=covar[1:init, ], tmt=sample(c(0, 1), init, replace=T)))
  design <- cbind(design, design[,2:(k+1)]*design[,(k+2)])   #append interaction columns
  row.names(design)<- NULL

  pi <- apply(design, 1, probi, true.beta)
  y <- as.numeric(rbinom(init, 1, pi))   #generate first observations based on true beta

  #find initial estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
    beta <- coef(bayesglm (y~design[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm (y~design[,-1], family=binomial(link="logit")))

  }
  all.beta <- beta


  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(design, beta) +diag(1/prior.scale, ncol(design)))%*% t(cr))))
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr

  for (i in (init+1):n){

    if (z.probs[1]=="learn"){
      z.probs <- learn.zprobs(design=design, all.z = expand.grid(rep(list(c(-1,1)),k)) , k)
    }


    #allocate treatment which minimizes expected loss
    eloss.p <- cr.exp.loss(z.now=as.numeric(covar[i,]), t.now=1, z.probs, N, design, k, beta, y, cr, wr, prior.scale, bayes,  ...) #expected loss for treatment 1
    eloss.m <- cr.exp.loss(z.now=as.numeric(covar[i,]), t.now=-1, z.probs, N, design, k, beta, y, cr, wr,  prior.scale, bayes, ...)#expected loss for treatment -1

    probs <- eloss.m/(eloss.p+eloss.m)
    if (is.na(p)){
      probs <- 0.5

    }

    if(stoc==T){
      opt.tmt <- sample(c(-1,1), 1, prob=c(probs, 1-probs))         #Assign treatments
    }else{
      if (eloss.p > eloss.m) {
        new.tmt <- -1
      } else if (eloss.p < eloss.m) {
        new.tmt <- 1
      } else if (eloss.p == eloss.m) {
        new.tmt <- sample(c(-1,1), 1)
      }
    }

    #new row for design matrix
    new.d <- as.numeric(cbind(1, data.frame(covar[i,]), opt.tmt))
    #new row if we allow for tmt-cov interaction:
    if (!is.null(int)){
      new.d <- c(new.d, new.d[2:(k+1)]*opt.tmt)   #append interaction columns
    }

    design <- as.matrix(rbind(design, as.numeric(new.d)))     #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    y <- c( y, rbinom(1, 1, pi))   #Simulate new observation

    if(bayes==T){
      beta <- coef(bayesglm (y~design[,-1], family=binomial(link="logit")))  #update beta

    }else{
      beta <- coef(glm (y~design[,-1], family=binomial(link="logit")))  #update beta

    }
    all.beta <- rbind(all.beta, beta)                          #Store all betas

  }

  design <- data.frame(design)

  results <- list(D=design, y=y, all.beta=all.beta, beta = beta)

  return(results)

}


