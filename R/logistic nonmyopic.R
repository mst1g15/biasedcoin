
#############################################################################################
#non-myopic logistic regression




#' Assuming the currrent response, future covariate value and future treatment,
#' Calculate optimality if horizon is 1. If not, iterate back to exp.loss function.
#'
#' @param z.next vector of covariate values for future unit
#' @param t.next treatment of future unit
#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param zp vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param int set to NULL if there are no interactions, set to T of there are interactions
#' @param lossfunc  a function for the optimality criterion to minimize
#' @param beta estimate of the regression coefficients
#' @param y responses that have been observed up until the current unit
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.

#' @return optimality assuming current response, future covariate value and future treatment
#'
#'
#' @export
#'
future.loss <- function(z.next, t.next, zp, N, design, int, lossfunc, beta, y, bayes, dyn=NULL, ...){

  if (!is.null(int)) {
    design.next <- rbind(design, c(1, z.next, t.next, z.next*t.next))
  }else{
    design.next <- rbind(design, c(1, z.next, t.next))
  }


  #if horizon is 1, calculate loss. else iterate.
  if (N==1){

    loss <- lossfunc(Imat.beta(design.next, beta), ...)

  } else{
    if (!is.null(dyn)){

      dyn <- dyn+1
    }

    loss <- exp.loss( z.now=z.next, t.now=t.next, zp, N-1, design, int, lossfunc, beta, y, bayes,dyn,  ...)



  }

  return(loss)

}






#' Assuming a current response, break down the expected future optimality by cases for every combination of:
#' 1) future possible covariate
#' 2) future possible treatment
#' Find a weighted average across all these cases

#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param zp vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param int set to NULL if there are no interactions, set to T of there are interactions
#' @param lossfunc  a function for the optimality criterion to minimize
#' @param beta estimate of the regression coefficients
#' @param y responses that have been observed up until the current unit
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#'
#'
#' @return expected optimality assuming a current response
#'
#'
#' @export
#'

future.y <- function(y.now, z.now, t.now, zp, N, design, int, lossfunc, beta, y, bayes, dyn=NULL, ...){

  if(!is.null(int)){
    d.now <- c(1, z.now, t.now, z.now*t.now)
  }else{
    d.now <- c(1, z.now, t.now)
  }

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

  if(!is.null(int)){
    j  <- (ncol(design)-2)/2
  }else{
    j  <- ncol(design)-2
  }
  all.z <- expand.grid(rep(list(c(-1,1)),j))  #grid for all possible covariates
  names(all.z) <- NULL

  #for each possible covariate value, calculate loss when t.next=1
  loss.p <- apply(all.z, 1, future.loss , t.next=1, zp, N, design, int, lossfunc, sim.beta, sim.y, bayes, dyn, ...)#loss for each covariate combination when the future patient has tmt=1

  #for each possible covariate value, calculate loss when t.next=1
  loss.m <- apply(all.z, 1, future.loss , t.next=-1, zp, N, design, int, lossfunc, sim.beta, sim.y, bayes, dyn, ...) #loss for each covariate combination when the future patient has tmt=-1

  #find loss for optimal treatment (t*) for each possible covariate value
  loss <- ifelse(loss.p < loss.m, loss.p, loss.m)

  #expected loss: weighed according to probabilities of each covariate value
  if (!is.null(dyn)){
    exploss <- sum(loss*zp[dyn,])
  }else{
    exploss <- sum(loss*zp)

  }


  return(exploss)

}





#' Break down expected future optimality into two components:
#' 1) assuming that the current response is 0
#' 2) assuming that the current response is 1
#' Find the weighted average of the two cases
#'
#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param zp vector of probabilities for each level of covariate z (needs to in the same order as all.z below)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param int set to NULL if there are no interactions, set to T of there are interactions
#' @param lossfunc  a function for the optimality criterion to minimize
#' @param beta estimate of the regression coefficients
#' @param y responses that have been observed up until the current unit
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#'
#' @return expected optimality one step ahead in the future
#'
#'
#' @export
#'
exp.loss <- function(z.now, t.now, zp, N, design, int, lossfunc, beta, y, bayes, dyn=NULL, ...){

  #probability that y.now=1
  if(!is.null(int)){
    pi <- probi(c(1, z.now, t.now, z.now*t.now), beta)
  }else{
    pi <- probi(c(1, z.now, t.now), beta)
  }

  #expected loss, given y.now=1
  loss1 <- future.y(y.now=1, z.now, t.now, zp, N, design, int, lossfunc, beta, y, bayes, dyn, ...)

  #expected loss, given y.now=1
  loss0 <- future.y(y.now=0, z.now, t.now, zp, N, design, int, lossfunc, beta, y,bayes, dyn, ...)

  #weighted by probability that P(y=1)
  exp.loss<- pi*loss1 + (1-pi)*loss0


  return(exp.loss)
}








#' Allocate treatments according to an information matrix based optimality criterion allowing for a non-myopic approach.
#' We assume a logistic model for the response and simulate responses sequentially.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param zprobs probabilities for each covariate value being 1
#' @param N natural number greater than 0 for horizon
#' @param lossfunc a function for the optimality criterion to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#' @return Design matrix D, all estimates of beta, final estimate of beta, responses y
#'
#'
#' @export
logit.nonmy <- function(covar,  true.beta, init, z.probs, N, int=NULL, lossfunc=calc.y.D, same.start=NULL, rand.start=NULL, stoc=T,
                        bayes=T, u=NULL,  true.bvcov=NULL, dyn=NULL, ...){

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe

  opt <- c()
  # randomly select treatment for first unit

  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }

  # starting design


  if (!is.null(same.start)) {
    design <-same.start
  }else if (!is.null(rand.start)) {
    design <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else if (!is.null(int)){
    design<-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
  } else{
    design <-    logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

  }

  rownames(design)<- NULL

  pi <- apply(design, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }


  #find initial estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
      beta <- coef(bayesglm (y~design[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm (y~design[,-1], family=binomial(link="logit")))
  }
  all.beta <- beta


  if (is.numeric(z.probs) | is.data.frame(z.probs)){
    zp <- zpfunc(z.probs)
  }


  for (i in (init+1):n){

    if (z.probs[1]=="learn"){
      zpl <- learn.zprobs(design, expand.grid(rep(list(c(-1,1)),j)) , j)
      if(!is.null(dyn)){


        Nprobs <- as.data.frame(matrix((rep(zpl, N)),  nrow=N))

        dyn=1
      }else{

        Nprobs <- zpl
      }

    }else{
      if(!is.null(dyn)){
      Nprobs <-data.frame(zp[(i+1):(i+N),])

      dyn=1
      }else{
        Nprobs <- zp
      }
    }


    #allocate treatment which minimizes expected loss
    eloss.p <- exp.loss(z.now=as.numeric(covar[i,]), t.now=1, zp=Nprobs, N, design, int, lossfunc, beta, y, bayes, dyn, ...) #expected loss for treatment 1
    eloss.m <- exp.loss(z.now=as.numeric(covar[i,]), t.now=-1,  zp=Nprobs, N, design, int , lossfunc, beta, y, bayes, dyn, ...)#expected loss for treatment -1

    probs <- (1/eloss.m)/(1/eloss.p+1/eloss.m)

    if(stoc==T){
      new.tmt <- sample(c(-1,1),  1, prob=c(probs, 1-probs))         #Assign treatments
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
    if(!is.null(int)){
      new.d <- as.numeric(cbind(1, covar[i,], new.tmt, covar[i,]*new.tmt ))
    }else{
          new.d <- as.numeric(cbind(1,covar[i,], new.tmt))

    }


    design <- as.matrix(rbind(design, as.numeric(new.d)))     #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }

   y <- c(y, new.y)
    if(bayes==T){
      beta <- coef(bayesglm (y~design[,-1], family=binomial(link="logit")))  #update beta
    }else{
      beta <- coef(glm (y~design[,-1], family=binomial(link="logit")))  #update beta
    }
    all.beta <- rbind(all.beta, beta)                          #Store all betas

  if (!is.null(true.bvcov)){
    opt <- c(opt, lossfunc(Imat.beta(design, true.beta), ...))
  }else{
    opt <- c(opt, lossfunc(Imat.beta(design, beta), ...))
  }
  }


  design <- data.frame(design)
  rownames(design) <- NULL
  results <- list(D=design, y=y, betas=all.beta, beta = beta, opt=opt)

  return(results)

}




