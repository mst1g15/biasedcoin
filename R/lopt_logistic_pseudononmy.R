#' Allocate treatments according to an information matrix based optimality criterion allowing for a pseudo-nonmyopic approach.
#' We assume a logistic model for the response.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param cr.lossfunc a function for the optimality criterion to minimize which is appropriate for linear combinations of parameters
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param sim number of trajectories to simulate
#' @param z.probs vector of probabilities for each level of covariate z
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param coordex set to T if the coordinate exchange algorithm is used to assign treatments in the trajectory. Sequential approach used otherwsise.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#' @param ... further arguments to be passed to <lossfunc>

#'

#' @return Design matrix D, values of the objective function
#'
#'
#' @export

cr.simfuture.logis <- function(covar, true.beta, threshold, kappa, init, int, cr.lossfunc, k, wt=T,  sim, z.probs,
                            same.start=NULL, rand.start=NULL, stoc=T, bayes=T, u=NULL, prior.default=T, coordex=NULL,
                            true.bvcov=NULL, ...){

  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL
  design <- opt <- yprop.tot <- c()
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  all.loss.p <-  all.loss.m <- c()
  all.lopt <- all.dawopt <-c()
  r <- nrow(cr)
  #initial guess for beta
  if (is.null(int)){
    beta <- rep(0, j+2)
  }else{
    beta <- rep(0, 2*j+2)
  }

  if (!is.null(same.start)) {
    D <-same.start
  }else {
    #starting design
    if (!is.null(same.start)) {
      D <-same.start
    }else if (!is.null(rand.start)) {
      D <- cbind(rep(1, init), covar[1:init,], sample(c(0,1), init, replace=T))
    }else if (!is.null(int)) {
      D <- as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=T, code=0, lossfunc=cr.lossfunc, cr,  prior.scale))
    }else{
      D <- as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=NULL, code=0, lossfunc=cr.lossfunc, cr,  prior.scale))

    }

  }


  pi <- apply( D , 1, probi, true.beta)

  if (!is.null(u)){
      y <- ifelse(u[1:init] < pi, 1, 0)  #generate first observation based on true beta

  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }


  #find new estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
    if(prior.default==FALSE){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale ))
    }else{
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.betas <- beta

  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  wr <- t(ifelse( pr >= kappa, pr, 0))
  true.pr <- pnorm(threshold, mean= cr %*% true.beta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, true.beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  true.wr <- t(ifelse( true.pr > kappa, true.pr, 0))

  all.wr <- wr


  yprop <- sum(y==1)/init
  tmtprop <- sum(D[,"tmt"]==1)/init

  allt.lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr)
  not.lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)




  if (!is.null(true.bvcov)){
    all.lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr)
    all.dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, true.wr)

  }else{
    all.lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
    all.dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

  }
  for (i in (init+1):n){
    #i=11

    if (!is.null(int)){
      design.p <- rbind(D, c(1, as.numeric(covar[i,]), 1, as.numeric(covar[i,])))
      design.m <- rbind(D, c(1, as.numeric(covar[i,]), 0,  rep(0,j)))

    }else{
      design.p <- rbind(D, c(1, as.numeric(covar[i,]), 1))
      design.m <- rbind(D, c(1, as.numeric(covar[i,]), 0))
    }


    if(!is.null(ncol(z.probs))){
      z.probs.i <- as.data.frame(z.probs[ (i+1):n,])
    }else{
      z.probs.i <- z.probs
    }


    if (!is.null(true.bvcov)){
      beta <- true.beta
    }

    if (i!=n){

     if(!is.null(coordex)){

        loss.p <- future.coordex.logis(design.p, n-i, n,  z.probs.i,   beta, k, int, sim, code=0, lossfunc=cr.lossfunc,  cr, prior.scale, wr )
       loss.m <- future.coordex.logis(design.m, n-i, n, z.probs.i, beta, k, int,  sim, code=0, lossfunc=cr.lossfunc, cr, prior.scale, wr)

      }else{

        loss.p <- future.logis(design.p, n-i,  z.probs.i,  beta,  int, sim, code=0, lossfunc=cr.lossfunc,  cr, prior.scale, wr )
        loss.m <- future.logis(design.m, n-i,  z.probs.i, beta, int,  sim, code=0, lossfunc=cr.lossfunc, cr, prior.scale, wr)
      }

    }else {
      loss.p <- cr.lossfunc(Imat.beta(design.p, beta), cr, prior.scale, wr)
      loss.m <- cr.lossfunc(Imat.beta(design.m, beta), cr, prior.scale, wr)
    }



    all.loss.p <- c( all.loss.p, loss.p)
    all.loss.m <- c( all.loss.m, loss.m )


    probs <- (1/loss.p)/(1/loss.p+1/loss.m)



    if(stoc==T){
      new.tmt <- sample(c(1, 0), 1, prob=c(probs, 1-probs))         #Assign treatments
    }else{
      if (loss.p > loss.m) {
        new.tmt <- 0
      } else if (loss.p < loss.m) {
        new.tmt <- 1
      } else if (loss.p == loss.m) {
        new.tmt <- sample(c(0,1), 1)
      }
    }

    #new row of design matrix
    if (!is.null(int)){
      new.d <- c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt)
    } else{
      new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt))
    }


    D <- as.matrix(rbind(D, new.d))


    pi <- probi(D[i,], true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }
    y <- c( y, new.y)                               #Simulate new observation
    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.betas <- rbind(all.betas, beta)
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    wr <- t(ifelse( pr > kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)
    true.pr <- pnorm(threshold, mean= cr %*% true.beta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, true.beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    true.wr <- t(ifelse( true.pr > kappa, true.pr, 0))



    allt.lopt <- c( allt.lopt, calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr) )
    not.lopt <- c( not.lopt, calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr))




    if (!is.null(true.bvcov)){
      lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr)
      dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, true.wr)

    }else{
     lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
      dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

    }
    all.lopt <- c(all.lopt,lopt)
    all.dawopt <- c(all.dawopt, dawopt)


    yprop <- c( yprop, (sum(y==1)/i))
    tmtprop <- c( tmtprop, (sum(D[,"tmt"]==1)/i))

  }
  return(list(D=D,  y=y, all.betas=all.betas, all.wr= all.wr, beta = beta,
              yprop=yprop, tmtprop=tmtprop, all.lopt=all.lopt, all.dawopt=all.dawopt,
              allt.lopt = allt.lopt,  not.lopt=not.lopt ,
              loss.p= all.loss.p, loss.m= all.loss.m))


}


#' Allocate treatments according to weighted L-optimal objective function allowing for a pseudo-nonmyopic approach.
#' We assume a logistic model for the response and continuous treatment.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param sim number of trajectories to simulate
#' @param z.probs vector of probabilities for each level of covariate z
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#' @param ... further arguments to be passed to <lossfunc>

#'

#' @return Design matrix D, value of weighted L-optimal objective function
#'
#'
#' @export



cr.simfuture.logis.cont <- function(covar, true.beta, threshold, kappa, init,  sim, z.probs,k=2, wt, prior.scale=100,
                               same.start=NULL, rand.start=NULL,  bayes=T, u=NULL, prior.default=T,
                               true.bvcov=NULL, ...){

  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL
  opt <- yprop.tot <- c()
  allt.lopt <- betat.lopt <- not.lopt <- c()
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  r <- nrow(cr)
  #initial guess for beta
  if (is.null(int)){
    beta <- rep(0, j+2)
  }else{
    beta <- rep(0, 2*j+2)
  }

  if (!is.null(same.start)) {
    D <-same.start
  }else {
    #starting design
    if (!is.null(same.start)) {
      D <-same.start
    }else if (!is.null(rand.start)) {
      D <- cbind(rep(1, init), covar[1:init,], sample(c(0,1), init, replace=T))
    }else if (!is.null(int)) {
      D <- as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=T, code=0, lossfunc=cr.lossfunc, cr,  prior.scale))
    }else{
      D <- as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=NULL, code=0, lossfunc=cr.lossfunc, cr,  prior.scale))

    }

  }


  pi <- apply( D , 1, probi, t(true.beta))

  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)  #generate first observation based on true beta

  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }


  #find new estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
    if(prior.default==FALSE){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale ))
    }else{
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.betas <- beta

  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  true.pr <- pnorm(threshold, mean=  cr %*% (true.beta), sd=sqrt(diag( cr%*%solve( Imat.beta(D, t(true.beta)) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  wr <- t(ifelse( pr >= kappa, true.pr, 0))
  true.wr <- t(ifelse( true.pr >= kappa, true.pr, 0))
  all.wr <- wr


  yprop <- sum(y==1)/init
  tmtprop <- sum(D[,"tmt"]==1)/init



  for (i in (init+1):n){


    if(!is.null(ncol(z.probs))){
      z.probs.i <- as.data.frame(z.probs[ (i+1):n,])
    }else{
      z.probs.i <- z.probs
    }

    n.r <- n-i+1




    new.tmt <- optim(par=runif(1, 0, 1), wLopt.pseudo.t, method="Brent", lower=0, upper=1, D=D, z=as.numeric(covar[i,]),  beta=beta, cr=cr, prior.scale=prior.scale, wr=wr, M= sim, n.r=n.r, z.probs=z.probs.i)$par


    if (!is.null(true.bvcov)){
      opt <- c(opt, wLopt.t(new.tmt, D, as.numeric(covar[i, ]), t(true.beta) ,cr=cr, prior.scale=prior.scale, wr=true.wr))


    }else{
      opt <- c(opt, wLopt.t(new.tmt, D, as.numeric(covar[i, ]), t(beta) ,cr=cr, prior.scale=prior.scale, wr=wr))

    }



    allt.lopt <- c( allt.lopt, wLopt.t(new.tmt, D, as.numeric(covar[i, ]), t(true.beta) ,cr=cr, prior.scale=prior.scale, wr=true.wr) )
    not.lopt <- c( not.lopt, wLopt.t(new.tmt, D, as.numeric(covar[i, ]), t(beta) ,cr=cr, prior.scale=prior.scale, wr=wr))



      new.d <- c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt)


    D <- as.matrix(rbind(D, new.d))


    pi <- probi(new.d, t(true.beta))       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }
    y <- c( y, new.y)                               #Simulate new observation
    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.betas <- rbind(all.betas, beta)
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    wr <- t(ifelse( pr > kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)




    yprop <- c( yprop, (sum(y==1)/i))
    tmtprop <- c( tmtprop, (sum(D[,"tmt"]==1)/i))

  }
  return(list(D=D,  y=y, all.betas=all.betas, all.wr= all.wr, beta = beta,
              yprop=yprop, tmtprop=tmtprop, all.lopt=opt,
              allt.lopt= allt.lopt, not.lopt =not.lopt ))


}







#' Calculate average L-optimal criterion assuming a logistic model with continuous treatment after M possible trajectories
#' @param t new treatment
#' @param D current deisgn matrix
#' @param z covariate values of new unit
#' @param beta current estimate of regression parameters
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#' @param M number of trajectories
#' @param n.r length of each trajectory
#' @param z.probs assumed covariate distribution
#'
#'
#' @return average value of L-optimal objective function across M trajectories
#'
#' @export



wLopt.pseudo.t <- function(t, D, z,  beta, cr, prior.scale=100, wr=NULL, M, n.r, z.probs){

    new.d <-c(1, z, t, z*t)

  #colnames(new.d) <- colnames(D)
  D <- rbind(D, new.d)
  all.opt <- rep(0, M)

  for (i in 1:M){


    simcov <- gencov(z.probs, n.r, code=0)

    optdes <- optim(runif(n.r, min=0, max=1), wLopt.t, method="L-BFGS-B", lower=0, upper=1, D=D, z=simcov, beta =beta, cr=cr, prior.scale=prior.scale, wr=wr)

    all.opt[i] <- optdes$value

  }


  av.opt <- mean(all.opt)

  return(av.opt)

}





