


#' Given a design, the current estimate of beta (or true value of beta), compute the bias matrix
#'
#' @param beta the current estimates of the parameter values
#' @param D current design matrix
#' @param sim number of simulated betas to generate
#' @param true.bvcov set to the true values of beta if the bias matrix is to be computed using the true values
#' @return bias matrix
#'
#'
#' @export

calc.bias <- function(beta,  D, sim, true.bvcov=NULL){

  #to store the betas

  if (!is.null(true.bvcov)){
   b <- true.bvcov
  }else{
    b <- beta
  }
  simbetas <- matrix(0, ncol=length(beta), nrow=sim)

  for (i in 1:sim){

    #generate a new set of responses under the design D, assuming the current estimate of beta
    pi <- apply(D, 1, probi, b)
    y <- as.numeric(rbinom(nrow(D), 1, pi))

    #given this new set of responses, find a new estimate of beta
    beta.boot <- coef(bayesglm (y~D[,-1], family=binomial(link="logit")))


    simbetas[i,] <- beta.boot



  #average across the sim new estimates of beta
  meanb <- colMeans(simbetas)

  #compute bias matrix

  if (!is.null(true.bvcov)){
    bias.mat<- (true.bvcov-meanb)%*%t(true.bvcov-meanb)
  }else{
    bias.mat<- (beta- meanb)%*%t(beta - meanb)
  }

  }

  return(bias.mat)
}




#' Given a design, the current estimate of beta (or true value of beta), compute the var-covar matrix
#'
#' @param beta the current estimates of the parameter values
#' @param D current design matrix
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @param true.bvcov set to the true values of beta if the var-covar matrix is to be computed using the true values
#' @return var-covar matrix
#'
#'
#' @export

calc.vcov <- function( beta, D, epsilon=0.00001, true.bvcov=NULL){

  if (!is.null(true.bvcov)){
    pi <- apply(D, 1, probi, true.bvcov)

  }else{
    pi <- apply(D, 1, probi, beta)
  }
  w <- pi*(1-pi)

  vcov <- solve(t(D)%*%diag(w)%*%(D)+epsilon*diag(ncol(D)))

  return(vcov)

}




#' Given a design, the current estimate of beta (or true value of beta), compute the MSE matrix
#'
#' @param beta the current estimates of the parameter values
#' @param D current design matrix
#' @param y vector of responses
#' @param sim number of simulated betas to generate
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @param true.bvcov set to the true values of beta if the var-covar matrix is to be computed using the true values
#' @return MSE matrix
#'
#'
#' @export


calc.mse <- function(beta, D, y, sim, epsilon=0.00001, true.bvcov=NULL){

  pi <- probi(D[nrow(D),],  beta)
  new.y <- as.numeric(rbinom(1, 1, pi))
  y <- c(y, new.y)

  model <- bayesglm (y~D[,-1], family=binomial(link="logit"))
  sim.beta <- coef(model)
  varcov <- calc.vcov(beta, D, epsilon, true.bvcov)
  bias  <- calc.bias(sim.beta, D, sim, true.bvcov)
  mse  <- varcov + bias

  return(mse)
}


if(FALSE){

#taken from sim nony func logis.R with some more bits added to output
###### edit this one in the logit file


logit.des.bayes <- function(covar,  true.beta, init, int=NULL, sim, lossfunc=calc.y.D,  same.start=NULL, rand.start=NULL, epsilon=0.00001, true.bvcov=NULL,...){
  #purpose: allocate treatments according to an information matrix based optimality criterion
  #         when you have a logistic model for the response. We simulate responses sequentially
  #input:   dataframe covar of covariate values
  #         beta0 a vector of initial parameter guesses
  #         init the number of observations to include in the preliminary design
  #         int indicates that you want to include all tmt-cov interactions
  #         lossfunc is the optimality criterion
  #output:  design matrix D
  #         responses y
  #         matrix of beta estimates, beta

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  all.probs <-rep(0, n-init)
  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }
  opt <- yprop.tot <- yprop <- tmtprop <- c()
  all.a.mse <- all.a.bias <- all.d.mse <- all.a.var <- all.beta <- c()
  all.bias <- all.var <- list()

  #starting design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else if ((!is.null(int))) {
    D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
  }else{
    D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

  }




  pi <- apply(D, 1, probi, true.beta)

  y <- as.numeric(rbinom(init, 1, pi))   #generate first observation based on true beta

  #find new estimate of beta by using logistic regression on the first init responses
  model <- bayesglm (y~D[,-1], family=binomial(link="logit"))

  beta <- coef(model)
  varcov <- calc.vcov(beta, D, epsilon, true.bvcov)
  a.var <- sum(diag(varcov))
  bias  <- calc.bias(beta, D, sim, true.bvcov)
  a.bias <- sum(diag(bias))
  mse  <- varcov + bias
  a.mse <- sum(diag(mse))
  d.mse <- det(mse)

  all.bias <- c(all.bias, list(bias))
  all.var <- c(all.var, list(varcov))


  for (i in (init+1):n){

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions
    if (!is.null(int)){
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))
    } else{
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
    }


    d.plus <-  lossfunc(Imat.beta(mat.plus,  beta), ...)                #Optimality criterion

    #design matrix if tmt -1 is assigned
    if (!is.null(int)){
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1, -as.numeric(covar[i, ]) ))
    } else{
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1)) #Design matrix with new treatment =-1
    }


    d.minus <- lossfunc(Imat.beta(mat.minus,  beta), ...)                 #Optimality criterion

    #probs <- (1/d.minus)/(1/d.plus+1/d.minus)
   # all.probs[i-init] <- probs


    #new.tmt <- sample(c(-1,1), 1, prob=c(probs, 1-probs))         #Assign treatments
    if (d.plus < d.minus ){
      new.tmt <- 1
    }else{
      new.tmt <- -1
    }
    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else{
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi#

    y <- c( y, rbinom(1, 1, pi))                               #Simulate new observation

    model <- bayesglm (y~D[,-1], family=binomial(link="logit"))
    beta <- coef(model)
    varcov <- calc.vcov(beta, D, epsilon, true.bvcov)
    a.var <- sum(diag(varcov))
    bias  <- calc.bias(beta, D, sim, true.bvcov)
    a.bias <- sum(diag(bias))
    mse  <- varcov + bias
    a.mse <- sum(diag(mse))
    d.mse <- det(mse)                  #Store all betas

    all.bias <- c(all.bias, list(bias))
    all.var <- c(all.var, list(varcov))

    if (!is.null(true.bvcov)){

    opt <- c(opt, lossfunc(Imat.beta(D, true.bvcov), ...))
    }else{
      opt <- c(opt, lossfunc(Imat.beta(D, beta), ...))

    }
    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    #yprop <- rbind(yprop, table(D[,2],y)[,2] /table(D[,2]))
    #tmtprop <- rbind(tmtprop, table(D[,2],D[,3])[,2] /table(D[,2]))
    all.beta <- rbind(all.beta, beta)                          #Store all betas
    yprop <- c(yprop, sum(y==1)/i)

    all.a.var <- c(all.a.var, a.var)

    all.a.bias <- c(all.a.bias, a.bias)



    all.a.mse <- c(all.a.mse, a.mse)
    all.d.mse <- c(all.d.mse, d.mse)

  }


  row.names(D) <-NULL

  D <- as.data.frame(D)
  #finalprop <- table(y, D[,"covar"])[2,]/table( D[,"covar"])


  results <- list(D=D, y=y, betas=all.beta, beta = beta, all.probs = all.probs,
                  #yprop=yprop,
                  yprop.tot=yprop.tot,
                  opt=opt,  all.a.bias=all.a.bias,  all.a.var = all.a.var,
                   all.a.mse = all.a.mse, all.d.mse = all.d.mse,
                  all.bias=all.bias, all.var =all.var
                  #, finalprop=finalprop
  )

  return(results)

}


}


#' Allocate treatments according to the MSE matrix when a logisic model for the response is assumed.
#' We simulate responses sequentially.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc a function for the objective function to minimize
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @param true.bvcov set to the true values of beta if the mse matrix is to be computed using the true values
#' @param same.start set to the intial design if desired or set to NULL otherwise
#' @param ... further arguments to be passed to <logit.coord> and <lossfunc>
#' @return the design matrix D, responses y, all estimates of beta, final estimate of beta, probabilities of treatment
#' assignment, proportion of favorable responses, value of objective function, trace of var-covar matrix, trace of bias matrix,
#' trace of mse matrix, determinant of mse matrix
#'
#' @export

logit.mse <- function(covar,  true.beta, init, int=NULL, lossfunc=calc.y.D, epsilon=0.00001,  same.start=NULL,true.bvcov=NULL,...){
  #purpose:
  #input:   covar: dataframe of of covariate values
  #         true.beta: true parameter values
  #         init: the number of observations to include in the preliminary design
  #         int: indicates that you want to include all tmt-cov interactions
  #         lossfunc: objective function


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  all.probs <-rep(0, n-init)
  tmt <- yprop <- opt <- all.betas <- all.a.mse <- all.a.bias <- all.a.var <- all.d.mse <- c()
 all.bias <- all.var <- list()
  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }

  if (!is.null(same.start)){
    D <- same.start
  }else{


     if (!is.null(int)) {
      D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
    }else{
      D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

    }


  }

  pi <- apply(D, 1, probi, true.beta)
  y <- as.numeric(rbinom(init, 1, pi))   #generate first observation based on true beta

  #find new estimate of beta by using logistic regression on the first init responses
  model <- bayesglm (y~D[,-1], family=binomial(link="logit"))

  #
  beta <- coef(model)
  varcov <- calc.vcov(beta, D, epsilon, true.bvcov)
  a.var <-  sum(diag(varcov))
  bias  <- calc.bias(beta, D, sim, true.bvcov)
  a.bias <- sum(diag(bias))

  all.bias <- c(all.bias, list(bias))
  all.var <- c(all.var, list(varcov))

  mse  <- varcov + bias
  a.mse <- sum(diag(mse))
  d.mse <- det(mse)

  for (i in (init+1):n){

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions
    if (!is.null(int)){
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1, as.numeric(covar[i, ]) ))
    } else{
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
    }


    d.plus <-  sum(diag(calc.mse(beta, mat.plus,   y, sim)))                #Optimality criterion

    #design matrix if tmt -1 is assigned
    if (!is.null(int)){
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1, -as.numeric(covar[i, ]) ))
    } else{
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1)) #Design matrix with new treatment =-1
    }


    d.minus <-  sum(diag(calc.mse(beta, mat.minus,  y, sim)))                  #Optimality criterion

    #probs <- d.minus/(d.plus+d.minus)
    #all.probs[i-init] <- probs


    #new.tmt <- sample(c(-1,1), 1, prob=c(probs, 1-probs))         #Assign treatments

    if (d.plus < d.minus){
      new.tmt <- 1
    } else{
      new.tmt <- -1
    }

    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else {
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }
    tmt <- c(tmt, new.tmt)

    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    y <- c( y, rbinom(1, 1, pi))                               #Simulate new observation

    model <- bayesglm (y~D[,-1], family=binomial(link="logit"))
    beta <- coef(model)
    varcov <- calc.vcov(beta, D, epsilon, true.bvcov)
    a.var <- sum(diag(varcov))
    bias  <- calc.bias(beta, D, sim, true.bvcov)
    a.bias <- sum(diag(bias))
    mse  <- varcov + bias
    a.mse <- sum(diag(mse))
    d.mse <- det(mse)

    all.bias <- c(all.bias, list(bias))
    all.var <- c(all.var, list(varcov))

    all.betas <- rbind(all.betas, beta)                          #Store all betas
    if (!is.null(true.bvcov)){

      opt <- c(opt, lossfunc(Imat.beta(D, true.bvcov), ...))
    }else{
      opt <- c(opt, lossfunc(Imat.beta(D, beta), ...))

    }
    yprop <- c(yprop, sum(y==1)/i)

    all.a.var <- c(all.a.var, a.var)

    all.a.bias <- c(all.a.bias, a.bias)


    all.a.mse <- c(all.a.mse, a.mse)
    all.d.mse <- c(all.d.mse, d.mse)

  }


 row.names(D) <-NULL

  D <- data.frame(D)

  results <- list(D=D, y=y, betas=all.betas, beta = beta, all.probs = all.probs,
                  yprop=yprop, opt=opt, all.a.var = all.a.var, all.a.bias=all.a.bias,
                   all.a.mse = all.a.mse, all.d.mse = all.d.mse, all.bias=all.bias, all.var=all.var)

  return(results)

}

