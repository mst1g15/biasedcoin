

#' Calculate L-optimal criterion assuming a logistic model (with or without weights)
#' @param D design matrix
#' @param beta current estimate of coefficients
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#'
#' @return L-optimal criterion
#'
#' @export
calc.logit.wL <- function(Imat, cr, prior.scale=100, wr=NULL){
  r <- nrow(cr)
  var <- rep(0, r)


  for (j in 1:r){
    if(is.null(wr)){
      var[j] <- t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]
    }else{
      var[j] <-  as.matrix(wr[j])%*%t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]

    }
  }
  opt <- sum(var)

  return(opt)
}







#' Calculate L-optimal criterion assuming a logistic model (with or without weights) assuming continuous treatment,
#' given one new patient
#' @param t new treatment
#' @param D current deisgn matrix
#' @param z covariate values of new unit
#' @param beta current estimate of coefficients
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#'
#' @return L-optimal criterion
#'
#' @export
wLopt.t <- function(t, D, z, beta, cr, prior.scale=100, wr=NULL){
  if (is.null(dim(z))){
      new.d <-c(1, z, t, z*t)

  }else{
    new.d <-cbind(1, z, t, t*z)

  }

  #rownames(new.d) <- rownames(D)
  if (is.data.frame(new.d)){
     colnames(new.d) <- colnames(D)

  }
  X <- as.matrix(rbind(D, new.d))
  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  Imat <- t(X) %*% W %*% X


  r <- nrow(cr)
  var <- rep(0, r)


  for (j in 1:r){
    if(is.null(wr)){
      var[j] <- t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]
    }else{
      var[j] <-  as.matrix(wr[j])%*%t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]

    }
  }
  opt <- sum(var)

  return(opt)
}




#' Calculate L-optimal criterion assuming a logistic model (with or without weights) assuming continuous treatment,
#' given a vector of treatment values and a vector of covariates values for additional patients
#' @param t new treatment
#' @param z covariate values of new unit
#' @param D design matrix
#' @param beta current estimate of coefficients
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#'
#' @return L-optimal criterion
#'
#' @export
wLopt.t.init <- function(t, z, beta, cr, prior.scale=100, wr=NULL){

  X <-cbind(1, z, t, z*t)
  X <- as.matrix(X)
  #rownames(new.d) <- rownames(D)
  #colnames(new.d) <- colnames(D)
  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  Imat <- t(X) %*% W %*% X


  r <- nrow(cr)
  var <- rep(0, r)


  for (j in 1:r){
    if(is.null(wr)){
      var[j] <- t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]
    }else{
      var[j] <-  as.matrix(wr[j])%*%t(cr[j,])%*% solve(Imat + diag(1/prior.scale ,ncol(Imat))) %*% cr[j,]

    }
  }
  opt <- sum(var)

  return(opt)
}



#' Calculate DA-optimal criterion assuming a logistic model (with or without weights)
#' @param D design matrix
#' @param beta current estimate of coefficients
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if weights equal
#'
#' @return DA-optimal criterion
#'
#' @export
calc.logit.wDA <- function(Imat, cr, prior.scale=100, wr=NULL){
  r <- nrow(cr)



  if(is.null(wr)){
    DA <- det( (cr)%*% solve(Imat+ diag(1/prior.scale ,ncol(D))) %*% t(cr))
  }else{

    var <-rep(0, r)

    inv <- solve(Imat + diag(1/prior.scale, ncol(D)))
    for (m in 1:r){
      var[m] <-   as.matrix( t(cr[m,]) %*% inv %*% as.matrix((cr[m,]))) ^ wr[m]
    }
    DA <- prod(var)

  }

  return(DA)

}



#' Assuming a logistic model for the response, allocate treatments using a coordinate exchange algorithm according
#' to an information matrix-based optimality criterion.
#' Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param beta an estimate for the true values of beta
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param cr.lossfunc loss function appropriate for optimality criteria which are for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise

#' @return design matrix D
#'
#'
#' @export
#'


cr.logit.coord <- function(covar, beta, threshold, kappa,cr.lossfunc=calc.logit.wL,  k, wt=NULL, ... ){
  
  covar <- as.data.frame(covar)
  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(beta)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  
  if(is.null(wt)){
    wr <- rep(1, nrow(cr))
  }else{
    wr <- wt
  }
  Djs <- list() #list to store m potential designs
  Ms <- list() #list to store the information matrices of the m potential designs
  for (q in 1:k){
    
    
    Dj <- as.matrix(cbind(rep(1, n), covar=covar, tmt=sample(c(0,1), n, replace=T))) #design matrix with random treatment assignment
    Dj <- cbind(Dj, Dj[,2:(j+1)]*Dj[,(j+2)])   #append interaction columns
    
    
    
    
    repeat{
      D <- Dj
      
      for (i in 1:n){
        D[i, "tmt"] <- 1
        D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]  #add interactino column
        plusD <- cr.lossfunc(Imat.beta(D, beta), cr, prior.scale, wr)  # ... needed in case lossfunc is DA
       
        
        
        D[i, "tmt"] <- 0
        D[i, (j+3):(2*j+2)] <-  D[i, "tmt"]*D[i, 2:(j+1)]
        
        minusD <- cr.lossfunc(Imat.beta(D, beta), cr, prior.scale, wr)  
                

        if (plusD < minusD){
          D[i,"tmt"] <- 1
          D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]
          
        } else {
          
          D[i,"tmt"] <- 0
          
          D[i, (j+3):(2*j+2)] <- D[i,"tmt"]* D[i, 2:(j+1)]
          
        }
      }
      
      if (identical( cr.lossfunc(Imat.beta(D, beta), cr, prior.scale, wr)  , cr.lossfunc(Imat.beta(Dj, beta), cr, prior.scale, wr)  )) break #if no change to resulting design matrix, terminate
      else Dj <- D
      
    }
    
    Djs[[q]] <- D
    
    Ms[[q]]  <- Imat.beta(D, beta)
    
  }
  
  mindet <- which.min(unlist(lapply(Ms, cr.lossfunc, cr, prior.scale, wr)))
  
  D <- Djs[mindet][[1]]
  
  
  return(D)
  
}


#' Assuming a logistic model for the response, allocate treatment sequentially based on an optimality criterion for
#' linear combinations of parameters. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param cr.lossfunc loss function appropriate for optimality criteria which are for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#'
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final estimate of beta, all weights, all estimates of standard deviation,
#'  beta, probabilities for treatment assignment, all values of optimalities (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#'
#' @export


cr.logit.des <- function(covar, true.beta, threshold, kappa, init, cr.lossfunc, k, wt, int=T, prior.scale=100,
                    same.start=NULL, rand.start=NULL, stoc=T, bayes=T,  u=NULL, prior.default=T,
                    true.bvcov=NULL, ... ){

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  all.lopt <- all.dawopt <- c()
  r <- nrow(cr)
  loss.p <- loss.m <- c()
  beta <- rep(0, p)
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    tmt <-  sample(c(0,1), init, replace=T)
    D <- cbind(rep(1, init), covar[1:init,], tmt, covar[1:init,]*tmt )
  }else{
    D <-  as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=T, code=0,lossfunc=cr.lossfunc, cr,  prior.scale))

  }


  pi <- apply(D, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)  #generate first observation based on true beta
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }



    #find new estimate of beta by using logistic regression on the first init responses
    if (bayes==TRUE){
      if(prior.default==F){
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale ))
      }else{
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
      }
    }else {
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }

  all.beta <- beta

  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  true.pr <- pnorm(threshold, mean= cr %*% true.beta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  true.wr <- t(ifelse( true.pr >= kappa, true.pr, 0))
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr


  yprop <- sum(y==1)/init
  tmtprop <- sum(D[,"tmt"]==1)/init

  if (!is.null(true.bvcov)){
    all.lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr)
   # all.dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, true.wr)

  }else{
    all.lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
    #all.dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

  }




  for (i in (init+1):n){

        #design matrix if tmt 1 is assigned - split into cases including/excluding interactions



    if (wt!=T){
      wr=NULL
    }

    if (!is.null(true.bvcov)){
      beta <- true.beta
    }

    mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))


    d.plus <-  cr.lossfunc(Imat.beta(mat.plus, beta), cr, prior.scale, wr)

    mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0, j) ))

    d.minus <-  cr.lossfunc(Imat.beta(mat.minus, beta), cr, prior.scale, wr)

    loss.p <- c(loss.p, d.plus)
    loss.m <- c(loss.m, d.minus)


    probs <- (1/d.plus)/sum(1/d.plus+1/d.minus)

    if (is.na(probs)){                                                #this is needed in case of separation occurring
      probs <- 0
    }


    if (stoc==TRUE){
      new.tmt <- sample(c(0,1), 1, prob=c(1-probs, probs))         #Assign treatments
    }else{
      if (probs > 0.5) {
        new.tmt <- 1
      } else {
        new.tmt <- 0
      }


    }
    #new row of design matrix

    new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt))

    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL
    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)

    }else{
      new.y <- rbinom(1, 1, pi)
    }

    y <- c(y, new.y) #Simulate new observation



      if(bayes==T){
        if(prior.default==F){
          beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale))
        }else{
          beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
        }
      }else{
        beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))                #Update beta
      }

    all.beta <- rbind(all.beta, beta)                          #Store all betas
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    true.pr <- pnorm(threshold, mean= cr %*% true.beta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    true.wr <- t(ifelse( true.pr > kappa, true.pr, 0))
    wr <- t(ifelse( pr >= kappa, pr, 0))

    all.wr <- rbind(all.wr, wr)



    if (!is.null(true.bvcov)){
      lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, true.wr)
      #dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, true.wr)

    }else{
      lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
      #dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

    }

    all.lopt <- c(all.lopt,lopt)
   # all.dawopt <- c(all.dawopt, dawopt)



    yprop <- c( yprop, (sum(y==1)/i))
    tmtprop <- c( tmtprop, (sum(D[,"tmt"]==1)/i))
  }



  D <- data.frame(D)

  results <- list(D=D, y=y, all.beta=all.beta, all.wr= all.wr, beta = beta,
                  yprop=yprop, tmtprop=tmtprop, all.lopt=all.lopt,
                  #all.dawopt= all.dawopt,
                  loss.p=loss.p, loss.m=loss.m)

  return(results)

}








#' Assuming a logistic model for the response, allocate treatment at random. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param cr.lossfunc loss function appropriate for optimality criteria which are for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'

#'
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final estimate of beta, all weights, all estimates of standard deviation,
#'  beta, probabilities for treatment assignment, all values of optimalities (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#'
#' @export


cr.rand <- function(covar, true.beta, threshold, kappa, init, cr.lossfunc, k, wt, int=T, prior.scale=100,
                         same.start=NULL, rand.start=NULL,  bayes=T,  u=NULL,prior.default=T,
                    true.bvcov=NULL, ... ){

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  all.lopt <- all.dawopt <- c()
  r <- nrow(cr)
  beta <- rep(0, p)
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    tmt <-  sample(c(0,1), init, replace=T)
    D <- cbind(rep(1, init), covar[1:init,], tmt, covar[1:init,]*tmt )
  }else{
    D <-  as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=T, code=0,lossfunc=cr.lossfunc, cr,  prior.scale))

  }


  pi <- apply(D, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)  #generate first observation based on true beta
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }



  #find new estimate of beta by using logistic regression on the first init responses
  if (bayes==TRUE){
    if(prior.default==F){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale ))
    }else{
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }
  }else {
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.beta <- beta

  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr


  yprop <- sum(y==1)/init
  tmtprop <- sum(D[,"tmt"]==1)/init


  if (!is.null(true.bvcov)){
    lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, wr)
    dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, wr)

  }else{
    lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
    dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

  }

  all.lopt <- c(all.lopt,lopt)
  all.dawopt <- c(all.dawopt, dawopt)


  for (i in (init+1):n){



     new.tmt <- sample(c(0,1), 1)


    new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt))

    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL
    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)

    }else{
      new.y <- rbinom(1, 1, pi)
    }

    y <- c(y, new.y) #Simulate new observation



    if(bayes==T){
      if(prior.default==F){
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale))
      }else{
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
      }
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))                #Update beta
    }

    all.beta <- rbind(all.beta, beta)                          #Store all betas
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    wr <- t(ifelse( pr > kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)


    if (!is.null(true.bvcov)){
     lopt <- calc.logit.wL(Imat.beta(D, true.beta), cr, prior.scale, wr)
     dawopt <- calc.logit.wDA(Imat.beta(D,true.beta),  cr, prior.scale, wr)

    }else{
      lopt <- calc.logit.wL(Imat.beta(D, beta), cr, prior.scale, wr)
      dawopt <- calc.logit.wDA(Imat.beta(D,beta),  cr, prior.scale, wr)

    }

    all.lopt <- c(all.lopt,lopt)
    all.dawopt <- c(all.dawopt, dawopt)


    yprop <- c( yprop, (sum(y==1)/i))
    tmtprop <- c( tmtprop, (sum(D[,"tmt"]==1)/i))
  }



  D <- data.frame(D)

  results <- list(D=D, y=y, all.beta=all.beta, all.wr= all.wr, beta = beta,
                  yprop=yprop, tmtprop=tmtprop, all.lopt=all.lopt,all.dawopt= all.dawopt)

  return(results)

}






#' Assuming a logistic model for the response, allocate treatment sequentially based on an optimality criterion for
#' linear combinations of parameters. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param cr.lossfunc loss function appropriate for optimality criteria which are for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#'
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final estimate of beta, all weights, all estimates of standard deviation,
#'  beta, probabilities for treatment assignment, all values of optimalities (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#'
#' @export


cr.logit.cont <- function(covar, true.beta, threshold, kappa, init,  k, wt,  prior.scale=100,
                          same.start=NULL, rand.start=NULL,  bayes=T,  u=NULL, prior.default=T,
                          true.bvcov=NULL, ... ){

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  all.lopt <-
    r <- nrow(cr)

  beta <- rep(0, p)
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    tmt <-  sample(c(0,1), init, replace=T)
    D <- cbind(rep(1, init), covar[1:init,], tmt, covar[1:init,]*tmt )
  }else{
    D <-  as.matrix(logit.coord(as.data.frame(covar[1:init,]), beta,  k, int=T, code=0,lossfunc=cr.lossfunc, cr,  prior.scale))

  }


  pi <- apply(D, 1, probi, t(true.beta))
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)  #generate first observation based on true beta
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }



  #find new estimate of beta by using logistic regression on the first init responses
  if (bayes==TRUE){
    if(prior.default==F){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale ))
    }else{
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }
  }else {
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.beta <- beta

  crbeta <-  cr %*% beta
  pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr


  yprop <- sum(y==1)/init
  tmtprop <- sum(D[,"tmt"]==1)/init




  for (i in (init+1):n){

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions



    if (wt!=T){
      wr=NULL
    }



    findopt <- optim(runif(1, min=0, max=1), fn=wLopt.t,  D=D, z=as.numeric(covar[i, ]), beta =beta, cr=cr, prior.scale=prior.scale, wr=wr, method="L-BFGS-B",lower=0, upper=1)

    new.tmt <- findopt$par





    if (!is.null(true.bvcov)){
      new.opt <- wLopt.t(t=new.tmt, D=D, z=as.numeric(covar[i, ]), beta =t(true.beta), cr=cr, prior.scale=prior.scale, wr=wr)
    }else{
      new.opt <- findopt$value
    }


    all.lopt <- c(all.lopt, new.opt)


    #new row of design matrix

    new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt))

    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL
    pi <- probi(new.d, t(true.beta))       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)

    }else{
      new.y <- rbinom(1, 1, pi)
    }

    y <- c(y, new.y) #Simulate new observation



    if(bayes==T){
      if(prior.default==F){
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit"), prior.df=Inf, prior.scale=prior.scale))
      }else{
        beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
      }
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))                #Update beta
    }

    all.beta <- rbind(all.beta, beta)                          #Store all betas
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag( cr%*%solve( Imat.beta(D, beta) +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    wr <- t(ifelse( pr > kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)


    yprop <- c( yprop, (sum(y==1)/i))
    tmtprop <- c( tmtprop, (sum(D[,"tmt"]==1)/i))
  }



  D <- data.frame(D)

  results <- list(D=D, y=y, all.beta=all.beta, all.wr= all.wr, beta = beta,
                  yprop=yprop, tmtprop=tmtprop, all.lopt=all.lopt)

  return(results)

}








