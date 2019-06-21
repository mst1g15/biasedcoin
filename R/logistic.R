library("arm")
library("prodlim")

#' Compute P(y=1) for logistic regression
#'
#' @param x a row of the design matrix
#' @param beta coefficients of the model

#' @return Probability P(y=1)'
#'
#' @export

probi <- function(x, beta){

  x <- as.matrix(x)
  eta <- beta %*% x
  p <- (exp(eta)/(1+exp(eta)))  #logistic function

  return(p)

}


#' Compute the information matrix, given beta, for logistic regression
#'
#' @param X the design matrix
#' @param beta coefficients of the model

#' @return Information matrix
#'
#' @export

Imat.beta <- function(X, beta){

 X <- as.matrix(X)

  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  I.mat <- t(X) %*% W %*% X

  return(I.mat)

}



#' Compute the D-optimality criterion of the information matrix for logistic regression
#'
#' @param I the information matrix
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return D-optimality of the information matrix
#'
#' @export
#'
calc.y.D <- function(I,  epsilon=0.00001){

  k <- ncol(I)
  I <- as.matrix(I)
  M <- I + diag(k)*epsilon
  D <- 1/det(M)

  return(D)
}



#' Compute the A-optimality criterion of the information matrix for logistic regression
#'
#' @param I the information matrix
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return A-optimality of the information matrix
#'
#' @export
#'
calc.y.A <- function(I,  epsilon=0.00001){

  k <- ncol(I)
  I <- as.matrix(I)
  M <- solve(I + diag(k)*epsilon)
  D <- sum(diag(M))

  return(D)
}




#' Compute the D-optimality criterion of the information matrix for logistic regression with an added row of the design matrix
#' (to be used in conjunction with optim)
#'
#' @param t new treatment
#' @param D current deisgn matrix
#' @param z covariate values of new unit
#' @param int set to TRUE if treatment-covariate interactions are included in the model
#' @param beta current estimate of regression parameters
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return D-optimality criterion
#'
#' @export
#'
Dopt.y.t <- function(t, D, z, int, beta, epsilon=0.00001){
  k <- ncol(D)
  if (!is.null(int)&length(z)==1){
    new.d <-c(1, z, t, z*t)
  }else if (!is.null(int)&length(z)>1){
    new.d <-cbind(1, z, t, t%*%as.matrix(z))
  }else {
    new.d <-as.matrix(cbind(1, z, t))
  }

  row.names(new.d) <- rownames(D) <- NULL
  colnames(new.d) <- colnames(D)
  X <- as.matrix(rbind(D, new.d))
  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  I <- t(X) %*% W %*% X
  opt <- 1/det(I + diag(k)*epsilon)

  return(opt)
}



#' Compute the D-optimality criterion of of an initial design for logistic regression
#' (to be used in conjunction with optim)
#'
#' @param t vector of treatments
#' @param vector of covariate values
#' @param int set to TRUE if treatment-covariate interactions are included in the model
#' @param beta current estimate of regression parameters
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return D-optimality criterion
#'
#' @export
Dopt.y.t.init <- function(t, z, int, beta, epsilon=0.00001){

  if (!is.null(int)){
    X <-cbind(1, z, t, z*t)
  }else{
    X <-cbind(1, z, t)
  }
  X <- as.matrix(X)
  k <- ncol(X)
  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  I <- t(X) %*% W %*% X

  opt <- 1/det(I + diag(k)*epsilon)

  return(opt)
}




#' Compute the G-optimality criterion of the information matrix for logistic regression
#' @param I the information matrix
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return G-optimality of the information matrix
#'
#' @export
#'

calc.y.G <- function(I, epsilon=0.00001){

  k <- ncol(I)
  I <- as.matrix(I)
  X <- as.matrix(cbind(1, expand.grid(rep(list(c(-1,1)),k-1))))
  G <- max(diag((X)%*%solve(I + diag(k)*epsilon)%*%t(X)))

  return(G)

}


#' Compute the DA-optimality criterion of the information matrix for logistic regression
#' @param I the information matrix
#' @param A a matrix where each column indicates the linear combination of parameters of interest
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @return DA-optimality of the information matrix
#'
#' @export
calc.y.DA <- function(I, A=t(matrix(c(0, 0, 1))), epsilon=0.00001){

  j <- ncol(I)
  I <- as.matrix(I)
  loss <- sum(diag(A%*%solve(I + epsilon*diag(j))%*%t(A)))

  return(loss)
}








#' Assuming a logistic model for the response, allocate treatment at random.
#' Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc a function for the optimality criterion to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#' @return design matrix D, responses y, all estimates of betas, final beta, probabilities for treatment assignment,
#' proportion of Y=1, proportion tmt=1, optimality
#'
#'
#' @export
#'
logit.rand <- function(covar, true.beta, init, int=NULL, lossfunc=calc.y.D,  same.start=NULL, rand.start=NULL, bayes=T, u=NULL, true.bvcov=NULL, ...){


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  all.probs <-rep(0, n-init)
  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }
  opt <- yprop.tot <- yprop <- tmtprop <- c()

  #starting design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else if (!is.null(int)) {
    D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
  }else{
    D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

  }




  pi <- apply(D, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }

  #find new estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
    beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.beta <- beta
  for (i in (init+1):n){


      new.tmt <- sample(c(-1,1), 1)


    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else{
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }


    y <- c( y, new.y)                               #Simulate new observation

    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.beta <- rbind(all.beta, beta)                          #Store all betas

    if (!is.null(true.bvcov)){
      opt <- c(opt, lossfunc(Imat.beta(D, true.beta), ...))
    }else{
      opt <- c(opt, lossfunc(Imat.beta(D, beta), ...))
    }



    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    #yprop <- rbind(yprop, table(D[,2],y)[,2] /table(D[,2]))


  }


  row.names(D) <-NULL

  D <- data.frame(D)


  results <- list(D=D, y=y, all.beta=all.beta, beta = beta, all.probs = all.probs, #yprop=yprop,
                  yprop.tot=yprop.tot,
                  opt=opt
  )

  return(results)

}










#' Assuming a logistic model for the response, allocate treatment sequentially based on an information matrix-based optimality criterion.
#' Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc a function for the optimality criterion to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final beta, probabilities for treatment assignment,
#' proportion of Y=1, proportion tmt=1, optimality
#'
#'
#' @export
#'
logit.des <- function(covar, true.beta, init, int=NULL, lossfunc=calc.y.D,  same.start=NULL, rand.start=NULL, stoc=T, bayes=T, u=NULL, true.bvcov=NULL, ...){


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  all.probs <-rep(0, n-init)
  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }
  loss.p <- loss.m <- c()
  opt <- yprop.tot <- yprop <- tmtprop <- c()

  #starting design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else if (!is.null(int)) {
    D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
  }else{
    D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

  }




  pi <- apply(D, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)
  }else{
      y <- as.numeric(rbinom(init, 1, pi))
  }

  #find new estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.beta <- beta
  for (i in (init+1):n){

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions
    if (!is.null(int)){
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))
    } else{
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
    }




    #design matrix if tmt -1 is assigned
    if (!is.null(int)){
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1, -as.numeric(covar[i, ]) ))
    } else{
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1)) #Design matrix with new treatment =-1
    }

    if (!is.null(true.bvcov)){
      d.plus <-  lossfunc(Imat.beta(mat.plus,  true.beta), ...)
      d.minus <- lossfunc(Imat.beta(mat.minus,  true.beta), ...)

    }else{
      d.plus <-  lossfunc(Imat.beta(mat.plus,  beta), ...)
      d.minus <- lossfunc(Imat.beta(mat.minus,  beta), ...)

    }


    #cat(d.plus, d.minus, "\n")


    loss.p <- c(loss.p, d.plus)
    loss.m <- c(loss.m, d.minus)
    probs <- (1/d.minus)/(1/d.plus+1/d.minus)
    all.probs[i-init] <- probs

    if(stoc==T){
          new.tmt <- sample(c(-1,1), 1, prob=c(probs, 1-probs))         #Assign treatments
    }else{
      if (d.plus > d.minus) {
        new.tmt <- -1
      } else if (d.plus < d.minus) {
        new.tmt <- 1
      } else if (d.plus == d.minus) {
        new.tmt <- sample(c(-1,1), 1)
      }
    }

    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else{
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }


    y <- c( y, new.y)                               #Simulate new observation

    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.beta <- rbind(all.beta, beta)                          #Store all betas
    if (!is.null(true.bvcov)){
      opt <- c(opt, lossfunc(Imat.beta(D, true.beta), ...))
    }else{
      opt <- c(opt, lossfunc(Imat.beta(D, beta), ...))
    }

    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    #yprop <- rbind(yprop, table(D[,2],y)[,2] /table(D[,2]))


  }


  row.names(D) <-NULL

  D <- data.frame(D)


  results <- list(D=D, y=y, all.beta=all.beta, beta = beta, all.probs = all.probs, #yprop=yprop,
                  yprop.tot=yprop.tot,
                  opt=opt,
                  loss.p=loss.p, loss.m=loss.m
  )

  return(results)

}









#' Assuming a logistic model for the response, allocate a continuous treatment sequentially based on an information matrix-based optimality criterion.
#' Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc a function for the optimality criterion to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final beta, probabilities for treatment assignment,
#' proportion of Y=1, proportion tmt=1, optimality
#'
#'
#' @export
#'
logit.cont <- function(covar, true.beta, init, int=NULL, lossfunc=Dopt.y.t,  same.start=NULL, rand.start=NULL,  bayes=T, u=NULL, true.bvcov=NULL, ...){


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  all.probs <-rep(0, n-init)
  if(!is.null(int)){
    beta <- rep(0, j+2+j)
  }else{
    beta <- rep(0, j+2)
  }
  opt <- yprop.tot <- yprop <- tmtprop <- c()

  #starting design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], runif(init, -1, 1))
  }else if (!is.null(int)) {
    D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
  }else{
    D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

  }




  pi <- apply(D, 1, probi, true.beta)
  if (!is.null(u)){
    y <- ifelse(u[1:init] < pi, 1, 0)
  }else{
    y <- as.numeric(rbinom(init, 1, pi))
  }

  #find new estimate of beta by using logistic regression on the first init responses
  if(bayes==T){
    beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.beta <- beta
  for (i in (init+1):n){


    findopt <- optim(par=runif(1, min=-1, max=1), fn=Dopt.y.t,  D=D, z=(covar[i, ]), int=int, beta =beta, method="L-BFGS-B",lower=-1, upper=1)

    new.tmt <- findopt$par




    if (!is.null(true.bvcov)){
      new.opt <- Dopt.y.t(t=as.numeric(new.tmt), D=D, z=(covar[i, ]), int=int, beta =(true.beta))
    }else{
      new.opt <- findopt$value
    }




    opt <- c(opt, new.opt)

    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else{
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

    pi <- probi(new.d, true.beta)       #Compute new pi
    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)
    }else{
      new.y <-  rbinom(1, 1, pi)
    }


    y <- c( y, new.y)                               #Simulate new observation

    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.beta <- rbind(all.beta, beta)                          #Store all betas




    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    #yprop <- rbind(yprop, table(D[,2],y)[,2] /table(D[,2]))


  }


  row.names(D) <-NULL

  D <- data.frame(D)


  results <- list(D=D, y=y, all.beta=all.beta, beta = beta,  #yprop=yprop,
                  yprop.tot=yprop.tot,
                  opt=opt)


  return(results)

}










#' Assuming a logistic model for the response, allocate treatments using a coordinate exchange algorithm according
#' to an information matrix-based optimality criterion.
#' Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param beta an estimate for the true values of beta
#' @param k an integer for the number of "outer loops"
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param code set to NULL if (-1,1) coding is used for the treatments. Set to 0 if (0, 1) is used.
#' @param lossfunc a function for the optimality criterion to minimize

#' @return design matrix D
#'
#'
#' @export
#'


logit.coord <- function(covar, beta, k, int=NULL, code=NULL, lossfunc=calc.y.D, ... ){

  covar <- as.data.frame(covar)
  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe


  Djs <- list() #list to store m potential designs
  Ms <- list() #list to store the information matrices of the m potential designs
  for (q in 1:k){

    if (!is.null(code)){
      Dj <- as.matrix(cbind(rep(1, n), covar=covar, tmt=sample(c(0,1), n, replace=T))) #design matrix with random treatment assignment

    }else{
          Dj <- as.matrix(cbind(rep(1, n), covar=covar, tmt=sample(c(-1,1), n, replace=T))) #design matrix with random treatment assignment


    }


    if (!is.null(int)){
      Dj <- cbind(Dj, Dj[,2:(j+1)]*Dj[,(j+2)])   #append interaction columns
    }



    repeat{
      D <- Dj

      for (i in 1:n){
        D[i, "tmt"] <- 1
        if (!is.null(int)){
          D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]  #add interactino column
        }
        plusD <- lossfunc(Imat.beta(D, beta), ...)  # ... needed in case lossfunc is DA
        if (!is.null(code)){
          D[i, "tmt"] <- 0
        }else{
                  D[i, "tmt"] <- -1 #calculate criterion for unit i assigned to tmt -1

        }


        if (!is.null(int)){
          D[i, (j+3):(2*j+2)] <-  D[i, "tmt"]*D[i, 2:(j+1)]
        }
        minusD <- lossfunc(Imat.beta(D, beta), ...)
        if (plusD < minusD){
          D[i,"tmt"] <- 1
          if (!is.null(int)){
            D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]
          }
        } else {

          if(!is.null(code)){
            D[i,"tmt"] <- 0

          }else{
                      D[i,"tmt"] <- -1

          }

          if (!is.null(int)){
            D[i, (j+3):(2*j+2)] <- D[i,"tmt"]* D[i, 2:(j+1)]
          }
        }
      }

      if (identical(lossfunc(Imat.beta(D, beta), ...), lossfunc(Imat.beta(Dj, beta), ...))) break #if no change to resulting design matrix, terminate
      else Dj <- D

    }

    Djs[[q]] <- D

    Ms[[q]]  <- Imat.beta(D, beta)

  }

  mindet <- which.min(unlist(lapply(Ms, lossfunc, ...)))

  D <- Djs[mindet][[1]]


  return(D)

}



