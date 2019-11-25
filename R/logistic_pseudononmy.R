



#' Calculate expected optimality for a given trajectory using sequential approach.
#' We assume a logistic model for the response.
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param z.probs vector of probabilities for each level of covariate z
#' @param beta current estimate of regression coefficients
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param sim number of trajectories to simulate
#' @param code set to NULL if (-1,1) coding is used for the treatments. Set to 0 if (0, 1) is used.
#' @param lossfunc the objective function to minimize
#' @param ... further arguments to be passed to <lossfunc>

#' @return loss of the design matrix which includes trajectory
#'
#'
#' @export


future.logis <- function(D.fix, n.r, z.probs,  beta, int, sim, code=0, lossfunc,  ...){


  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)
      losses <- rep(0, sim)


  for (q in 1:sim){

    covar <- gencov(z.probs,n.r, code)


    D <- D.fix

    for (i in 1:n.r){

      if (!is.null(int)){
        mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))
      } else{
        mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
      }

      d.plus <-  lossfunc(Imat.beta(mat.plus,  beta), ...)                #Optimality criterion

      #design matrix if tmt -1 is assigned
      if (!is.null(int)){

        if (!is.null(code)){
          mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0, ncol(covar)) ))

        }else{
                  mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0 ))

        }
      } else{
        if (!is.null(code)){

        mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 1)) #Design matrix with new treatment =-1
        }else{
          mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0)) #Design matrix with new treatment =-1

        }
      }


      d.minus <- lossfunc(Imat.beta(mat.minus,  beta), ...)                 #Optimality criterion
      if (!is.null(code)){
        new.tmt <- ifelse(d.plus < d.minus, 1, 0)         #Assign treatments

      }else{
              new.tmt <- ifelse(d.plus < d.minus, 1, 0)         #Assign treatments

      }

      #new row of design matrix
      if (!is.null(int)){

        new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt))
      } else {
        new.d <- as.numeric(c(1, covar[i, ], new.tmt))
      }


      D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

      losses[q] <- lossfunc(Imat.beta(D, beta), ... )

    }


  }

  mean.loss <- mean(losses)


  return(mean.loss)



}






#' Calculate expected optimality for a given trajectory using the coordinate exchange approach
#' We assume a logistic model for the response.
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param n total number of patients in the experiment
#' @param z.probs vector of probabilities for each level of covariate z
#' @param beta current estimate of regression coefficients
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param sim number of trajectories to simulate
#' @param code set to NULL if (-1,1) coding is used for the treatments. Set to 0 if (0, 1) is used.
#' @param lossfunc the objective function to minimize
#' @param ... further arguments to be passed to <lossfunc>

#' @return loss of the design matrix which includes trajectory
#'
#' @export

future.coordex.logis <- function(D.fix, n.r, n, z.probs,  beta, k, int, sim, code=NULL, lossfunc,  ...){

  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)
    losses <- rep(0, sim)

  for (q in 1:sim){


    covar <- gencov(z.probs,n.r, code)

    Dms <- Ms <- list() #list to store k potential designs

    for (m in 1:k){
      if (!is.null(code)){
        D.r <- as.matrix(cbind(rep(1, n.r), covar, tmt=sample(c(0,1), n.r, replace=T))) #design matrix with random treatment assignment

      }else{
              D.r <-as.matrix(cbind(1, covar, sample(c(-1,1), n.r, replace=T)))

      }

      if (!is.null(int)){
        if (dim(D.r)[1]==1){
        D.r <- cbind(D.r, t(D.r[,-c(1, j+2)]*D.r[,(j+2)]))
      }else{
        D.r <- cbind(D.r, D.r[,-c(1, j+2)]*D.r[,(j+2)])

      }
      }

      colnames(D.r) <-  colnames(D.fix)

      Dm <- rbind(D.fix, D.r)


      #print(D.fix)
      #print(D.r)
            repeat{
        D <- Dm
        for (i in (n.fix+1):n){

          D.new <- D




          if (D[i, "tmt"] == 1){
            if (!is.null(code)){
              D.new[i, "tmt"] <- 0
              if (!is.null(int)){
                D.new[i, (j+3):D.c] <-rep(0,j)
              }

            }else{
              D.new[i, "tmt"] <- -1
              if (!is.null(int)){
                D.new[i, (j+3):D.c] <- -D.new[i, 2:(j+1)]
              }
            }


          }else{
            D.new[i, "tmt"] <- 1

            if (!is.null(int)){
              D.new[i, (j+3):D.c] <- D.new[i, 2:(j+1)]
            }


          }

          oldD.l <- lossfunc(Imat.beta(D, beta) ,  ... )

          newD.l <- lossfunc(Imat.beta(D.new, beta) ,  ... )

        #print(oldD.l)
        #print(newD.l)


          if (oldD.l <=  newD.l){
            D <- D
          } else {
            D <- D.new
          }



        }


        l1 <- lossfunc(Imat.beta(D, beta) , ...)
        l2 <- lossfunc(Imat.beta(Dm, beta) , ...)

       # cat(q,  l1, l2, "\n")

        if(l1==l2){
          break
        }else Dm <- D
        #if no change to resulting design matrix, terminate
      }
         Dms[[m]] <- D
         Ms[[m]]  <- Imat.beta(D, beta)
    }


    mindes <-  which.min(unlist(lapply(Ms, lossfunc, ...)))  #find the optimum design matrix

    result <-Dms[mindes][[1]]
    rownames(result) <- NULL
    result <- as.data.frame(result)
    losses[q] <- lossfunc(Imat.beta(result, beta),  ...  )


    }

  mean.loss <- mean(losses)

  return(mean.loss)

}





#' Allocate treatments according to an information matrix based optimality criterion allowing for a pseudo-nonmyopic approach.
#' We assume a logistic model for the response.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param sim number of trajectories to simulate
#' @param z.probs vector of probabilities for each level of covariate z
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc the objective function to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param coordex set to T if coordinate exchange algorithm is used to allocate treatments in the trajectory, set to NULL for sequential approach.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return Design matrix D, the value of value of objective function
#'
#'
#' @export

simfuture.logis <- function(covar, true.beta, init, k, sim, z.probs, int=NULL, lossfunc, same.start=NULL, rand.start=NULL, stoc=T, bayes=T,  coordex=NULL,u=NULL, true.bvcov=NULL, ...){

  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL
  design <- opt <- yprop.tot <- c()
  all.loss.p <- all.loss.m <- c()


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
      D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
    }else if (!is.null(int)) {
      D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
    }else{
      D <-  logit.coord(covar[1:init,], beta, 2, int=NULL, lossfunc, ...)

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
    beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.betas <- beta


  for (i in (init+1):n){
  #i=11

    if (!is.null(int)){
      design.p <- rbind(D, c(1, as.numeric(covar[i,]), 1, as.numeric(covar[i,])))
      design.m <- rbind(D, c(1, as.numeric(covar[i,]), -1, - as.numeric(covar[i,])))

    }else{
      design.p <- rbind(D, c(1, as.numeric(covar[i,]), 1))
      design.m <- rbind(D, c(1, as.numeric(covar[i,]), -1))
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

        loss.p <- future.coordex.logis(design.p, n-i, n,  z.probs.i,   beta, k, int, sim, code=NULL, lossfunc,  ... )
        loss.m <- future.coordex.logis(design.m, n-i, n, z.probs.i, beta, k, int,  sim, code=NULL, lossfunc, ...)

      }else{

      loss.p <- future.logis(design.p, n-i,  z.probs.i,  beta,  int, sim, code=NULL, lossfunc,  ... )
      loss.m <- future.logis(design.m, n-i,  z.probs.i, beta, int,  sim, code=NULL, lossfunc, ...)
      }

    }else {
      loss.p <- lossfunc(Imat.beta(design.p, beta), ...)
      loss.m <- lossfunc(Imat.beta(design.m, beta), ...)
    }


    all.loss.p <- c(all.loss.p, loss.p)
    all.loss.m <- c(all.loss.m, loss.m)


    probs <- (1/loss.p)/(1/loss.p+1/loss.m)



    if(stoc==T){
      new.tmt <- sample(c(-1,1), 1, prob=c(probs, 1-probs))         #Assign treatments
    }else{
      if (loss.p > loss.m) {
        new.tmt <- -1
      } else if (loss.p < loss.m) {
        new.tmt <- 1
      } else if (loss.p == loss.m) {
        new.tmt <- sample(c(-1,1), 1)
      }
    }

    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else{
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))


    pi <- probi(new.d, true.beta)       #Compute new pi

    if (!is.null(u)){
          new.y <- ifelse(u[i] < pi, 1, 0)

    }else{
      new.y <- rbinom(1, 1, pi)
    }




    y <- c( y, new.y)                               #Simulate new observation
    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    if(bayes==T){
      beta <- coef(bayesglm(y~D[,-1], family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
    }
    all.betas <- rbind(all.betas, beta)


    if (!is.null(true.bvcov)){
      opt <- c(opt, lossfunc(Imat.beta(D, true.beta), ...))
    }else{
      opt <- c(opt, lossfunc(Imat.beta(D, beta), ...))
    }

  }
  return(list(D=D, y=y, beta=beta, all.betas = all.betas, yprop.tot=yprop.tot, opt=opt,

              loss.p=all.loss.p, loss.m=all.loss.m))


}



#' Calculate expected optimality for a given trajectory using sequential approach when the treatment is continuous
#' We assume a logistic model for the response.
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param z.probs vector of probabilities for each level of covariate z
#' @param beta current estimate of regression coefficients
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param sim number of trajectories to simulate
#' @param lossfunc the objective function to minimize
#' @param ... further arguments to be passed to <lossfunc>

#' @return loss for the design matrix which includes trajectory
#'
#'
#' @export

future.logis.cont <- function(D.fix, n.r, z.probs,  beta, int, sim, lossfunc,  ...){


  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)
  losses <- rep(0, sim)


  for (q in 1:sim){

    covar <- gencov(z.probs,n.r)


    D <- D.fix

    for (i in 1:n.r){




      findopt <- optim(t, lossfunc,  method="Brent",lower=-1, upper=1, D=D, z=as.numeric(covar[i, ]), int=int, beta =beta)

      new.tmt <- findopt$par



      #new row of design matrix
      if (!is.null(int)){

        new.d <- as.numeric(c(1, as.numeric(covar[i, ]), new.tmt, as.numeric(covar[i, ])*new.tmt))
      } else {
        new.d <- as.numeric(c(1, covar[i, ], new.tmt))
      }


      D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

      losses[q] <- findopt$value

    }


  }

  mean.loss <- mean(losses)


  return(mean.loss)



}




#' Allocate continuous treatments according to an information matrix based optimality criterion allowing for a pseudo-nonmyopic approach.
#' We assume a logistic model for the response.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the data generating mechanism
#' @param init the number of units in the initial design
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param sim number of trajectories to simulate
#' @param z.probs vector of probabilities for each level of covariate z
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc the objective function to minimize
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param true.bvcov if set to T, use the true parameter values to compute obejctive function. If set to NULL, use estimated parameter values.
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return Design matrix D, value of the objective function
#'
#'
#' @export

simfuture.logis.cont <- function(covar, true.beta, init, k, sim, z.probs, int=NULL, lossfunc, same.start=NULL, rand.start=NULL, bayes=T,  u=NULL, true.bvcov=NULL, ...){

  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL
  opt <- yprop.tot <- c()


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
      D <- cbind(rep(1, init), covar[1:init,], runif(init, -1,1))
    }else if (!is.null(int)) {
      D <-  logit.coord(covar[1:init,], beta, 2, int=T, lossfunc, ...)
    }else{
      t <-optim(par = runif(init, -1,1), Dopt.y.t.init, method="L-BFGS-B", lower=-1, upper=1, z=as.numeric(covar[1:init,]), int=int, beta=c(rep(0, 3)), epsilon=0.00001)$par
      D <-cbind(1, covar[1:init,], t)
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
    beta <- coef(bayesglm(y~as.matrix(D[,-1]), family=binomial(link="logit")))
  }else{
    beta <- coef(glm(y~D[,-1], family=binomial(link="logit")))
  }

  all.betas <- beta

  for (i in (init+1):(n)){


    if (z.probs[1]=="learn"){
      z.probs.i <- learn.zprobs(D, expand.grid(rep(list(c(-1,1)),j)) , j)[2]
    } else if(!is.null(ncol(z.probs))){
      z.probs.i <- as.data.frame(z.probs[ (i+1):n,])
    }else{
      z.probs.i <- z.probs
    }

    n.r <- nrow(z.probs.i)
    if (is.null(n.r)==T){
      n.r <- 1
    }

    new.tmt <- as.matrix(optim(runif(1, -1, 1), Dopt.pseudo.t, method="Brent", lower=-1, upper=1, D=D, z=covar[i,], int=int, beta=beta, M=sim, n.r=n.r, z.probs=z.probs.i)$par)


    if (!is.null(true.bvcov)){
      opt <- c(opt, Dopt.y.t(t=new.tmt, D=D, z=covar[i,], int=int, beta=true.beta ))


    }else{
          opt <- c(opt, Dopt.y.t(t=new.tmt, D=D, z=covar[i,], int=int, beta=beta ))

    }

   z.now <- as.matrix(covar[i,])

    #new row of design matrix
    if (!is.null(int)){
      new.d <- cbind(1, z.now, new.tmt, new.tmt%*%z.now)
    } else{
      new.d <-cbind(1, z.now, new.tmt)
    }

    D <- rbind(D, as.numeric(new.d))

    pi <- probi(c(new.d), true.beta)       #Compute new pi

    if (!is.null(u)){
      new.y <- ifelse(u[i] < pi, 1, 0)

    }else{
      new.y <- rbinom(1, 1, pi)
    }




    y <- c( y, new.y)                               #Simulate new observation
    yprop.tot <- c(yprop.tot, sum(y==1)/i)
    if(bayes==T){
      beta <- coef(bayesglm(y~as.matrix(D[,-1]), family=binomial(link="logit")))
    }else{
      beta <- coef(glm(y~as.matrix(D[,-1]), family=binomial(link="logit")))
    }
    all.betas <- rbind(all.betas, beta)



  }
  return(list(D=D, y=y, beta=beta, all.betas = all.betas, yprop.tot=yprop.tot, opt=opt))


}






#' Calculate average D-optimal criterion assuming a logistic model with continuous treatment after M possible trajectories
#' @param t new treatment
#' @param D current deisgn matrix
#' @param z covariate values of new unit
#' @param int set to TRUE if treatment-covariate interactions are included in the model
#' @param beta current estimate of regression parameters
#' @param M number of trajectories
#' @param n.r length of each trajectory
#' @param z.probs assumed covariate distribution
#'
#'
#' @return average value of D-optimal objective function across M trajectories
#'
#' @export


Dopt.pseudo.t <- function(t, D, z, int=NULL, beta, M, n.r, z.probs){

  row.names(z) <-NULL
  if (!is.null(int)&length(z)==1){
    new.d <-c(1, z, t, z*t)
  }else if (!is.null(int)&length(z)>1){
    new.d <-cbind(1, z, t, t%*%as.matrix(z))
  }else {
    new.d <-as.matrix(cbind(1, z, t))
  }



  colnames(new.d) <- colnames(D)
  D <- rbind(D, new.d)
  all.dopt <- rep(0, M)

  for (i in 1:M){


    simcov <- gencov(z.probs, n.r, code=NULL)

    optdes <- optim(runif(n.r, min=-1, max=1), Dopt.y.t, method="L-BFGS-B", lower=-1, upper=1, D=D, z=simcov, int=int, beta =beta)

    all.dopt[i] <- optdes$value

  }


  av.dopt <- mean(all.dopt)

  return(av.dopt)

}


