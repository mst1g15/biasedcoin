#' Calculate (weighted) L-optimal expected optimality for a given trajectory using sequential approach.
#' We assume a linear model for the response.
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param sim number of trajectories to simulate
#' @param z.probs vector of probabilities for each level of covariate z
#' @param lossfunc a function for the optimality criterion to minimize

#' @return loss of the design matrix which includes trajectory
#'
#'
#' @export

cr.future <- function(D.fix, n.r, sim, z.probs,lossfunc,  cr, prior.scale, wr, ...){


  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)
  losses <- rep(0, sim)


  for (q in 1:sim){

    covar <- gencov(z.probs,n.r, code=0)


    D <- D.fix

    for (i in 1:n.r){

        mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))


      d.plus <-  cr.lossfunc(mat.plus, cr, prior.scale, wr)                #Optimality criterion


      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0,ncol(covar)) ))


      d.minus <- cr.lossfunc(mat.minus, cr, prior.scale, wr)                 #Optimality criterion

      new.tmt <- ifelse(d.plus < d.minus, 1, 0)         #Assign treatments

      #new row of design matrix

        new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))



      D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

      losses[q] <- cr.lossfunc(D, cr, prior.scale, wr)

    }


  }

  mean.loss <- mean(losses)


  return(mean.loss)



}









#' Allocate treatments according to a (weighted) L-optimal criterion allowing for a pseudo-nonmyopic approach.
#' We assume a linear model for the response.
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param true.sigma the true parameter value for the standard deviation
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param cr.lossfunc loss function appropriate for optimality criteria which are for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param coordex set to T if the coordinate exchange algorithm is used to assign treatments in the trajectory. Sequential approach used otherwsise.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#'
#'
#' @return design matrix D, responses y, all estimates of betas, final estimate of beta, all weights, all estimates of standard deviation,
#'  beta, probabilities for treatment assignment, all values of optimalities (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#'
#' @export
cr.des.pseudo <- function(covar, true.beta, true.sigma, threshold,  kappa, init,  cr.lossfunc, sim, z.probs, k, wt=T,  prior.scale=100,
                   same.start=NULL, rand.start=NULL, stoc=T, bayes=T, prior.default=T, coordex=NULL, u=NULL, ... ){

  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)

  r <- nrow(cr)


  # use L-optimal design with equal weights for initial design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else{
    D <-  as.matrix(coord.cr(as.data.frame(covar[1:init,]), k, cr, cr.lossfunc, prior.scale))

  }


  #record proportion of patients in new treatment
  if (j==1){
    prop <-c(table(D[,2:(2+j)])[3]/ (table(D[,2:(2+j)])[1] + table(D[,2:(2+j)])[3]),
             table(D[,2:(2+j)])[4]/ (table(D[,2:(2+j)])[2] + table(D[,2:(2+j)])[4]))
  } else if (dim ( table(D[,2:(2+j)]))[3]<2 | is.na(dim ( table(D[,2:(2+j)]))[3]<2 )) {
    prop <- rep(0, r)
  } else {  prop <- as.vector(table(D[,2:(2+j)])[, ,2])/as.vector(table(D[,2:(1+j)]))
  }
  #set up vector
  all.prop <- prop



  D <- as.matrix(D)

  #generate init observations

  if (!is.null(u)){
    y <- D%*%true.beta+ + u[1:init]   #generate first observation based on true beta

  }else{
    y <- rnorm( n=init, mean=D%*%true.beta, sd=true.sigma)
  }



  #find first estimate of beta
  if(bayes!=T){
    model <- glm(y~D[,-1], family=gaussian)
  } else if(prior.default==T){
    model <- bayesglm(y~D[,-1], family=gaussian)
  }else{
    model <- bayesglm(y~D[,-1], family=gaussian,prior.df=Inf, prior.scale=prior.scale)
  }
  beta <- model$coefficients

  emp.sd <- sqrt(sum((y-D[1:init,]%*%beta)^2)/(init-p))
  sd <- sigma(model)

  #append to matrix of all estimates
  all.beta <- beta
  all.sd <- sd
  all.emp.sd <- sd
  #calculate weights
  crbeta <- cr %*% beta
  sdest <- sqrt(diag(sd^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))
  pr <- pnorm(threshold, mean= crbeta, sdest)
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr

  #calculate alpha
  alpha <- cr %*% beta < threshold
  all.alpha <- t(alpha)

  #calculate power
  tvalue <- qt(0.05, init-p)
  power.emp <- (cr %*% beta - threshold) /(sqrt(diag(sd^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))) < tvalue
  all.power.emp <- as.numeric(power.emp)

  power.th <- pt((tvalue + (threshold -cr %*% (true.beta)) /(sqrt(diag(sd^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr))))), init-p)
  all.power.th <- as.numeric(power.th)



  #calculate weighted optimality of each hypothesis
  all.lopt <- calc.wL(D, cr, prior.scale, wr)
  all.dawopt <- calc.wDA(D, cr, prior.scale, wr)
  #all.daopt <- calc.wDA(D, (matrix(cr[1:(r-1),])), prior.scale, wr=NULL)
  #if full cr matrix is used, we have numerical issues due to matrix being full column rank




  #initializing done

  #iterate
  for (i in (init+1):n){


    design.p<- rbind(D, c(1, as.numeric(covar[i, ]), 1, as.numeric(covar[i, ]) ))


    if (wt!=T){
      wr=NULL
    }

      design.m <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0,j) ))



    if(!is.null(ncol(z.probs))){
      z.probs.i <- as.data.frame(z.probs[ (i+1):n,])
    }else{
      z.probs.i <- z.probs
    }



    if (i!=n){

      if(!is.null(coordex)){

        loss.p <- future.coordex(design.p, n-i, n, k, sim, int=T, z.probs.i, code=0, cr.lossfunc, cr, prior.scale, wr)
        loss.m <- future.coordex(design.m, n-i, n, k, sim, int=T, z.probs.i, code=0, cr.lossfunc, cr, prior.scale, wr)

      }else{

        loss.p <- cr.future(design.p, n-i, sim,  z.probs.i,   cr.lossfunc,  cr, prior.scale, wr)
        loss.m <- cr.future(design.m, n-i, sim,  z.probs.i,  cr.lossfunc,  cr, prior.scale, wr)
      }


    }else {
      loss.p <- cr.lossfunc(design.p,   cr, prior.scale, wr)
      loss.m <- cr.lossfunc(design.m,   cr, prior.scale, wr)
    }




    if (stoc==T){

      if (loss.p==loss.m){
        probs = 0.5
      } else {
        probs <- (1/loss.p)/(1/loss.p+1/loss.m)
      }

      new.tmt <- sample(c(0,1), 1, prob=c(1-probs, probs))         #Assign treatments

    }else{
      if(loss.p < loss.m){
        new.tmt <- 1
      }else if (loss.p > loss.m){
        new.tmt <- 0
      }else{
        new.tmt <- sample(c(0,1), 1)
      }
    }


    #new row of design matrix

      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL


    if (!is.null(u)){
      new.y <-  new.d%*%true.beta + u[i]
    }else{
      new.y <- rnorm(1, new.d%*%true.beta, true.sigma)                               #Simulate new observation

    }

    y <- c(y, new.y)


    if (bayes!=T){
      model <- glm(y~D[,-1], family=gaussian)
    }else if(prior.default==T){
      model <- bayesglm(y~D[,-1], family=gaussian)
    }else{
      model <- bayesglm(y~D[,-1], family=gaussian,prior.df=Inf, prior.scale=prior.scale)
    }

    #model parameters
    beta <- coef(model)
    sd <- sigma(model)
    all.beta <- rbind(all.beta, beta)                          #Store all betas
    all.sd <- c(all.sd, sd)

    emp.sd <- sqrt(sum((y-D[1:i,]%*%beta)^2)/(i-p))

    all.emp.sd <- c(all.emp.sd, emp.sd)

    #calculate weights
    crbeta <- cr %*% beta
    pr <- pnorm(threshold, mean= crbeta, sd=sqrt(diag(sd^2 * cr%*%solve( t(D) %*%D +diag(1/prior.scale, ncol(D)))%*% t(cr))))
    wr <- t(ifelse( pr > kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)

    #calculate alpha
    alpha <- cr %*% beta < threshold
    all.alpha <- rbind(all.alpha, t(alpha))

    tvalue <- qt(0.05, i-p)
    power.emp <- (cr %*% beta - threshold) /(sqrt(diag(sd^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))) < tvalue
    all.power.emp <- rbind(all.power.emp, as.numeric(power.emp))

    power.th <- pt(tvalue + (threshold -cr %*% true.beta) /(sqrt(diag(true.sigma^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))), i-p)
    all.power.th <- rbind(all.power.th, as.numeric(power.th))

    D <- as.data.frame(D)

    #proportion of patients with new treatment
    if (j==1){
      prop <-c(table(D[,2:(2+j)])[3]/ (table(D[,2:(2+j)])[1] + table(D[,2:(2+j)])[3]),
               table(D[,2:(2+j)])[4]/ (table(D[,2:(2+j)])[2] + table(D[,2:(2+j)])[4]))
    } else if (dim ( table(as.data.frame(D)[,2:(2+j)]))[3]<2 | is.na(dim( table(as.data.frame(D)[,2:(2+j)]))[3])) {
      prop <- rep(0, r)
    } else {  prop <- as.vector(table(as.data.frame(D)[,2:(2+j)])[, ,2])/as.vector(table(as.data.frame(D)[,2:(1+j)]))
    }
    all.prop <- rbind(all.prop, prop)


    D <- as.matrix(D)

    #calculate weighted optimality of each hypothesis
    lopt <- calc.wL(D, cr, prior.scale, wr)
    dawopt <- calc.wDA(D, cr, prior.scale, wr)
   # daopt <- calc.wDA(D, t(matrix(cr[1:(r-1),])), prior.scale, wr=NULL)

    all.lopt <- c(all.lopt,lopt)
    all.dawopt <- c(all.dawopt, dawopt)
    #all.daopt <- c(all.daopt, daopt)
    cat(i, "\n")
    print(D)
    cat(beta, "\n")
    cat(wr, "\n")




  }


  D <- data.frame(D)


  results <- list(D=D, y=y, all.beta=all.beta, all.wr=all.wr,  all.sd = all.sd, all.emp.sd= all.emp.sd, beta = beta,
                  all.lopt=all.lopt, all.dawopt=all.dawopt, #
                  #all.daopt = all.daopt,
                  all.prop=all.prop, all.alpha=all.alpha,
                  all.power.emp = all.power.emp, all.power.th = all.power.th)

  return(results)

}


