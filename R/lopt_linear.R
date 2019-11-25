library("MASS")
library("arm")
library("mvtnorm")

#' Compute the L-optimal objective function assuming a linear model (with or without weights)
#' @param D design matrix
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#'
#' @return Value of the L-optimal objective function
#'
#' @export
calc.wL <- function(D, cr, prior.scale=100, wr=NULL){
  r <- nrow(cr)
  var <- rep(0, r)

  for (j in 1:r){
    if(is.null(wr)){
      var[j] <- t(cr[j,])%*% solve(t(D)%*%D + diag(1/prior.scale ,ncol(D))) %*% cr[j,]
    }else{
      var[j] <-  as.matrix(wr[j])%*%t(cr[j,])%*% solve(t(D)%*%D + diag(1/prior.scale ,ncol(D))) %*% cr[j,]

    }
  }
  opt <- sum(var)

  return(opt)
}







#' Compute the L-optimal objective function assuming a linear model (with or without weights)
#' @param D design matrix
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if equal weights
#'
#' @return Value of the L-optimal objective function
#'
#' @export
calc.vars <- function(D, cr, prior.scale=100){
  r <- nrow(cr)
  var <- rep(0, r)
  
  for (j in 1:r){
      var[j] <- t(cr[j,])%*% solve(t(D)%*%D + diag(1/prior.scale ,ncol(D))) %*% cr[j,]

    }
  
 
  
  return(var)
}






#' Compute the DA-optimal objective function assuming a linear model (with or without weights)
#' @param D design matrix
#' @param cr  matrix of contrasts
#' @param prior.scale prior scale parameter
#' @param wr matrix of weights, set to NULL if weights equal
#'
#' @return Value of the DA-optimal objective function
#'
#' @export
calc.wDA <- function(D, cr, prior.scale=100, wr=NULL){
  r <- nrow(cr)

  if(is.null(wr)){
      DA <- det( (cr)%*% solve(t(D)%*%D + diag(1/prior.scale ,ncol(D))) %*% t(cr))
  }else{

    var <-rep(0, r)

    inv <- solve(t(D)%*% D + diag(1/prior.scale, ncol(D)))
    for (m in 1:r){
      var[m] <-   as.matrix( t(cr[m,]) %*% inv %*% as.matrix((cr[m,]))) ^ wr[m]
    }
   DA <- prod(var)

  }

  return(DA)

}


#' Allocate treatment using a coordinate exchange algorithm with an optimality criterion for linear combinations of parameters.
#' Contrasts are assumed to be differences in treatment effect for each subgroup.
#' @param covar dataframe of covarite values
#' @param k integer for the number of "outer loop"
#' @param cr  matrix of contrasts
#' @param cr.lossfunc a loss function appropriate for linear combinations of parameters.
#' @param prior.scale prior scale parameter
#'
#' @return Design matrix
#'
#' @export
coord.cr <- function(covar, k, cr, cr.lossfunc, prior.scale=100, wr=NULL){

  n <- nrow(covar) #covar  must be a dataframe
  j <- ncol(covar)

  Dms <- list() #list to store k potential designs

  for (m in 1:k){

    Dm <- as.matrix(cbind(rep(1, n), covar, tmt=sample(c(0,1), n, replace=T))) #design matrix with random treatment assignment
    Dm <- cbind(Dm, Dm[, 2:(j+1)]*Dm[,(j+2)])


    repeat{
      D <- Dm
      for (i in 1:n){
        D[i, "tmt"] <- 1   #calculate criterion for unit i assigned to tmt 1
        D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]

        plusD <- cr.lossfunc(D, cr, prior.scale, wr)

        D[i, "tmt"] <- 0 #calculate criterion for unit i assigned to tmt -1
        D[i, (j+3):(2*j+2)] <- rep(0, j)

        minusD <- cr.lossfunc(D, cr, prior.scale, wr)

        if (plusD < minusD){
          D[i,"tmt"] <- 1
          D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]
        } else {
          D[i,"tmt"] <- 0
          D[i, (j+3):(2*j+2)] <- rep(0, j)
        }
      }

      if (identical(D, Dm)) break #if no change to resulting design matrix, terminate
      else Dm <- D

    }

    Dms[[m]] <- D




  }


  mindes <- which.min(unlist(lapply(Dms, cr.lossfunc, cr, prior.scale)))  #find the optimum design matrix

  result <- as.data.frame(Dms[mindes][[1]])


  return(result)
}



#' Assuming a linear model for the response, allocate treatment sequentially based on an optimality criterion for
#' linear combinations of parameters. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param true.sigma the true parameter value for the standard deviation
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param cr.lossfunc loss function appropriate for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param u vector of uniform random numbers for generating responses. If set to NULL, responses generated from the binomial distribution.
#' @param ... further arguments to be passed to <lossfunc>
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
cr.des <- function(covar, true.beta, true.sigma, threshold,  kappa, init, cr.lossfunc, k, wt, int=T, prior.scale=100,
                     same.start=NULL, rand.start=NULL, stoc=T, bayes=T, prior.default=T, u=NULL, ... ){

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
 


  D <- as.matrix(D)


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
  #ll.daopt <- calc.wDA(D, t(matrix(cr[1:(r-1),])), prior.scale, wr=NULL)
  #if full cr matrix is used, we have numerical issues due to matrix being full column rank




  #initializing done

  #iterate
  for (i in (init+1):n){

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions
    if (!is.null(int)){
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1, as.numeric(covar[i, ]) ))
    } else{
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
    }


    if (wt!=T){
      wr=NULL
    }
    d.plus <-  cr.lossfunc(mat.plus, cr, prior.scale, wr)

    #design matrix if tmt 0 is assigned
    if (!is.null(int)){
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0,j) ))
    } else{
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0)) #Design matrix with new treatment =0
    }

    d.minus <- cr.lossfunc(mat.minus, cr, prior.scale, wr)

    if (d.plus <0){
      d.plus=0
    }

    if (d.minus <0){
      d.plus=0
    }


    if (stoc==T){

      if (d.plus==d.minus){
        probs = 0.5
      } else {
        probs <- (1/d.plus)/(1/d.plus+1/d.minus)
      }

      new.tmt <- sample(c(0,1), 1, prob=c(1-probs, probs))         #Assign treatments

    }else{
      if(d.plus < d.minus){
        new.tmt <- 1
      }else if (d.plus > d.minus){
        new.tmt <- 0
      }else{
        new.tmt <- sample(c(0,1), 1)
      }
    }


    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else {
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL


    if (!is.null(u)){
      new.y <- new.d%*%true.beta + u[i]
    }else{
      new.y <-  rnorm(1, new.d%*%true.beta, true.sigma)                               #Simulate new observation

    }

    y <- c(y, new.y)                               #Simulate new observation


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
    #daopt <- calc.wDA(D, t(matrix(cr[1:(r-1),])), prior.scale, wr=NULL)

    all.lopt <- c(all.lopt,lopt)
    all.dawopt <- c(all.dawopt, dawopt)
    #all.daopt <- c(all.daopt, daopt)

  }


  D <- data.frame(D)


  results <- list(D=D, y=y, all.beta=all.beta, all.wr=all.wr,  all.sd = all.sd, all.emp.sd= all.emp.sd, beta = beta,
                  all.lopt=all.lopt, all.dawopt=all.dawopt,
                  #all.daopt = all.daopt,
                  all.prop=all.prop,all.alpha=all.alpha,
                  all.power.emp = all.power.emp, all.power.th = all.power.th)

  return(results)

}





#' Assuming a Bayesian linear model for the response with a normal inverse gamma prior,
#' allocate treatment sequentially based on an optimality criterion for
#' linear combinations of parameters. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param sigma the true parameter value for the standard deviation
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param a shape parameter for the inverse gamma distribution for the prior, by default set to 2
#' @param d scale parameter for the inverse gamma distribution for the prior, by default set to NULL for a zero vector
#' @param M mean vector for the normal distribution NULL,
#' @param cr.lossfunc loss function appropriate for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param rand.start If set to T, function generates an initial design randomly. Else, coordinate exchange is used.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param bayes set to T if bayesglm is used instead of glm. Default prior assumed.
#' @param prior.default set to T if default priors for bayesglm is used. If set to False and bayes=T, normal priors used.
#' @param ... further arguments to be passed to <lossfunc>
#'
#'
#'
#'
#' @return design matrix D, responses y, all estimates of M, all estimates of V, all estimates of a, all estimates of d,
#' all weights, determinant of V, trace of V,
#'  probabilities for treatment assignment, all values of the objective functions (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#
#' @export

cr.bayes.des <- function(covar, true.beta, sigma, threshold,  kappa, init, cr.lossfunc, k, wt=T, int=T, a=2, d=2, M=NULL, prior.scale=100,
                         same.start=NULL, rand.start=NULL, stoc=T, bayes=T, prior.default=T, ... ){


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  tvalue <- qt(0.95, n-p-1)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)

  r <- nrow(cr)

  if(is.null(M)){
    M <- matrix(rep(0, length(true.beta)))
  }

  V <- diag(prior.scale, length(true.beta))


  inv.V <- solve(V)

  # use L-optimal design with equal weights for initial design
  # use L-optimal design with equal weights for initial design
  if (!is.null(same.start)) {
    D <-same.start
  }else if (!is.null(rand.start)) {
    D <- cbind(rep(1, init), covar[1:init,], sample(c(-1,1), init, replace=T))
  }else{
    D <-  as.matrix(coord.cr(as.data.frame(covar[1:init,]), k, cr, cr.lossfunc, prior.scale))

  }

  D <- as.data.frame(D)

  #record proportion of patients in new treatment
  if (j==1){
    prop <-c(table(D[,2:(2+j)])[3]/ (table(D[,2:(2+j)])[1] + table(D[,2:(2+j)])[3]),
             table(D[,2:(2+j)])[4]/ (table(D[,2:(2+j)])[2] + table(D[,2:(2+j)])[4]))
  }else if (dim ( table(D[,2:(2+j)]))[3]<2) {
    prop <- rep(0, r)
  } else {  prop <- as.vector(table(D[,2:(2+j)])[, ,2])/as.vector(table(D[,2:(1+j)]))
  }
  #set up vector
  all.prop <- prop

  covnum <- table(covar[1:init,])
  all.covnum <- covnum

  D <- as.matrix(D)
  Dms <- list()

  y <- rnorm( n=init, mean=D%*%true.beta, sd=sigma) #generate init observations




  M.new <- solve(inv.V + t(D)%*%D)%*%(inv.V%*%M + t(D)%*%y)
  V.new <- solve(inv.V + t(D)%*%D)
  a.new <- a + t(M)%*%inv.V%*%M + t(y)%*%y-t(M.new)%*%solve(V.new)%*%M.new
  d.new <- d + length(y)


  M <- M.new
  V <- V.new
  a <- as.numeric(a.new)
  d <- as.numeric(d.new)
  #append to matrix of all estimates
  all.M <- t(M)
  all.V <- V
  all.a <- a
  all.d <- d

  trace.V <- sum(diag(V))
  det.V <- det(V)

  n.sim <- 1000

  #calculate weights
  sim <- rmvt(n=n.sim, sigma=a*V/d, df=d) + matrix(rep(t(M),each=n.sim),nrow=n.sim)
  crbeta <- sim %*% t(cr)
  pr <- colSums(crbeta < matrix(rep(t(threshold),each=n.sim),nrow=n.sim))/n.sim
  wr <- t(ifelse( pr >= kappa, pr, 0))
  all.wr <- wr


  #calculate alpha
  alpha <- cr %*% M < threshold
  all.alpha <- t(alpha)
  #calculate power
  tvalue <- qt(0.05, init-p)
  power.emp <- (cr %*% M - threshold) /(sqrt(diag( (a/(d-2)) * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))) < tvalue
  all.power.emp <- as.numeric(power.emp)

  power.th <- pt(tvalue + (threshold -cr %*% true.beta) /(sqrt(diag(sigma^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))), init-p)
  all.power.th <- as.numeric(power.th)



  #calculate weighted optimality of each hypothesis
  all.lopt <- calc.wL(D, cr, prior.scale, wr)
  #all.dawopt <- calc.wDA(D, cr, prior.scale, wr)
  #all.daopt <- calc.wDA(D, t(as.matrix(cr[1:(r-1),])), prior.scale, wr=NULL)
  #if full cr matrix is used, we have numerical issues due to matrix being full column rank

 vars <- calc.vars(D, cr)
  #initializing done

  #iterate
  for (i in (init+1):n){
    var.plus <- var.minus <- rep(0, r)

    #design matrix if tmt 1 is assigned - split into cases including/excluding interactions
    if (!is.null(int)){
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1, as.numeric(covar[i, ]) ))
    } else{
      mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
    }


    if (wt!=T){
      wr=NULL
    }
    d.plus <-  cr.lossfunc(mat.plus, cr, prior.scale, wr)


    #design matrix if tmt 0 is assigned
    if (!is.null(int)){
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0, rep(0,j) ))
    } else{
      mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), 0)) #Design matrix with new treatment =0
    }

    d.minus <- cr.lossfunc(mat.minus, cr, prior.scale, wr)

    if (d.plus <0){
      d.plus=0
    }

    if (d.minus <0){
      d.plus=0
    }


    if (stoc==T){

      if (d.plus==d.minus){
        probs = 0.5
      } else {
        probs <- d.minus/(d.plus+d.minus)
      }

      new.tmt <- sample(c(0,1), 1, prob=c(probs, 1-probs))         #Assign treatments

    }else{
      if(d.plus < d.minus){
        new.tmt <- 1
      }else if (d.plus > d.minus){
        new.tmt <- 0
      }else{
        new.tmt <- sample(c(0,1), 1)
      }
    }


    #new row of design matrix
    if (!is.null(int)){
      new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
    } else {
      new.d <- as.numeric(c(1, covar[i, ], new.tmt))
    }


    D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix
    row.names(D)<- NULL
    y <- c( y, rnorm(1, new.d%*%true.beta, sigma))                               #Simulate new observation



    inv.V <- solve(V)
    M.new <- solve(inv.V + t(D)%*%D)%*%(inv.V%*%M + t(D)%*%y)
    V.new <- solve(inv.V + t(D)%*%D)
    a.new <- a + t(M)%*%inv.V%*%M + t(y)%*%y-t(M.new)%*%solve(V.new)%*%M.new
    d.new <- d + length(y)


    M <- M.new
    V <- V.new
    a <- as.numeric(a.new)
    d <- as.numeric(d.new)
    #append to matrix of all estimates
    all.M <- rbind(all.M, t(M))
    all.V <- list(all.V, V)
    all.a <- c(all.a, a)
    all.d <- c(all.d , d)

    trace.V <- c(trace.V, sum(diag(V)))
    det.V <- c(det.V, det(V))

    #calculate weights
    sim <- rmvt(n=n.sim, sigma=a*V/d, df=d) + matrix(rep(t(M),each=n.sim),nrow=n.sim)
    crbeta <- sim %*% t(cr)
    pr <- colSums(crbeta < matrix(rep(t(threshold),each=n.sim),nrow=n.sim))/n.sim
    wr <- t(ifelse( pr >= kappa, pr, 0))
    all.wr <- rbind(all.wr, wr)

    #calculate alpha
    alpha <- cr %*% M < threshold
    all.alpha <- rbind(all.alpha, t(alpha))

    tvalue <- qt(0.05, i-p)
    power.emp <- (cr %*% M - threshold) /(sqrt(diag((a/(d-2)) * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))) < tvalue
    all.power.emp <- rbind(all.power.emp, as.numeric(power.emp))

    power.th <- pt(tvalue + (threshold -cr %*% true.beta) /(sqrt(diag(sigma^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))), i-p)
    all.power.th <- rbind(all.power.th, as.numeric(power.th))

    D <- as.data.frame(D)
    #proportion of patients with new treatment
    if (j==1){
      prop <-c(table(D[,2:(2+j)])[3]/ (table(D[,2:(2+j)])[1] + table(D[,2:(2+j)])[3]),
               table(D[,2:(2+j)])[4]/ (table(D[,2:(2+j)])[2] + table(D[,2:(2+j)])[4]))
    } else if (dim ( table(as.data.frame(D)[,2:(2+j)]))[3]<2) {
      prop <- rep(0, r)
    } else {  prop <- as.vector(table(as.data.frame(D)[,2:(2+j)])[, ,2])/as.vector(table(as.data.frame(D)[,2:(1+j)]))
    }
    all.prop <- rbind(all.prop, prop)
    D <- as.matrix(D)


    covnum <- table(covar[1:i,])
    all.covnum <- rbind(all.covnum, covnum)

    lopt <- calc.wL(D, cr, prior.scale, wr)
    #dawopt <- calc.wDA(D, cr, prior.scale, wr)
    #daopt <- calc.wDA(D, t(as.matrix(cr[1:(r-1),])), prior.scale, wr=NULL)
    
    vars <- rbind(vars, calc.vars(D, cr))

    all.lopt <- c(all.lopt,lopt)
    #all.dawopt <- c(all.dawopt, dawopt)
    #all.daopt <- c(all.daopt, daopt)

  }


  D <- data.frame(D)


  results <- list(D=D, Dms=Dms, y=y, all.M=all.M, all.V=all.V, all.a=all.a, all.d=all.d,
                  all.wr=all.wr,   all.prop=all.prop, all.covnum, det.V=det.V, trace.V =trace.V,
                  all.alpha=all.alpha,     all.lopt=all.lopt, 
                  #all.dawopt=all.dawopt, all.daopt = all.daopt,
                  all.power.emp = all.power.emp, all.power.th = all.power.th, vars=vars)

  return(results)

}




#' Assuming a Bayesian linear model for the response with a normal inverse gamma prior,
#' allocate treatment sequentially based on an optimality criterion for
#' linear combinations of parameters. Responses are simulated assuming the true parameter values.
#'
#' @param covar a dataframe for the covariates
#' @param true.beta the true parameter values of the regression coefficients
#' @param sigma the true parameter value for the standard deviation
#' @param threshold the cut-off value for hypothesis tests
#' @param kappa the value of probability at which weights are set at zero
#' @param init the number of units in the initial design
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param a shape parameter for the inverse gamma distribution for the prior, by default set to 2
#' @param d scale parameter for the inverse gamma distribution for the prior, by default set to NULL for a zero vector
#' @param M mean vector for the normal distribution NULL,
#' @param cr.lossfunc loss function appropriate for linear combinations of parameters
#' @param k the number of "outer loops" in the coordinate exchange algorithm
#' @param wt set to T if the above lossfunction is weighted, NULL otherwise
#' @param prior.scale the prior scale parameter
#' @param ... further arguments to be passed to <lossfunc>
#'
#'
#'
#'
#' @return design matrix D, responses y, all estimates of M, all estimates of V, all estimates of a, all estimates of d,
#' all weights, determinant of V, trace of V,
#'  probabilities for treatment assignment, all values of the objective functions (weighted L, DA, weighted DA), proportion of treatment=1,
#'  proportion of covariate in each group
#'  Type 1 error, true value of power, empirical value of power.
#
#' @export

nonseq.cr.bayes.des <- function(covar, true.beta, sigma, threshold,  kappa, init, cr.lossfunc, k, wt=T, int=T, a=2, d=2, M=NULL, prior.scale=100,
                         prior.default=T, ... ){
  
  
  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  p <- length(true.beta)
  tvalue <- qt(0.95, n-p-1)
  cr1 <- matrix(0, 2^j, j+1)
  cr2 <- expand.grid(rep(list(c(0,1)),j))
  cr <- cbind(cr1, 1, cr2)
  cr <- as.matrix(cr)
  
  r <- nrow(cr)
  
  if(is.null(M)){
    M <- matrix(rep(0, length(true.beta)))
  }
  
  V <- diag(prior.scale, length(true.beta))
  
  
  inv.V <- solve(V)
  
  
  n.sim <- 1000
  
  #calculate weights
  sim <- rmvt(n=n.sim, sigma=a*V/d, df=d) + matrix(rep(t(M),each=n.sim),nrow=n.sim)
  crbeta <- sim %*% t(cr)
  pr <- colSums(crbeta < matrix(rep(t(threshold),each=n.sim),nrow=n.sim))/n.sim
  wr.initial <- t(ifelse( pr >= kappa, pr, 0))
  
  
  D <- coord.cr(covar, k, cr, cr.lossfunc, prior.scale=100, wr=wr.initial)
  D <- as.matrix(D)
  
    row.names(D)<- NULL
    y <-  rnorm(nrow(D), D%*%true.beta, sigma)                              #Simulate new observation
    
    
    
    inv.V <- solve(V)
    M.new <- solve(inv.V + t(D)%*%D)%*%(inv.V%*%M + t(D)%*%y)
    V.new <- solve(inv.V + t(D)%*%D)
    a.new <- a + t(M)%*%inv.V%*%M + t(y)%*%y-t(M.new)%*%solve(V.new)%*%M.new
    d.new <- d + length(y)
    
    
    M <- M.new
    V <- V.new
    a <- as.numeric(a.new)
    d <- as.numeric(d.new)
    
    
    trace.V <- sum(diag(V))
    det.V <- det(V)
    
    #calculate weights
    sim <- rmvt(n=n.sim, sigma=a*V/d, df=d) + matrix(rep(t(M),each=n.sim),nrow=n.sim)
    crbeta <- sim %*% t(cr)
    pr <- colSums(crbeta < matrix(rep(t(threshold),each=n.sim),nrow=n.sim))/n.sim
    wr <- t(ifelse( pr >= kappa, pr, 0))

    #calculate alpha
    alpha <- cr %*% M < threshold

    tvalue <- qt(0.05, n-p)
    power.emp <- (cr %*% M - threshold) /(sqrt(diag((a/(d-2)) * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))) < tvalue

    power.th <- pt(tvalue + (threshold -cr %*% true.beta) /(sqrt(diag(sigma^2 * cr%*% solve(t(D) %*%D+diag(1/prior.scale, ncol(D))) %*% t(cr)))), n-p)

    lopt <- calc.wL(D, cr, prior.scale, wr)
    
    vars <- calc.vars(D, cr)
    
    
    #dawopt <- calc.wDA(D, cr, prior.scale, wr)
    #daopt <- calc.wDA(D, t(as.matrix(cr[1:(r-1),])), prior.scale, wr=NULL)
   
    #proportion of patients with new treatment
    if (j==1){
      prop <-c(table(D[,2:(2+j)])[3]/ (table(D[,2:(2+j)])[1] + table(D[,2:(2+j)])[3]),
               table(D[,2:(2+j)])[4]/ (table(D[,2:(2+j)])[2] + table(D[,2:(2+j)])[4]))
    } else if (dim ( table(as.data.frame(D)[,2:(2+j)]))[3]<2) {
      prop <- rep(0, r)
    } else {  prop <- as.vector(table(as.data.frame(D)[,2:(2+j)])[, ,2])/as.vector(table(as.data.frame(D)[,2:(1+j)]))
    }
  
  
  
  D <- data.frame(D)
  
  
  results <- list(D=D, y=y, M=M, V=V, a=a, d=d,
                  wr=wr,   prop=prop, det.V=det.V, trace.V =trace.V,
                  alpha=alpha,     lopt=lopt, 
                  #dawopt=dawopt, daopt = daopt,
                  power.emp = power.emp, power.th = power.th, vars=vars)
  
  return(results)
  
}



