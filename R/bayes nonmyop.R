#non-myopic logistic regression for weighted L-optimal design 

probi <- function(x, beta){
  #purpose: compute P(Y=1) for logistic regression
  #inputs: x is a row of the design matrix
  #output: beta are the coefficients
  
  x <- as.matrix(x)
  eta <- beta %*% x
  p <- (exp(eta)/(1+exp(eta)))  #logistic function
  
  return(p)
  
}






Imat.beta <- function(X, beta){
  #purpose: compute information matrix, given beta 
  #inputs: X is the design matrix 
  #        beta is the current estimate of the regression coefficients 
  #output: I.mat is the information matrix 
  
  X <- as.matrix(X)
  
  p <- apply(X, 1, probi, beta=beta)  #calculate p(Y=1) for each row of the design matrix
  W <- diag(p * (1 - p))   #weights
  I.mat <- t(X) %*% W %*% X
  
  return(I.mat)
  
}




future.loss <- function(z.next, t.next, z.probs, N, design, k,  beta, y, cr, wr, prior.scale, dyn=NULL, Nprobs=NULL, ...){
  #purpose: calculate loss for one timepoint in the future 
  #inputs:  
  #         z.next: vector covariate values for future patient being considered 
  #         t.next: treatment for future patient being considered  
  #         z.now: vector covariate values for current patient j+1
  #         t.now: proposed treatment for current patient   
  #         z.probs: vector of probabilities for each level of z (needs to in the same order as all.z below)
  #         N: natural number greater than 0 for horizon  
  #         design: design matrix constructed with j patients who have been allocated treatments 
  #         lossfunc: optimality criterion 
  #         beta: latest estimate of regression coefficients 
  #         y: responses observed or simulated
  #output:  calculated loss 
  
  #add new point to design matrix
  
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
  
  loss <- exp.loss( z.now=z.next, t.now=t.next, z.probs, N-1, design, k, beta, y, cr, wr, dyn, prior.scale, Nprobs, ...)
  
  
  
  }
  
  return(loss)
  
}


future.y <- function(y.now, z.now, t.now, z.probs, N, design, k, beta, y, cr, wr,  prior.scale, dyn=NULL, Nprobs=NULL, ...){
  #purpose: calculate loss given current design, current patient covariate, proposed treatment and N possible trajectory in the future 
  #         assuming model for response is logistic regression 
  #inputs:  
  #         z.next: vector covariate values for future patient being considered 
  #         t.next: treatment for future patient being considered  
  #         z.now: vector covariate values for current patient j+1
  #         t.now: proposed treatment for current patient   
  #         z.probs: vector of probabilities for each level of z (needs to in the same order as all.z below)
  #         N: natural number greater than 0 for horizon  
  #         design: design matrix constructed with j patients who have been allocated treatments 
  #         lossfunc: optimality criterion 
  #         beta: latest estimate of regression coefficients 
  #         y: responses observed or simulated
  #output:  calculated loss 
  
  
  
  d.now <- c(1, z.now, t.now, z.now*t.now)
  row.names(design) <-NULL
  
  design <- rbind(design, d.now)
  
  #append y.now to seen data
  sim.y <- c(y, y.now)
  
  
  
  
  
  
  sim.M <- solve(inv.V + t(D)%*%D)%*%(inv.V%*%M + t(D)%*%sim.y)
  sim.V <- solve(inv.V + t(D)%*%D)
  a.new <- as.numeric(a + t(M)%*%inv.V%*%M + t(sim.y)%*%sim.y-t(sim.M)%*%solve(sim.V)%*%sim.M)
  d.new <- as.numeric(d + length(sim.y))
  
  
 
    
  
  all.z <- expand.grid(rep(list(c(-1,1)),k))  #grid for all possible covariates
  names(all.z) <- NULL
  
  #for each possible covariate value, calculate loss when t.next=1
  loss.p <- apply(all.z, 1, future.loss , t.next=1, z.probs, N, design, k, sim.a, sim.d, sim.M, sim.V, sim.y, cr, wr, prior.scale, dyn, Nprobs, ...)#loss for each covariate combination when the future patient has tmt=1
  
  #for each possible covariate value, calculate loss when t.next=1
  loss.m <- apply(all.z, 1, future.loss , t.next=-1, z.probs, N, design, k,  sim.a, sim.d, sim.M, sim.V, sim.y, cr, wr, prior.scale, dyn, Nprobs, ...) #loss for each covariate combination when the future patient has tmt=-1
  
  #find loss for optimal treatment (t*) for each possible covariate value 
  loss <- ifelse(loss.p < loss.m, loss.p, loss.m)   
  
  #expected loss: weighed according to probabilities of each covariate value
  exploss <- sum(loss*z.probs)  
  
  
  return(exploss) 
  
}



exp.loss <- function(z.now, t.now, z.probs, N, design, k, beta, y, cr, wr,  prior.scale, dyn=NULL, Nprobs=NULL,...){
  #purpose: calculate expected loss given current design, current patient covariate and horizon 
  #inputs:  
  #         z.now: value of covariate Z for current patient j+1
  #         t.now: proposed treatment for current patient 
  #         z.probs: vector of probabilities for each level of z (needs to in the same order as all.z below)
  #         N: natural number greater than 0 for horizon  
  #         design: design matrix constructed with j patients who have been allocated treatments 
  #         lossfunc: optimality criterion 
  #         beta: latest estimate of regression coefficients 
  #         y: responses observed or simulated
  #output:  expected loss 
  
  
  #probability that y.now=1
  
  pi <- probi(as.numeric(c(1, z.now, t.now, z.now*t.now)), M) 
  
  #expected loss, given y.now=1
  loss1 <- future.y(y.now=1, z.now, t.now, z.probs, N, design, k, a, d, M, V, y,  cr, wr,  prior.scale,dyn, Nprobs, ...) 
  
  #expected loss, given y.now=1
  loss0 <- future.y(y.now=0, z.now, t.now, z.probs, N, design, k,  a, d, M, V, y, cr, wr,  prior.scale, dyn, Nprobs, ...) 
  
  #weighted by probability that P(y=1)
  exp.loss<- pi*loss1 + (1-pi)*loss0
  
  
  return(exp.loss)
}







logit.Lbnon <- function(covar,  true.beta, threshold, kappa, init, z.probs, N,  a=2, d=2, M=NULL, prior.scale =100, dyn=NULL, Nprobs=NULL, ...){
  #purpose: allocate treatments according to an information matrix based optimality criterion 
  #         when you have a logistic model for the response. We simulate responses sequentially.
  #         Allow for nonmyopic approach 
  #input:   dataframe covar of covariate values 
  #         beta0 a vector of initial parameter guesses
  #         init the number of observations to include in the preliminary design 
  #         N is the horizon 
  #         lossfunc is the optimality criterion 
  #         m: number of simulations of y to be generated in the expectation 
  #output:  design matrix D
  #         responses y 
  #         matrix of beta estimates, beta
  
  n <- nrow(covar)
  k <- ncol(covar) #covar must be a dataframe 
  
  
  
  cr1 <- matrix(0, 2^k, k+1)  
  cr2 <- expand.grid(rep(list(c(0,1)),k))
  cr <- as.matrix(cbind(cr1, 1, cr2))
  
  
  r <- nrow(cr)
  
  
  
  
  if(is.null(M)){
    M <- matrix(rep(0, length(true.beta)))
  }
  
  V <- diag(prior.scale, length(true.beta))
  
  
  inv.V <- solve(V)
  
  
  # randomly select treatment for first unit
  design <- as.matrix(cbind(1, covar=covar[1:init, ], tmt=sample(c(0, 1), init, replace=T))) 
  design <- cbind(design, design[,2:(k+1)]*design[,(k+2)])   #append interaction columns 
  row.names(design)<- NULL
  design <- as.matrix(design)
  Dms <- list()
  
  y <- rnorm( n=init, mean=design%*%true.beta, sd=sigma) #generate init observations
  
  
  
  
  M.new <- solve(inv.V + t(design)%*%design)%*%(inv.V%*%M + t(design)%*%y)
  V.new <- solve(inv.V + t(design)%*%design)
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
  
  
  #calculate weighted optimality of each hypothesis  
  inv <- solve(t(design)%*% design+ diag(prior.scale, length(true.beta)))
  opt <- rep(0, r)
  for (m in 1:r){
    opt[m] <- as.matrix(wr[m]%*% t(cr[m,]) %*% inv %*% as.matrix((cr[m,])))
  }
  all.opt <- sum(opt)
  
  for (i in (init+1):n){
    
    if (z.probs[1]=="learn"){
      z.probs <- learn.zprobs(design=design, all.z = expand.grid(rep(list(c(-1,1)),k)) , k)
    }
    
    
    #allocate treatment which minimizes expected loss
    eloss.p <- exp.loss(z.now=as.numeric(covar[i,]), t.now=1, z.probs, N, design, k, a, d, M, V, y, cr, wr, prior.scale,  ...) #expected loss for treatment 1 
    eloss.m <- exp.loss(z.now=as.numeric(covar[i,]), t.now=-1, z.probs, N, design, k, a, d, M, V, y, cr, wr,  prior.scale, ...)#expected loss for treatment -1
    
    probs <- eloss.m/(eloss.p+eloss.m) 
    
    #diagnostics
    
    #if (probs<0){
    # probs <-0
    
    #} else if (probs >1){
    # probs <- 1
    
    #}
    
    opt.tmt <- sample(c(-1,1), 1, prob=c(1-probs, probs))         #Assign treatments 
    
    
    #new row for design matrix 
    new.d <- as.numeric(cbind(1, data.frame(covar[i,]), opt.tmt))
    #new row if we allow for tmt-cov interaction:
    if (!is.null(int)){
      new.d <- c(new.d, new.d[2:(k+1)]*opt.tmt)   #append interaction columns 
    }
    
    design <- as.matrix(rbind(design, as.numeric(new.d)))     #Add the new row to the design matrix 
    
    pi <- probi(new.d, true.beta)       #Compute new pi
    y <- c( y, rbinom(1, 1, pi))   #Simulate new observation 
    
    beta <- coef(bayesglm (y~design[,-1], family=binomial(link="logit")))  #update beta
    all.beta <- rbind(all.beta, beta)                          #Store all betas 
    
  }
  
  design <- data.frame(design)
  
  results <- list(D=design, y=y, all.beta=all.beta, beta = beta)
  
  return(results)
  
}


