



#' Find the empirical distribution of the covariates
#' @param D design matrix
#' @param all.z dataframe containing all covariate combinations
#' @param j number of covariates
#'
#' @return dataframe z.probs with distribution of covariates
#'
#'
#' @export
learn.zprobs <- function(D, all.z, j){
  occ <- row.match(as.data.frame(D[,2:(j+1)]), all.z)
  freq <- data.frame(table(occ)/length(occ))
  z.probs <- cbind(1:(2^j), rep(0, 2^j))
  z.probs[as.numeric(levels(freq$occ)), 2] <- freq[,2]
  z.probs <- z.probs[,2]
   #only work a single dynamic covar
  return(z.probs)

}



#' Assuming the currrent response, future covariate value and future treatment,
#' Calculate optimality if horizon is 1. If not, iterate back to exp.loss function.
#'
#' @param z.next vector of covariate values for future unit
#' @param t.next treatment of future unit
#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param zp vector of probabilities for each level of covariate z (needs to in the same order as all.z)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param int set to NULL if there are no interactions, set to T of there are interactions
#' @param lossfunc  the objective function to minimize
#' @param dyn set to NULL of there are no dynamic covariates, set to a real number if there are dynamic covariates
#' @param ... further arguments to be passed to <lossfunc>

#' @return value of objective function one step ahead in the future, assuming z.next and t.next
#'
#'
#' @export
#'


future.loss.k <- function(z.next, t.next, z.now, t.now, zp, N, design, int, lossfunc, dyn=NULL,  ...){

  if (is.null(int)){
    # for the case of no interactions
    # add row for current and future design points to the design matrix
    design <- rbind(design, c(1, z.now, t.now), c(1, z.next, t.next))
  } else {
    #assuming interactions
    design <- rbind(design, c(1, z.now, t.now, z.now*t.now), c(1, z.next, t.next, z.next*t.next))  # add row for current and future design points to the design matrix
  }

  #if k==1, calculate the of the above design matrix.
  #Else, calculate expected loss going one step ahead into the future. #regularization used
 if (N==1) { loss <- lossfunc(design, ...)
 } else {

    if (!is.null(dyn)){

      dyn <- dyn+1
    }


    loss <- exp.loss.k(z.now=z.next, t.now=t.next, zp, N-1, design=rbind(design, c(1, z.now, t.now)), int, lossfunc, dyn, ...)
    }

  return(loss)
}




#' Break down the expected future objective function by cases for every combination of:
#' 1) future possible covariate
#' 2) future possible treatment
#' Find a weighted average across all these cases

#' @param z.now vector of covariate values for current unit
#' @param t.now treatment of current unit
#' @param zp vector of probabilities for each level of covariate z (needs to in the same order as all.z)
#' @param N natural number greater than 0 for horizon
#' @param design design matrix constructed for all units up until the current unit
#' @param int set to NULL if there are no interactions, set to T of there are interactions
#' @param lossfunc  the objective function to minimize
#' @param dyn set to NULL of there are no dynamic covariates, set to T if there are dynamic covariates
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return design matrix D
#'
#'
#' @export
#'

exp.loss.k <- function(z.now, t.now, zp, N, design, int, lossfunc, dyn=NULL,  ...){

  #all possible combinations of levels for k binary covariates
  if(!is.null(int)){
    j  <- (ncol(design)-2)/2
  }else{
    j  <- ncol(design)-2
  }
  all.z <- expand.grid(rep(list(c(-1,1)),j))
  names(all.z) <- NULL


  loss.p <- apply(all.z, 1, future.loss.k, t.next=1, z.now, t.now, zp, N, design, int, lossfunc , dyn,  ...)#loss for each covariate combination when the future patient has tmt=1
  loss.m <- apply(all.z, 1, future.loss.k, t.next=-1, z.now, t.now, zp, N, design,int, lossfunc , dyn,  ...) #loss for each covariate combination when the future patient has tmt=-1

  loss <- ifelse(loss.p < loss.m, loss.p, loss.m)   #optimal value loss for each covariate combination when future patient has z=z.now

  if (!is.null(dyn)){
    exploss <- sum(loss*zp[dyn,])
  }else{
    exploss <- sum(loss*zp)

  }

 #expected loss: weighed according to probabilities of each level of z

  return(exploss)
}



#' Allocate treatments according to an information matrix based optimality criterion allowing for a non-myopic approach.
#' We assume a linear model for the response and simulate responses sequentially.
#' @param covar a dataframe for the covariates
#' @param init the number of units in the initial design
#' @param z.probs probability of each covariate being equal to 1
#' @param k integer for number of "outer" loops in coordinate exchange algorithm for initial design
#' @param N natural number greater than 0 for horizon
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc the objective function to minimize
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return Design matrix D, all estimates of beta, final estimate of beta, responses y
#'
#'
#' @export
linear.nonmyop <- function(covar, init, z.probs, k=NULL, N, int=NULL , lossfunc=calc.D, stoc, ... ){


  n <- nrow(covar)
  k <- ncol(covar)
  names(covar) <- NULL
  opt <- c()
  #initial starting design with init patients. Constructed using coordinate exchange algorithm with D-A optimality

  design <- coordex(as.data.frame(covar[1:init,]), k, lossfunc, int, ...)

  if (!is.character(z.probs)){
      zp <- zpfunc(z.probs)
  }

  for (i in (init+1):(n-1)){

    if (z.probs[1]=="learn"){
      z.probs <- learn.zprobs(D=design, all.z = expand.grid(rep(list(c(-1,1)),k)) , k)
      zp <- zpfunc(z.probs)
    }


    #allocate n-init patients by choosing treatment which minimizes expected loss
    eloss.p <- exp.loss.k(z.now=as.numeric(covar[i,]), t.now=1, zp=zp, N, design, int, lossfunc, ...) #expected loss for treatment 1
    eloss.m <- exp.loss.k(z.now=as.numeric(covar[i,]), t.now=-1, zp=zp, N, design, int, lossfunc,  ...)#expected loss for treatment -1


    p <- (1/eloss.p)/sum(1/eloss.p+1/eloss.m)

    if (is.na(p)){
      p <- 0.5

    }

    if(stoc==T){
      new.tmt <- sample(c(-1,1), 1, prob=c(1-p, p))         #Assign treatments
    }else{
      if (eloss.p > eloss.m) {
        new.tmt <- -1
      } else if (eloss.p < eloss.m) {
        new.tmt <- 1
      } else if (eloss.p == eloss.m) {
        new.tmt <- sample(c(-1,1), 1)
      }
    }


    if (is.null(int)){
      design <- rbind(design, c(1, as.numeric(covar[i,]), new.tmt)) #return design
    } else {
      design <- rbind(design, c(1, as.numeric(covar[i,]), new.tmt,as.numeric(covar[i,])*new.tmt))
    }


    opt <- c(opt, lossfunc(design , ...))

  }

  D <- design

  colnames(D) <- c("intercept", rep("covar", k), "tmt")
  return(list(D=D, opt=opt))

}




#' Allocate treatments according to an information matrix based optimality criterion allowing for a non-myopic approach.
#' Allow for dynamic covariates. We assume a linear model for the response and simulate responses sequentially.
#' @param covar a dataframe for the covariates
#' @param init the number of units in the initial design
#' @param z.probs vector of probabilities for each level of covariate z
#' @param k integer for number of "outer" loops in coordinate exchange algorithm for initial design
#' @param N natural number greater than 0 for horizon
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param lossfunc the objective function to minimize
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return Design matrix D
#'
#' @export

linear.nonmyop.dyn <- function(covar, init, z.probs, k, N, int=NULL , lossfunc=calc.D, stoc=F, ... ){


  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL
  opt <- c()
  #initial starting design with init patients. Constructed using coordinate exchange algorithm with D-A optimality

  design <- coordex(as.data.frame(covar[1:init,]), k, lossfunc, int, ...)

  if (!is.character(z.probs)){
    zp <- zpfunc(z.probs)
  }


  for (i in (init+1):(n)){


    if (z.probs[1]=="learn"){
      zp <- learn.zprobs(D=design, all.z = expand.grid(rep(list(c(-1,1)),j)) , j)

    Nprobs <- as.data.frame(matrix((rep(zp, N)), byrow=T, nrow=N))
    }else{

    Nprobs <-data.frame(zp[(i+1):(i+N),])
    }


    dyn=1




    #allocate n-init patients by choosing treatment which minimizes expected loss
    eloss.p <- exp.loss.k(z.now=as.numeric(covar[i,]), t.now=1, zp=Nprobs, N, design, int, lossfunc, dyn, ...) #expected loss for treatment 1
    eloss.m <- exp.loss.k(z.now=as.numeric(covar[i,]), t.now=-1, zp=Nprobs, N, design, int, lossfunc, dyn, ...)#expected loss for treatment -1

    #opt.tmt <- ifelse(round(eloss.p, 5) < round(eloss.m, 5),  1, -1) #pick best, get rid of variability

    p <- (1/eloss.p)/sum(1/eloss.p+1/eloss.m)

    if (is.na(p)){
      p <- 0.5

    }

    if(stoc==T){
      new.tmt <- sample(c(-1,1), 1, prob=c(1-p, p))         #Assign treatments
    }else{
      if (eloss.p > eloss.m) {
        new.tmt <- -1
      } else if (eloss.p < eloss.m) {
        new.tmt <- 1
      } else if (eloss.p == eloss.m) {
        new.tmt <- sample(c(-1,1), 1)
      }
    }

    if (is.null(int)){
      design <- rbind(design, c(1, as.numeric(covar[i,]), new.tmt)) #return design
    } else {
      design <- rbind(design, c(1, as.numeric(covar[i,]), new.tmt,as.numeric(covar[i,])*new.tmt))
    }


    opt <- c(opt, lossfunc(design , ...))
  }

  D <- design

  colnames(D) <- c("intercept", rep("covar", k), "tmt")
  return(list(D=D, opt=opt))

}

