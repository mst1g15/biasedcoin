#design criteria


#' Given a design matrix for a linear model, compute the loss
#'
#' @param mat matrix or dataframe for the design matrix
#' @param epsilon a small real number used for regularization. If set to zero,
#' no regularization takes place
#' @param int set to NULL if there are no interactions, set to TRUE of there are interactions
#' @return loss of the matrix
#'
#' @export
calc.loss <- function(mat, epsilon=0.00001, int=NULL){

  n <- nrow(mat)
  j <- ncol(mat)
  if (exists("mat$tmt")){ tmt <-mat$tmt
  } else tmt <- mat[,j]
  Z<- as.matrix(mat[,-j])
  noise <- epsilon*diag(j-1)
  if (!is.null(int)){
    int <- matrix(rep(tmt, j-1), ncol=j-1, byrow=F)*Z
    noise <- epsilon*diag(j-1)
    loss <- det(solve(t(Z)%*%Z - t(int)%*%Z%*%solve(t(Z)%*%Z + noise)%*%t(Z)%*%int+ noise))

  } else{
    loss <- t(tmt)%*%Z%*%solve(t(Z)%*%Z + noise)%*%t(Z)%*%tmt
  }

  return(loss)
}


#' Given a design matrix for a linear model, compute the D-A optimal objective function
#'
#' @param mat matrix or dataframe for the design matrix
#' @param A a matrix where each column indicates the linear combination of parameters of interest
#' @param epsilon a small real number used for regularization
#' @return value of the D-A optimal objective function
#'
#' @export

calc.DA <- function(mat, A=t(c(0, 0, 1, 0)), epsilon=0.00001){
  j <- ncol(mat)
  mat <- as.matrix(mat)
  loss <- sum(diag(A%*%solve(t(mat)%*%mat + epsilon*diag(j))%*%t(A)))

  return(loss)
}





#' Given a design matrix for a linear model, compute the A-optimal objective function
#'
#' @param mat matrix or dataframe for the design matrix
#' @param epsilon a small real number used for regularization
#' @return value of the A optimal objective function
#'
#' @export

calc.A <- function(mat, epsilon=0.00001){

  j <- ncol(mat)
  mat <- as.matrix(mat)
  loss <- sum(diag(solve(t(mat)%*%mat + diag(j)*epsilon)))

  return(loss)
}




#' Given a design matrix for a linear model, compute the D-optimal objective function
#'
#' @param mat matrix or dataframe for the design matrix
#' @param epsilon a small real number used for regularization
#' @return value of the D optimal objective function
#'
#' @export
calc.D <- function(mat, epsilon=0.00001){


  j <- ncol(mat)
  mat <- as.matrix(mat)
  M <- t(mat)%*%mat + diag(j)*epsilon
  loss <- 1/det(M)

  return(loss)
}




#' Given a design matrix for a linear model, compute the G-optimal objective function
#'
#' @param mat matrix or dataframe for the design matrix
#' @param epsilon a small real number used for regularization
#' @return value of the G optimal objective function
#'
#' @export


calc.G <- function(mat, epsilon=0.00001){

  k <- ncol(mat)
  mat <- as.matrix(mat)
  X <- as.matrix(cbind(1, expand.grid(rep(list(c(-1,1)),k-1))))
  all <- diag((X)%*%solve(t(mat)%*%mat + diag(k)*epsilon)%*%t(X))
  loss <- max(diag((X)%*%solve(t(mat)%*%mat + diag(k)*epsilon)%*%t(X)))


  #return(list(loss=loss, all=all))

  return(loss)
}






#' Given a design matrix for a linear model wuth treatment-covariate interaction included, compute the G-optimal objective function
#'
#' @param mat matrix or dataframe for the design matrix
#' @param epsilon a small real number used for regularization
#' @return value of the G optimal objective function
#'
#' @export


calc.Gint <- function(mat, epsilon=0.00001){

  k <- ncol(mat)
  ncovar <- (k-2)/2
  mat <- as.matrix(mat)
  main <- as.matrix(expand.grid(rep(list(c(-1,1)),ncovar +1)))
  X <- cbind(1, main, main[,1:(ncovar)]*main[,-c(1:ncovar)])
  all <- diag((X)%*%solve(t(mat)%*%mat + diag(k)*epsilon)%*%t(X))
  loss <- max(diag((X)%*%solve(t(mat)%*%mat + diag(k)*epsilon)%*%t(X)))

    return(loss)

}


#' Given a design matrix for a linear model, calculate the covariate balance
#'
#' @param mat matrix or dataframe for the design matrix
#' @return covariate balance of the matrix
#'
#' @export

calc.tz <- function(mat){

  n <- nrow(mat)
  j <- ncol(mat)
  tmt <- as.matrix(mat$tmt)
  Z<- as.matrix(mat[,-j])

  result <- sum(abs(t(tmt)%*%Z))

  return(result)
}




#' Allocate treatment using coordinate exchange algorithm, assuming a linear model for the response
#' which has an intercept, an effect for each covariate, an effect for the treatment, and allows for
#' covariate-treatment interactions
#'
#'
#' @param covar a dataframe for the covariates
#' @param k the number of random starting designs to use
#' @param lossfunc the objective function to minimize
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param ... Further arguments to be passed to <lossfunc>

#'
#' @return design matrix
#'
#' @export
#'

coordex <- function(covar, k, lossfunc, int=NULL, ... ){

  #find the number of patients, n
  n <- nrow(covar) #covar  must be a dataframe
  j <- ncol(covar)

  Dms <- list() #list to store k potential designs

  for (m in 1:k){

    Dm <- as.matrix(cbind(rep(1, n), covar, tmt=sample(c(-1,1), n, replace=T))) #design matrix with random treatment assignment
    if (!is.null(int)){
      Dm <- cbind(Dm, Dm[,-c(1, j+1)]*Dm[,(j+2)])
    }

    repeat{
      D <- Dm
      for (i in 1:n){
        D[i, "tmt"] <- 1   #calculate criterion for unit i assigned to tmt 1

        if (!is.null(int)){
          D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]
        }

        plusD <- lossfunc(D, ...)

        D[i, "tmt"] <- -1 #calculate criterion for unit i assigned to tmt -1

        if (!is.null(int)){
          D[i, (j+3):(2*j+2)] <- -D[i, 2:(j+1)]
        }


        minusD <- lossfunc(D, ...)

        if (plusD < minusD){
          D[i,"tmt"] <- 1
          if (!is.null(int)){
            D[i, (j+3):(2*j+2)] <- D[i, 2:(j+1)]   }
        } else {
          D[i,"tmt"] <- -1
          if (!is.null(int)){
            D[i, (j+3):(2*j+2)] <- -D[i, 2:(j+1)]
          }
        }
      }

      if (lossfunc(D, ...) == lossfunc(Dm, ...)) break #if no change to resulting design matrix, terminate
      else Dm <- D

    }

    Dms[[m]] <- D




  }


  mindes <- which.min(unlist(lapply(Dms, lossfunc, ...)))  #find the optimum design matrix

  result <- as.data.frame(Dms[mindes][[1]])


  return(result)
}




#' Assuming a linear model for the response, allocate treatments randomly.
#'
#'
#' @param covar a dataframe for the covariates
#' @param lossfunc the objective function to minimize
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param init the number of units in the initial design, by default set to 1
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param ... further arguments to be passed to <lossfunc>

#'
#' @return design matrix and value of objective function at each sample size
#'
#'
#' @export
#'

linear.rand <- function(covar, lossfunc=calc.D, int=NULL, init=1, same.start=NULL,  ...){
  #purpose: use Atkinson approach to determine treatment assignment, allowing for any kind of optimality criterion,
  #         and also allowing for interactions
  #input: dataframe covar of covariate values
  #       lossfunc is the optimality criterion
  #       int is set to TRUE if you want to include the treatment-covariate interaction
  #output: design matrix


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  opt <- c()

  if (!is.null(same.start)){
    D <- same.start
  }else{
    D <- as.matrix(cbind(1, covar[1, ], sample(c(-1, 1), 1, replace=T))) # randomly select treatment for first unit
    if (!is.null(int)){
      D <- cbind( D, t(D[,c(2:(j+1))]*D[, (j+2)]))
    }

  }



  for (i in (init+1):n){


      new.tmt <- sample(c(-1,1), size=1)


    if (is.null(int)){
      new.d <- c(1, as.numeric(covar[i, ]), new.tmt)
      D <- as.matrix(rbind(D, new.d)) #add the new row to the design matrix
    } else {
      new.d <- c(1, as.numeric(covar[i, ]), new.tmt, t(as.numeric(covar[i, ])*new.tmt))
      D <- as.matrix(rbind(D, new.d ))
    }


    opt <- c(opt, lossfunc(D , ...))

  }
  row.names(D) <- NULL
  D <- as.data.frame(D)
  if (is.null(int)){
    colnames(D) <- c("intercept", rep("covar", j), "tmt")
  } else {
    colnames(D) <- c("intercept", rep("covar", j), "tmt",  rep("int", j))
  }



  results <- list(D=D, opt=opt)
  return(results)

}





#' Assuming a linear model for the response, allocate treatment sequentially based on an optimal design criterion
#'
#'
#' @param covar a dataframe for the covariates
#' @param lossfunc the objective function to minimize
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param init the number of units in the initial design, by default set to 1
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return design matrix and value of objective function at each sample size
#'
#'
#' @export
#'

atkins <- function(covar, lossfunc=calc.D, int=NULL, init=1, same.start=NULL, stoc=T, ...){


  n <- nrow(covar)
  j <- ncol(covar) #covar must be a dataframe
  opt <- c()
  all.probs <- c()
  loss.p <- loss.m <-c()

  if (!is.null(same.start)){
    D <- same.start
  }else{
    D <- as.matrix(cbind(1, covar[1, ], sample(c(-1, 1), 1, replace=T))) # randomly select treatment for first unit
  if (!is.null(int)){
     D <- cbind( D, t(D[,c(2:(j+1))]*D[, (j+2)]))
  }

  }



  for (i in (init+1):n){

    if (is.null(int)){
      #assuming no interactions
      #current row of design matrix, assuming tmt is 1
      d.plus <- lossfunc(rbind( D, c(1, as.numeric(covar[i, ]), 1)) , ...)
    } else {
      #assume tmt-covariate interactions
      #current row of design matrix, assuming tmt is -1
      d.plus <- lossfunc(rbind( D, c(1, as.numeric(covar[i, ]), 1, t(as.numeric(covar[i, ])*1 )))  , ...)
    }

    if (is.null(int)){
      d.minus <- lossfunc(rbind( D, c(1, as.numeric(covar[i, ]), -1)) , ...)
    }  else {
      d.minus <- lossfunc(rbind( D, c(1, as.numeric(covar[i, ]), -1, t(as.numeric(covar[i, ])*-1 ))) , ...)
    }


    probs <- (1/d.minus)/(1/d.plus + 1/d.minus)

    if (stoc==T){
      new.tmt <- sample(c(-1,1), size=1,  prob=c(probs, 1-probs))
    } else{

      if (d.plus > d.minus) {
        new.tmt <- -1
      } else if (d.plus < d.minus) {
        new.tmt <- 1
      } else if (d.plus == d.minus) {
        new.tmt <- sample(c(-1,1),size=1)
      }
    }

    loss.p <- c(loss.p, d.plus)
    loss.m <- c(loss.m, d.minus)

    if (is.null(int)){
      new.d <- c(1, as.numeric(covar[i, ]), new.tmt)
      D <- as.matrix(rbind(D, new.d)) #add the new row to the design matrix
    } else {
      new.d <- c(1, as.numeric(covar[i, ]), new.tmt, t(as.numeric(covar[i, ])*new.tmt))
      D <- as.matrix(rbind(D, new.d ))
    }


    opt <- c(opt, lossfunc(D , ...))

  }
  row.names(D) <- NULL
  D <- as.data.frame(D)
  if (is.null(int)){
    colnames(D) <- c("1", rep("covar", j), "tmt")
  } else {
    colnames(D) <- c("1", rep("covar", j), "tmt",  rep("int", j))
  }



  results <- list(D=D, opt=opt, loss.m=loss.m, loss.p=loss.p)
  return(results)

}
