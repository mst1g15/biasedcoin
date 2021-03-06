

#' Create a set of binary covariates with a specific covariate distribution
#'
#' @param z.probs if a scalar, assume there is one binary covariate z and
#' the scalar is P(z=1).
#' if z.probs is a vector of length k, there are k static covariates, z1, ..., zk, and
#' the the jth element of the vector is  P(zj=1).
#' if z.probs is a matrix with k columns and n.r rows, there are k dynamic covariates z1, ...zk.
#' the (i,j)th element is the probability that the jt unit of the ith covariate is equal to 1
#' @param n.r number of covariate values to generate
#' @param code set to NULL if (-1,1) coding is used for the treatments. Set to 0 if (0, 1) is used.

#' @return dataframe with covariates
#'
#' @export

gencov <- function(z.probs, n.r, code=NULL){

  if(!is.null(dim(z.probs))){  #there are dynamic covariates
    cov <- matrix(0, nrow=nrow(z.probs), ncol=ncol(z.probs))
    for (j in 1:ncol(z.probs)){
      for (i in 1:nrow(z.probs)){

        if (!is.null(code)){
          cov[i,j] <- sample(c(1,0),  1,  prob=c(z.probs[i,j], 1-z.probs[i,j]))

        }else{
          cov[i,j] <- sample(c(1,-1),  1,  prob=c(z.probs[i,j], 1-z.probs[i,j]))

        }



      }
    }
    covar <- as.data.frame(cov)

  }else if(length(z.probs)>1){  #there are multiple covariates but they are not dynamic
    #multiple covariates
    cov <- c()
    for (i in 1:length(z.probs)){

      if (!is.null(code)){
        cov <- cbind(cov,sample(c(1,0), n.r, replace=T, prob=c(z.probs[i], 1-z.probs[i])))

      }else{
        cov <- cbind(cov,sample(c(1,-1), n.r, replace=T, prob=c(z.probs[i], 1-z.probs[i])))

      }


    }

    covar <- as.data.frame(cov)
  }else if(length(z.probs)==1){  #there is one static covariate

    if (!is.null(code)){
      covar <- as.data.frame(sample(c(1,0), n.r, replace=T, prob=c(z.probs, 1-z.probs)))

    }else{
      covar <- as.data.frame(sample(c(1,-1), n.r, replace=T, prob=c(z.probs, 1-z.probs)))

    }

  }


  return (covar)
}








#' Given z.probs which takes a specific form, convert it to an nxp dataframe where the (i,j)th entry is the probability
#' that the ith patient's jth covariate value is 1
#'
#' @param z.probs if a scalar, assume there is one binary covariate z and
#' the scalar is P(z=1).
#' if z.probs is a vector of length k, there are k static covariates, z1, ..., zk, and
#' the the jth element of the vector is  P(zj=1).
#' if z.probs is a matrix with k columns and n.r rows, there are k dynamic covariates z1, ...zk.
#' the (i,j)th element is the probability that the jt unit of the ith covariate is equal to 1
#'
#' @return nxp dataframe where the (i,j)th entry is the probability
#' that the ith patient's jth covariate value is 1
#'
#' @export


zpfunc <- function(z.probs){

  if(!is.null(dim(z.probs))){  #there are dynamic covariates
   if (ncol(z.probs)==1){
     zp <- cbind(1-z.probs, z.probs)

   }else{
      a <- cbind(1-z.probs[,1], z.probs[,1])
    b <- cbind(1-z.probs[,2], z.probs[,2])
    q <- nrow(z.probs)
    zp <- c()
    for (i in 1:q){
      zp <- rbind(zp, c(outer(a[i,], b[i,])))
    }

   }


  }else if(length(z.probs)>1){  #there are multiple covariates but they are not dynamic
    #multiple covariates
    a <- c(1-z.probs[1], z.probs[1])
    b <- c(1-z.probs[2], z.probs[2])

    zp <- c(outer(a,b))

  }else if(length(z.probs)==1){  #there is one static covariate

    zp <- c(1-z.probs, z.probs)


  }


  return (zp)

}



if(FALSE){




gencov.cont <- function(covdist, para1, para2 , n){
  #purpose: create a set of covariates with a specific covariate distribution
  #inputs: z.probs: if a scalar, P(z=1). If a vector, c(P(z1=1).=, P(z2=1), P(z3=1), ..., P(zk=1))
  #        n.r: the number of patients



  if(length(para1)>1){  #dynamic covariates
    cov <- c()
    for (i in 1:length(para1)){
      cov <- c(cov,covdist(1, para1[i], para2[i]))
    }

    covar <- as.data.frame(cov)
  }else if(length(para1)==1){  #there is one static covariate
    covar <- as.data.frame(covdist(n, para1, para2))
  }


  return (covar)
}

}

