

#' Randomly allocate treatment
#'
#' @param covar a dataframe for the covariates
#'
#' @return design matrix
#' @export
#'
rand <- function(covar){

  n <- nrow(covar)

  tmt <- sample(c(1, -1), n, replace=T) #sample randomly

  #design matrix
  D <- cbind(1, covar=covar, tmt=tmt)
  return(D)

}



#' Allocate treatment basead on Efron's biased coin
#'
#' @param covar a dataframe for the covariates
#' @param p the probability for selecting treatment 1 when score1 > score-1
#'
#' @return design matrix
#'
#' @export
#'
#'
efron <- function(covar, p=2/3){


  n <- nrow(covar)
  tmt <- c()

  for (i in 1:n){
    d <- sum(tmt==1)-sum(tmt==-1) #find score

    if (d>0){
      newtmt <- sample(c(-1, 1), 1, prob=c(p, 1-p))
    } else if (d<0){
      newtmt <- sample(c(-1, 1), 1, prob=c(1-p, p))
    }  else{
      newtmt <- sample(c(-1,1), 1)
    }

    tmt <- c(tmt, newtmt)
  }

  #design matrix
  D <- cbind(1, covar, tmt)

  return(D)
}



#' Allocate treatment based on an imbalance criterion - most general form of minimization
#'
#' @param covar a dataframe for one covariate ONLY (not applicable for multiple covariates)
#' @param imb is a measure of imbalance. It can be set to "efron", "ks" or "max.imb"
#' @param p takes value between 0.5 and 1 and determines randomness in treatment assignment
#'
#' @return design matrix
#'
#' @export
#'
#'
min.gen <- function(covar, imb, p=2/3){


  n <- nrow(covar) #covar must be a dataframe
  covar <- as.vector(covar)
  label <- 1:n
  tmt <- rep(NA, n)  #tmts not assigned
  dat <- as.data.frame(cbind(label=label, covar=covar, tmt=tmt))
  colnames(dat)[2] <- "covar"

  for (i in 1:n){

    current <- dat[1:i,]
    #suppose we assign tmt 1 to unit i
    current$tmt[i] <- 1
    #sort according to z
    T.1 <- current[order(current$covar),] #this is needed for Max.imb

    #find measure of imbalance
    if (imb=="max.imb") {
      #find the position of tmt of interest in sorted list
      k <- match(i, T.1$label)
      D.1 <- max.imb(T.1$tmt, k)  #imbalance measure proposed in Hu and Hu (2012)
    } else if (imb=="efron") {
      D.1 <- abs(sum(T.1$tmt==1)-sum(T.1$tmt==0))  #basically randomization, "biased coin"
    } else if (imb=="ks") {
      if (length(T.1$covar[T.1$tmt==1])==0 | length(T.1$covar[T.1$tmt==-1])==0) { D.1 <- 1
      } else {D.1 <- ks.test(T.1$covar[T.1$tmt==1], T.1$covar[T.1$tmt==-1])$statistic}  #kolmogorov-smirnov distance
    }

    #repeat supposing we assign tmt 0 to unit i
    current$tmt[i] <- -1
    T.0 <- current[order(current$covar),]

    if (imb=="max.imb") {
      k <- match(i, T.0$label)
      D.0 <- max.imb(T.0$tmt, k)
    } else if (imb=="efron") {
      D.0 <- abs(sum(T.0$tmt==1)-sum(T.0$tmt==0))
    } else if (imb=="ks") {
      if (length(T.1$covar[T.1$tmt==1])==0 | length(T.1$covar[T.1$tmt==-1])==0) { D.0 <- 1
      } else {D.0 <- ks.test(T.0$covar[T.0$tmt==1], T.0$covar[T.0$tmt==-1])$statistic}
    }


    #assign
    if (D.1 < D.0){
      dat$tmt[i] <- rbinom(1, 1, p)
    } else if (D.1 > D.0) {
      dat$tmt[i] <- rbinom(1, 1, (1-p))
    } else if (D.1 == D.0) {
      dat$tmt[i] <- rbinom(1, 1, 1/2)
    }
  }

  dat$tmt[dat$tmt==0]<- -1 #we want to code the treatments 1 and -1, not 1 and 0
  dat$label <- rep(1, n)  #intercept
  D <- as.data.frame(dat)  #design matrix
  colnames(D)<- c("1", "covar", "tmt")
  return(D)

}







#' Allocate treatment based on Hu and Hu (2012)'s version of minimization.
#' Hu and Hu (2012) imbalance measure (only applicable for one continuous covariate)
#' Compute imbalance measure proposed in Hu and Hu (2012):
#' sort list of currently assigned treatments according to the covariate
#' consider all subsets containing the current patient k
#' Find difference in number of patients who have been given treatment 1 vs treatment 0
#' Find the maximum difference across all subsets
#'
#' @param T01 vector of assigned treatments, sorted according to the covariate
#' @param k position of the current patient in the sorted list T01
#'
#' @return max.imb
#'
#' @export
#'

max.imb <- function(T01, k){

  n.T01 <- length(T01)
  Diff <- c()

  #we only want to consider intervals that include k
  for (i in 1:k){
    for (j in k:n.T01){
      diff <- sum(T01[i:j]==1) - sum(T01[i:j]==0)  # difference in patient numbers
      Diff <- c(Diff, as.numeric(abs(diff)))
    }
  }

  sup.d <- max(Diff)

  return(sup.d)

}



#' Allocate treatment using Senn's version of Minimization, appropriate for binary covariates only.
#  Aim is to "minimize the sum of the absolute values of the first row of XtX".
#'
#' @param covar a dataframe for covariates  (must be binary)
#' @param p takes value between 0.5 and 1 and determines randomness in treatment assignment
#'
#' @return design matrix
#'
#' @export
#'
#'
min.sen <- function(covar, p=2/3){

  #find the number of patients, n
  n <- nrow(covar) #covar must be a dataframe

  covar <- cbind(rep(1,n), covar)
  tmt <-c() # treatment values stored here
  q.list <- c() #store difference in q between tmt 1 and tmt -1


  for (i in 1:n){

    covi <-as.matrix(covar[1:i,])  #pick the first i units
    t.plus <- as.matrix(c(tmt, 1)) #assume current treatment is 1
    q.plus <- sum(abs(t(t.plus) %*% covi)) # calculate TtZ

    t.minus <- c(tmt, -1) #assume current treatment is -1
    q.minus <- sum(abs(t(t.minus) %*% covi)) #calculate TtZ
    q <- q.plus-q.minus #find difference

    #assign treatment - remember we want to minimize q
    if (q>0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(p, 1-p))  #if q>0, then qplus>qminus so qminus leads to minimum value
    } else if (q<0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(1-p, p))
    } else {
      newtmt <-  sample(c(-1, 1), 1) #random
    }

    tmt <- c(tmt, newtmt)
    q.list <- c(q.list, q) #store the qs

  }

  D <- cbind(covar, tmt)  #Design matrix

  return(list(q=q.list,D=D))

}


#' Classic version of Minimization, appropriate for binary covariates only
#'
#' @param covar a dataframe for covariates  (must be binary)
#' @param p takes value between 0.5 and 1 and determines randomness in treatment assignment
#'
#' @return design matrix
#'
#' @export
#'
#'
min.classic <- function(covar, p=2/3){
  # purpose: use minimization to allocate treatments where we aim to
  #          minimize the classic form of imbalance
  # input: a dataframe covar with binary covariate values
  # output: design matrix

  #find the number of patients and the number of covariates

    n <- nrow(covar)    #covar must be a dataframe
    k <- ncol(covar)

    tmt <- sample(c(-1,1), 1) #randomly allocate first treatment
    q.list <- c(0)  #first "difference" measure is 0 since we make a random dicision

  for (i in 2:n){
    comp<- c()
    for (j in 1:k){
         #the j-loop creates an indicator matrix for whether the previous units had the same covariate level
         #as the current unit, for each covariate
        comp <- cbind(comp, covar[1:(i-1),j]==rep(covar[i,j], (i-1)))

    }

    n.plus <- sum(tmt*comp==1)  #the number of units who have the same covariate level and have tmt 1
    n.minus <- sum(tmt*comp==-1) #the number of units who have the same covariate level and have tmt -1
    q <- n.plus-n.minus #difference in the number
    q.list <- c(q.list, q) #store the differences

    #assign treatments ... remember that we want to minimize q.
    if (q>0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(p, 1-p))  #if q>0, then qplus>qminus so qminus leads to minimum value
    } else if (q<0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(1-p, p))
    } else {
      newtmt <-  sample(c(-1, 1), 1) #random
    }

    tmt <- c(tmt, newtmt)

  }

  D <- cbind(rep(1,n), covar=covar, tmt=tmt) #design matrix

  return(list(q=q.list, D=D))

}





#' Allocate treatment using a version of minimization, appropriate for continuous covariates
#'
#' @param covar a dataframe for covariates  (must be continuous)
#' @param p takes value between 0.5 and 1 and determines randomness in treatment assignment
#'
#' @return design matrix
#'
#' @export
#'
min.cont <- function(covar, p=2/3){
  # purpose: use minimization to allocate treatments where covariates are continuous
  # input: dataframe Z covar with continuous covariate values
  # output: design matraix


  n <- nrow(covar) #covar must be a dataframe


  covar <- cbind(rep(1,n), covar) #add the intercept


  tmt <- sample(c(-1, 1), 1) #first treatment randomly allocated
  q.list <- c(0) #q is zero for first tmt since random allocation

  for (i in 2:n){

    covi <- as.matrix(apply(covar[1:i,], 2, discrete)) # send first i units to discretize function
    #we are dynamically discretizing by splitting at the median for each pass through this loop

    tmt.plus <- c(tmt, 1)
    q.plus <- sum(abs(t(tmt.plus) %*% covi))

    tmt.minus <- c(tmt, -1)
    q.minus <- sum(abs(t(tmt.minus) %*% covi))

    q <- q.plus-q.minus

    #assign treatments ... remember that we want to minimize q.
    if (q>0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(p, 1-p))  #if q>0, then qplus>qminus so qminus leads to minimum value
    } else if (q<0) {
      newtmt <- sample(c(-1, 1), 1, prob=c(1-p, p))
    } else {
      newtmt <-  sample(c(-1, 1), 1) #random
    }

    tmt <- c(tmt, newtmt)

    q.list <- c(q.list, q)

  }

  D <- cbind( covar, tmt) #design matrix

  return(list(q=q.list,D=D))

}







#' discretize a vector of continuous values into binary (0 and 1) by using the
#  median as a cut-off point
#' @param z a vector of continuous values
#'
#' @return a vector of descretized values
#'
#' @export
#'

discrete <- function(z){
  #purpose: discretize a vector of continuous values into binary (0 and 1) by using the
  # median as a cut-off point
  #input: vector z of continuous values
  #output: vector z discretized

  med <- median(z)
  n <- length(z)

  for (i in 1:n){
    if (z[i]< med) {
      z[i] <- 0
    } else if (z[i] > med) {
      z[i] <- 1
    } else {
      z[i] <- sample(c(0, 1), 1)
    }

  }

  return(z)
}


