
#' Calculates expected optimality for a given trajectory after allocating treatments to the trajectory with
#' the coordinate exchange algorithm
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param n number of total patients
#' @param k number of "outer loops" in the coordinate exchange algorithm
#' @param sim number of trajectories to simulate
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param z.probs vector of probabilities for each level of covariate z
#' @param code set to NULL if (-1,1) coding is used for the treatments. Set to 0 if (0, 1) is used.
#' @param lossfunc the objective function to minimize
#' @param ... further arguments to be passed to <lossfunc>


#' @return loss of the design matrix which includes trajectory
#'
#'
#' @export
future.coordex <- function(D.fix, n.r, n,  k, sim, int, z.probs,  code=NULL, lossfunc, ...){


  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)

    losses <- rep(0, sim)

  for (q in 1:sim){
    covar <- gencov(z.probs,n.r, code)
    j <- ncol(covar)

      Dms <- list() #list to store k potential designs

      for (m in 1:k){

        if (!is.null(code)){
          D.r <- as.matrix(cbind(rep(1, n.r), covar, tmt=sample(c(0,1), n.r, replace=T))) #design matrix with random treatment assignment

        }else{
                  D.r <- as.matrix(cbind(rep(1, n.r), covar, tmt=sample(c(-1,1), n.r, replace=T))) #design matrix with random treatment assignment

        }
        if (!is.null(int)){

          if (dim(D.r)[1]==1){
            D.r <- cbind(D.r, t(D.r[,-c(1, j+2)]*D.r[,(j+2)]))
          }else{
                        D.r <- cbind(D.r, D.r[,-c(1, j+2)]*D.r[,(j+2)])

          }


        }
        #print(D.fix)
        #print(D.r)
        colnames(D.r) <-  colnames(D.fix)
        Dm <- rbind(D.fix, D.r)

                #count <- 1

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

              oldD.l <- lossfunc(D, ...)
              newD.l <- lossfunc(D.new, ...)




            if (oldD.l <=  newD.l){
             D <- D
            } else {
              D <- D.new
            }

            }

          l1 <- lossfunc(D, ...)
          l2 <- lossfunc(Dm,  ...)

          cat(q, i, l1, l2, "\n")
          if (l1==l2) {
            break
          }else
          Dm <- D

        }

        Dms[[m]] <- D




      }

      all.losses <- unlist(lapply(Dms, lossfunc,  ... ))
      mindes <- which.min(all.losses)  #find the optimum design matrix

      losses[q] <- all.losses[mindes]
      all.losses[mindes]
      print(Dms[mindes])
  }

  mean.loss <- mean(losses)


  return(mean.loss)



}









#' Calculate expected optimality for a given trajectory using sequential approach.
#' We assume a linear model for the response.
#' @param D.fix Design matrix constructed using the true covariates in the experiment so far
#' @param n.r length of trajectory to be simulated
#' @param sim number of trajectories to simulate
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param z.probs vector of probabilities for each level of covariate z
#' @param lossfunc the objective function to minimize
#' @param ... further arguments to be passed to <lossfunc>

#' @return loss of the design matrix which includes trajectory
#'
#'
#' @export


future <- function(D.fix, n.r,   sim, int, z.probs,lossfunc,  ...){


  n.fix <- nrow(D.fix)
  D.c <- ncol(D.fix)
  losses <- rep(0, sim)


  for (q in 1:sim){

    covar <- gencov(z.probs,n.r)


    D <- D.fix

    for (i in 1:n.r){

      if (!is.null(int)){
        mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1,as.numeric(covar[i, ]) ))
      } else{
        mat.plus <- rbind(D, c(1, as.numeric(covar[i, ]), 1))
      }

      d.plus <-  lossfunc(mat.plus,  ...)                #Optimality criterion

      #design matrix if tmt -1 is assigned
      if (!is.null(int)){
        mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1, -as.numeric(covar[i, ]) ))
      } else{
        mat.minus <- rbind(D, c(1, as.numeric(covar[i, ]), -1)) #Design matrix with new treatment =-1
      }


      d.minus <- lossfunc(mat.minus,  ...)                 #Optimality criterion

      new.tmt <- ifelse(d.plus < d.minus, 1, -1)         #Assign treatments

      #new row of design matrix
      if (!is.null(int)){
        new.d <- as.numeric(c(1, covar[i, ], new.tmt, covar[i, ]*new.tmt))
      } else {
        new.d <- as.numeric(c(1, covar[i, ], new.tmt))
      }


      D <- as.matrix(rbind(D, as.numeric(new.d)))                #Add the new row to the design matrix

      losses[q] <- lossfunc(D,  ...)

    }


  }

  mean.loss <- mean(losses)


  return(mean.loss)



}









#' Allocate treatments according to an information matrix based optimality criterion allowing for a pseudo-nonmyopic approach.
#' We assume a linear model for the response.
#' @param covar a dataframe for the covariates
#' @param sim number of trajectories to simulate
#' @param int set to T if you allow for treatment-covariate interactions in the model, NULL otherwise
#' @param z.probs vector of probabilities for each level of covariate z
#' @param lossfunc the objective function to minimize
#' @param stoc set to T if treatments are allocated using a stochastic method where the probability is
#' determined by the optimality crtierion. Set to F if treatments are allocated deterministically.
#' @param same.start the design matrix to be used for the initial design. If set to NULL, function generates initial design.
#' @param coordex set to T if coordinate exchange algorithm is used to allocate treatments in the trajectory, set to NULL for sequential approach.
#' @param ... further arguments to be passed to <lossfunc>
#'
#' @return Design matrix D, value of objective function
#'
#'
#' @export

simfuture <- function(covar, sim,  int=NULL,  z.probs, lossfunc, stoc=T, same.start=NULL, coordex=NULL, ... ){

  n <- nrow(covar)
  j <- ncol(covar)
  names(covar) <- NULL

  all.loss.p <- all.loss.m <-c()

  if (!is.null(same.start)){
    design <- same.start
    init <- nrow(same.start)+1
  }else{
    design <- c()
    init <- 1
  }
  opt <-c()



  for (i in init:n){


    if (!is.null(int)){
      design.p <- rbind(design, c(1, as.numeric(covar[i,]), 1, as.numeric(covar[i,])))
      design.m <- rbind(design, c(1, as.numeric(covar[i,]), -1, - as.numeric(covar[i,])))

    }else{
        design.p <- rbind(design, c(1, as.numeric(covar[i,]), 1))
        design.m <- rbind(design, c(1, as.numeric(covar[i,]), -1))
    }



    if (i!=n){

        if(!is.null(ncol(z.probs))){
          z.probs.i <- as.data.frame(z.probs[ (i+1):n,])
        }else{
          z.probs.i <- z.probs
        }


      if(!is.null(coordex)){

        loss.p <- future.coordex(design.p, n-i, n, k, sim, int, z.probs.i, code=NULL, lossfunc, ...)
        loss.m <- future.coordex(design.m, n-i, n, k, sim, int, z.probs.i, code=NULL, lossfunc, ...)

        }else{

          loss.p <- future(design.p, n-i,  sim, int, z.probs.i, lossfunc, ...)
          loss.m <- future(design.m, n-i,  sim, int, z.probs.i, lossfunc, ...)

       }

    }else {
      loss.p <- lossfunc(design.p, ...)
      loss.m <- lossfunc(design.m, ...)
    }


    all.loss.p <- c(all.loss.p, loss.p)
    all.loss.m <- c(all.loss.m, loss.m)


        probs <- (1/loss.p)/(1/loss.p+1/loss.m)

    if (stoc==T){
          new.tmt <- sample(c(-1,1), 1, prob=c(1-probs, probs))

    }else{

      if (loss.p > loss.m) {
        new.tmt <- -1
      } else if (loss.p < loss.m) {
        new.tmt <- 1
      } else if (loss.p == loss.m) {
        new.tmt <- sample(c(-1,1), 1)
      }

    }

    if(new.tmt == 1){
      design <- design.p
    }else{
      design <- design.m
    }

    opt <- c(opt, lossfunc(design, ...))
  }



  return(list(design=design, opt=opt, loss.m=all.loss.m, loss.p=all.loss.p))

}




