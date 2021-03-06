% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic_nonmyopic.R
\name{exp.loss}
\alias{exp.loss}
\title{Break down expected future optimality into two components:
1) assuming that the current response is 0
2) assuming that the current response is 1
Find the weighted average of the two cases}
\usage{
\method{exp}{loss}(z.now, t.now, zp, N, design, int, lossfunc, beta, y, bayes,
  dyn = NULL, ...)
}
\arguments{
\item{z.now}{vector of covariate values for current unit}

\item{t.now}{treatment of current unit}

\item{zp}{vector of probabilities for each level of covariate z (needs to in the same order as all.z below)}

\item{N}{natural number greater than 0 for horizon}

\item{design}{design matrix constructed for all units up until the current unit}

\item{int}{set to NULL if there are no interactions, set to T of there are interactions}

\item{lossfunc}{the objective function to minimize}

\item{beta}{estimate of the regression coefficients}

\item{y}{responses that have been observed up until the current unit}

\item{bayes}{set to T if bayesglm is used instead of glm. Default prior assumed.}

\item{dyn}{set to T if there is a dynamic covariate}

\item{...}{further arguments to be passed to <lossfunc>}
}
\value{
expected value of objective function one step ahead in the future
}
\description{
Break down expected future optimality into two components:
1) assuming that the current response is 0
2) assuming that the current response is 1
Find the weighted average of the two cases
}
