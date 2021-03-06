% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic.R
\name{Dopt.y.t.init}
\alias{Dopt.y.t.init}
\title{Given an initial design and assuming logistic regression, compute the D-optimal objective function
(to be used in conjunction with optim)}
\usage{
Dopt.y.t.init(t, z, int, beta, epsilon = 1e-05)
}
\arguments{
\item{t}{vector of treatments}

\item{z}{vector of covariate values}

\item{int}{set to TRUE if treatment-covariate interactions are included in the model}

\item{beta}{current estimate of regression parameters}

\item{epsilon}{a small real number used for regularization. If set to zero,
no regularization takes place}
}
\value{
value of the D-optimal objective function
}
\description{
Given an initial design and assuming logistic regression, compute the D-optimal objective function
(to be used in conjunction with optim)
}
