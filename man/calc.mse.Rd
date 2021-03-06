% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mse.R
\name{calc.mse}
\alias{calc.mse}
\title{Given a design, the current estimate of beta (or true value of beta), compute the MSE matrix}
\usage{
calc.mse(beta, D, y, sim, epsilon = 1e-05, true.bvcov = NULL)
}
\arguments{
\item{beta}{the current estimates of the parameter values}

\item{D}{current design matrix}

\item{y}{vector of responses}

\item{sim}{number of simulated betas to generate}

\item{epsilon}{a small real number used for regularization. If set to zero,
no regularization takes place}

\item{true.bvcov}{set to the true values of beta if the var-covar matrix is to be computed using the true values}
}
\value{
MSE matrix
}
\description{
Given a design, the current estimate of beta (or true value of beta), compute the MSE matrix
}
