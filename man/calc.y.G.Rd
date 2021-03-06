% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic.R
\name{calc.y.G}
\alias{calc.y.G}
\title{Given an information matrix and assuming logistic regression, compute the G-optimal objective function}
\usage{
calc.y.G(I, epsilon = 1e-05)
}
\arguments{
\item{I}{the information matrix}

\item{epsilon}{a small real number used for regularization. If set to zero,
no regularization takes place}
}
\value{
value of the G-optimal objective function
}
\description{
Given an information matrix and assuming logistic regression, compute the G-optimal objective function
}
