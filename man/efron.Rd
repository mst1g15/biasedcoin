% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minimization.R
\name{efron}
\alias{efron}
\title{Allocate treatment basead on Efron's biased coin}
\usage{
efron(covar, p = 2/3)
}
\arguments{
\item{covar}{a dataframe for the covariates}

\item{p}{the probability for selecting treatment 1 when score1 > score-1}
}
\value{
design matrix
}
\description{
Allocate treatment basead on Efron's biased coin
}
