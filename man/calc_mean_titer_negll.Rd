% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_mean_titer_negll}
\alias{calc_mean_titer_negll}
\title{Calculate the total negative log-likelihood of a mean titer}
\usage{
calc_mean_titer_negll(predicted_mean, max_titers, min_titers, titer_sd)
}
\arguments{
\item{predicted_mean}{The predicted mean titer.}

\item{max_titers}{Numeric vector of the upper bounds of the measured titers.}

\item{min_titers}{Numeric vector of the lower bounds of the measured titers.}

\item{titer_sd}{The expected standard deviation of titers.}
}
\description{
This is a base function to sum the total _negative_ log likelihood of a mean titer.
}
\details{
This function calculates the total negative log-likelihood of a predicted mean 
titer given a set of titers. The main assumption is that both measurement error and 
titer variation are normally distributed. Note that the argument \code{titer_sd} is the 
total expected standard deviation of the titer set, i.e. measurement error plus titer 
variation.
}
