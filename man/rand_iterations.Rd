% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{rand_iterations}
\alias{rand_iterations}
\title{Run random iterations on BN}
\usage{
rand_iterations(
  bn,
  target,
  query,
  data,
  evidence = NULL,
  num_iter = 3,
  store_bns = FALSE
)
}
\arguments{
\item{bn}{gRain BN object}

\item{target}{array of strings. Names of target variables}

\item{query}{array of strings. Names of observed variables which the data will be entered.}

\item{data}{dataframe row. Columns are variables, values are states of variables.}

\item{evidence}{list of keys: variables, values: variable states. Initial evidence for the BN.}

\item{num_iter}{integer. Maximum number of iterations to run CAT}

\item{store_bns}{logical. Stores BNs at each iterations}
}
\description{
Iteratively calculates the BN with randomly selected variables
}
