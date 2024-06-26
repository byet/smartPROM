% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{performance}
\alias{performance}
\title{AUCs of results}
\usage{
performance(results, dataset, targetVars)
}
\arguments{
\item{results}{list. result objects of cat_iterations for multiple lines of data}

\item{dataset}{dataframe that contains the true values for the target variables. Number of rows of dataset should be the same a number of objects in the results}

\item{targetVars}{vector of strings. Names of target variables.}
}
\value{
dataframe of AUCs. Columns are target variables, rows are CAT iterations
}
\description{
AUCs of results
}
