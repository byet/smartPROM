#library(dplyr)
#library(ROI)
#library(ROI.plugin.glpk)
#library(ompr)
#library(ompr.roi)


#' Discretize numeric vector
#' Discretize a numeric vector x with cutoff points.
#' 
#' @param x vector of numbers. Numeric vector
#' @param cutoff vector of numbers. Cutoff points.
#' @param states vector of strings. State labels. The number of state labels should be 1 more than number of cutoff points
#'
#' @return vector of factors
#' @export
#'
#' @examples discretize(c(1,1,1,2,3,4,5,6,7),c(3.5,6.5),c("Low","Moderate","High"))
discretize <- function(x, cutoff, states = NULL) {
  min_d <- min(x)
  max_d <- max(x)
  
  if (is.null(states)) {
    states <- list()
    for (a in 1:(length(cutoff) + 1)) {
      states[[a]] <- paste("State", a, sep = "_")
      
    }
    states <- unlist(states)
  }
  
  cut(
    x,
    breaks = c(min_d, cutoff, max_d),
    include.lowest = TRUE,
    labels = states
  )
  
}


# Match Normal Distributions
# Finds the point where two normal distributions have equal densities when weights are equal


#' Match Normal Distributions
#' Finds the point where two normal distributions have equal densities when weights are equal. Weights can be assigned.
#' 
#'
#' @param weights vector of two probabilities. Weights of distributions
#' @param mu vector of two numbers. Means of normal distributions
#' @param sd vector of two numbers. SDs of normal distributions
#'
#' @return number
#' @export
#'
#' @examples match_normal(c(0.5,0.5), c(3,2), c(1,1.5))
match_normal <- function(weights, mu, sd){
  require(nleqslv)
  dn <- function(x) {
    
    f <- weights[1]*dnorm(x, mu[1], sd[1]) - weights[2]*dnorm(x, mu[2], sd[2])
    
  }
  
  x <- (max(mu) - min(mu)) / 2 + min(mu)
  res <- nleqslv::nleqslv(x,dn)
  res$x
}





#' MIPM intervals
#'
#' @param x vector of numbers. Numeric vector
#' @param threshold vector of numbers. Threshold values to determine the threshold values.
#' @param MCID number. MCID number as a lower bound for the intervals
#' @param n_interval number of intervals
#' @param states state names
#' @param Big_M Big M value for the MILP
#' @param Lower_B lower bound for the intervals
#' @param Upper_B upper bound for the intervals
#'
#' @return vector of discretized factors
#' @importFrom ompr add_variable
#' @importFrom ompr set_objective
#' @importFrom ompr add_constraint
#' @importFrom ompr solve_model
#' @importFrom ompr MIPModel
#' @importFrom ompr.roi with_ROI
#' @export
#'
#' @examples
mipm_intervals <- function(x, threshold, MCID, n_interval, states = NULL, Big_M = 200, Lower_B, Upper_B){
  
  result <- MIPModel() |>
    add_variable(d[i], i = 1:n_interval, type = "continuous", lb = Lower_B, ub = Upper_B) |>
    add_variable(y[i], i = 1:n_interval, type = "binary") |>
    add_variable(b[i, j], i = 1:n_interval, j = 1:n_interval, type = "binary") |>
    add_variable(FU[i, j], i = 1:n_interval, j = 1:n_interval, type = "continuous") |>
    add_variable(FL[i, j], i = 1:n_interval, j = 1:n_interval, type = "continuous") |>
    set_objective(sum_expr(FU[i, j] + FL[i, j], i = 1:n_interval, j = 1:n_interval), "min") |>
    add_constraint(sum_expr(d[i], i = 1:n_interval) == Upper_B) |>
    add_constraint(d[i] >= MCID, i = 1:n_interval) |>
    add_constraint(sum_expr(d[l], l = 1:i) >= (threshold - ((1-y[i])*Big_M)), i = 1:n_interval) |>
    add_constraint(sum_expr(d[l], l = 1:i) <= (threshold + ((1-y[i])*Big_M)), i = 1:n_interval) |>
    add_constraint(sum_expr(y[i], i = 1:n_interval) == 1) |>
    add_constraint(d[i] - d[j] == FU[i, j] - FL[i, j], i = 1:n_interval, j = 1:n_interval) |>
    add_constraint(FU[i, j] >= 0, i = 1:n_interval, j = 1:n_interval) |>
    add_constraint((b[i,j]*Big_M) >= FU[i, j], i = 1:n_interval, j = 1:n_interval) |>
    add_constraint(FL[i, j] >= 0, i = 1:n_interval, j = 1:n_interval) |>
    add_constraint(((1-b[i,j])*Big_M) >= FL[i, j], i = 1:n_interval, j = 1:n_interval) |>
    solve_model(with_ROI(solver = "glpk"))
  
  names <- c()
  for (i in 1:n_interval) {
    names[i] <- paste("d","[",i,"]", sep = "")
  }
  
  sol <- result$solution[names]
  intervals <- cumsum(as.vector(sol))
  intervals <- intervals[-length(intervals)]
  discretize(x, cutoff = intervals, states = states)
}
