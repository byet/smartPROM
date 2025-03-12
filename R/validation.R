#' Learns a BN model from data using either a bnlearn algorithm object (for both structure and parameters) or a dataframe of arcs (only parameters)
#'
#' @param train_data dataframe of BN variables. Each column is factor.
#' @param algorithm bnlearn algorithm object
#' @param arc_set bnlearn BN structure object
#' @param na_omit Logical
#'
#' @return gRain model
#' @importFrom bnlearn hc
#' @importFrom bnlearn arcs<-
#' @importFrom bnlearn bn.fit
#' @importFrom bnlearn as.grain
#' @export
#'
#' @examples
learn_bn <-
  function(train_data,
           algorithm = hc,
           arc_set = NULL,
           na_omit = TRUE) {
    if (na_omit) {
      train_data <- na.omit(train_data)
    }
    
    train_data <- apply(train_data, 2, as.factor)
    train_data <- as.data.frame(train_data, stringsAsFactors = TRUE)
    if (is.null(arc_set)) {
      bnfit_graph <- algorithm(train_data, restart=20, perturb=3)
    } else {
      bnfit_graph <- empty.graph(names(train_data))
      arcs(bnfit_graph) <- arc_set
    }
    bnfit_model <- bn.fit(bnfit_graph, train_data, method = "hard-em", fit = "bayes")
    as.grain(bnfit_model)
    
  }


#' Learn models for every row of data by leaving that row of data out. Used for leave one out validation.
#'
#' @param train_data dataframe of BN variables. Each column is factor.
#' @param algorithm bnlearn algorithm object
#' @param arc_set dataframe of 2 columns 'from' and 'to' to define BN arcs
#' @param na_omit Logical
#'
#' @return list of bn models learned from leave one out validation
#' @export
#'
#' @examples
learnOneOut <- function(train_data,
                        algorithm = hc,
                        arc_set = NULL,
                        na_omit = TRUE) {
  n <- nrow(train_data)
  bns <- vector(mode = "list", length = n)
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  
  for (i in 1:n) {
    setTxtProgressBar(pb, i)
    dataInstance <- train_data[i, ]
    index <- row.names(dataInstance)
    learnData <- train_data[-i, ]
    bn <- learn_bn(learnData, algorithm, arc_set, na_omit)
    bns[[i]] <- list(index = index, bn = bn)
  }
  bns
}


#' Naive BN arcs
#' Defines the arcs for naive BN in a dataframe.
#'
#' @param factors vector of factor variable names
#' @param questions vector of question variable names
#'
#' @return dataframe of two columns "from" "to" that describes the naive BN structure
#' @export
#'
#' @examples naivebn_arcs(c("L1","L2"),c("Q1","Q2","Q3"))
naivebn_arcs <- function(factors, questions) {
  bnarcs <-
    data.frame(matrix(vector(), 0, 2, dimnames = list(c(), c("from", "to"))), stringsAsFactors = F)
  for (factor in factors) {
    arcs_for_factor <-
      expand.grid(factor, questions, stringsAsFactors = FALSE)
    names(arcs_for_factor) <- c("from", "to")
    bnarcs <- rbind(bnarcs, arcs_for_factor)
  }
  bnarcs
}


#' Builds an arcset dataframe with a naive BN structure.
#'
#' @param vars array of variable names. Latent variable in the naive BN structure is the last element(s) of the array
#' @param nfac integer. The number of latent variables
#' @param factor_names optional array of variable names. Names of latent variables
#'
#' @return bnlearn BN structure object
#' @importFrom bnlearn empty.graph
#' @importFrom bnlearn arcs
#' @export
#'
#' @examples
#'
#'
#'
naivebn_str <- function(vars,
                            nfac = 1,
                            factor_names = NULL) {
  bnstr = empty.graph(vars)
  factors <- tail(vars, nfac)
  if (is.null(factor_names)) {
    questions <- vars[!vars %in% factors]
  } else {
    questions <- factor_names
  }
  bnarcs <- naivebn_arcs(factors, questions)
  arcs(bnstr) <- bnarcs
  bnstr
}

#' Distance of two marginal distributions from
#'
#' @param compared_dist gRain posterior marginal distribuiton object
#' @param reference_dist gRain posterior marginal distribuiton object
#' @param distance "kld" or "hellinger"
#'
#' @return vector of distances
#' @export
#'
#' @examples
marginals_distance <- function(compared_dist, reference_dist, distance =
                                 "kld") {
  if (distance == "hellinger") {
    distMeasure <- held
  } else {
    distMeasure <- kld
  }
  
  if (is.list(compared_dist) & is.list(reference_dist)) {
    mapply(distMeasure, reference_dist, compared_dist)
  } else {
    distMeasure(reference_dist, compared_dist)
  }
}

#' Computes the distance of target variables of each CAT iteration to those from last iteration
#'
#' @param cat_iter_results result object from cat_iterations
#' @param distance "kld" or "hellinger"
#'
#' @return dataframe of distances where rows are outcome variables columns are iterations
#' @export
#'
#' @examples
distance_to_all_queries <- function(cat_iter_result, distance = "kld") {
  final_posterior <- cat_iter_result$posterior
  iteration_posteriors <- cat_iter_result$iter_posteriors
  sapply(iteration_posteriors, function(x) {
    marginals_distance(x, final_posterior, distance)
  })
}

#' Calculates average CAT target variable distance results in cross validation
#'
#' @param cat_iter_results result objects of cat_iterations from multiple models
#' @param distance "kld" or "hellinger"
#'
#' @return dataframe of average distances where rows are outcome variables columns are iterations
#' @export
#'
#' @examples
average_dist_models <- function(cat_iter_results, distance = "kld") {
  distances <- lapply(cat_iter_results, function(cat_iter_result)
    distance_to_all_queries(cat_iter_result, distance))
  Reduce("+", distances) / length(distances)
}


#' Returns maximum a posteriori of nodes from BN
#'
#' @param bn gRain bn Model
#' @param nodes array of string. Node names
#' @param excludeObs logical. Exclude obeserved variables from target
#'
#' @return named vector
#' @importFrom gRain querygrain
#' @export
#'
#' @examples
pred <- function(bn, nodes, excludeObs = FALSE) {
  posterior <- querygrain(bn, nodes, exclude = excludeObs)
  map <- lapply(posterior, function(x) {
    names(which.max(x))
  })
  setNames(unlist(map, use.names = FALSE), names(map))
}


#' Compute accuracy and mean absolute error for questoins that are not answered for a patient and model
#'
#' @param data dataframe row. Columns are variables, values are states of variables.
#' @param bn gRain BN object
#' @param questions vector of strings. Names of variables that will be compared
#' @param numericStates logical. Are states numeric
#'
#' @return vector of accuracy and mean absolute error
#' @export
#'
#' @examples
questionErrors <- function(bn, data, questions, numericStates = TRUE) {
  evNodes <- bn$evidence$nodes
  unknownNodes <- questions[!questions %in% evNodes]
  map <- pred(bn, questions, excludeObs  = FALSE)
  
  if (numericStates){
    # Numeric comparison based on factor levels
    combined <- rbind(data[,questions], map[questions])
    map <- combined[nrow(combined), questions]
    mode(map) <- "numeric"
    mode(data) <- "numeric"
    trueVals <- data[,questions]
    predictions <- map[,questions]
  } else{
    trueVals <- data[questions]
    predictions <- map[questions]
  }
  # Ignore NA values
  trueVals_complete <- trueVals[!is.na(trueVals)]
  predictions_complete <- predictions[!is.na(trueVals)]
  
  classAccuracy <- mean(trueVals_complete == predictions_complete)
  if (numericStates) {
    absError <- mean(as.matrix(abs(trueVals_complete - predictions_complete)))
    # Percentage of answers 1 away
    dist1 <- mean(abs(trueVals_complete - predictions_complete) <= 1)
    # Percentage of answers 2 away
    dist2 <- mean(abs(trueVals_complete - predictions_complete) <= 2)
  } else{
    absError <- NA
  }
  setNames(c(classAccuracy, absError, dist1, dist2), c("Accuracy", "MAE", "1_away", "2_away"))
}


#' Compute error for CAT iterations
#'
#' @param bns list of gRain BN models for each iteration
#' @param data dataframe row. Columns are variables, values are states of variables.
#' @param questions vector of strings. Variable names that will be compared
#'
#' @return dataframe where rows are Accuracy and MAE, columns are questions
#' @export
#'
#' @examples
errorAllQuestions <- function(bns, data, questions, numericStates = TRUE) {
  sapply(bns, function(bn)
    questionErrors(bn, data, questions, numericStates = numericStates))
}



#' Posterior probabilities of all CAT iterations and all target variables for multiple cat_result objects
#'
#' @param results list. result objects of cat_iterations for multiple lines of data
#' @param targetVars vector of strings. Names of target variable
#' @param niter integer. Number of iterations to calcualte predictiosns for. If NULL it will be calculated for all iterations
#'
#' @return list of posterior predictions for each iteration from each model 
#' @export
#'
#' @examples
predictIterations <- function(results, targetVars, niter = NULL){
  if(is.null(niter)){
    niter <- length(results[[1]]$iter_posteriors)
  }
  iterPreds <- list()
  for(i in 1:niter){
    iterPreds[[i]] <- predictTargets(results, targetVars, i)
  }
  iterPreds
}

#' Posterior probabilities of all target variable in a CAT iteration for multiple cat_result objects
#'
#' @param results list. result objects of cat_iterations for multiple lines of data
#' @param targetVars vector of strings. Names of target variable
#' @param iteration integer. Cat iteration number predictions will be calculated
#'
#' @return list of dataframes of posterior probability distributions of target variables for an iteration from each model
#' @export
#'
#' @examples
predictTargets <- function(results, targetVars, iteration=1){
  targetsPreds <- list()
  for(target in targetVars){
    targetsPreds[[target]] <- probPredictions(results, target, iteration)
  }
  targetsPreds
}


#' Posterior probabilities for a target variable in a CAT iteration for multiple cat_result objects
#'
#' @param results list. result objects of cat_iterations for multiple lines of data
#' @param targetVar string. Name of target variable
#' @param iteration integer. Cat iteration number predictions will be calculated
#'
#' @return dataframe of posterior probability distributions of target variable for an iteration from each model
#' @importFrom dplyr bind_rows
#' @export
#'
#' @examples
probPredictions <- function(results, targetVar, iteration){
  probPredictions <- list()
  for(i in 1:length(results)){
    probPredictions[[i]] <- results[[i]]$iter_posteriors[[iteration]][[targetVar]]
  }
  as.data.frame(bind_rows(probPredictions))
}



#' AUCs of results
#'
#' @param results list. result objects of cat_iterations for multiple lines of data
#' @param dataset dataframe that contains the true values for the target variables. Number of rows of dataset should be the same a number of objects in the results
#' @param targetVars vector of strings. Names of target variables.
#'
#' @return dataframe of AUCs. Columns are target variables, rows are CAT iterations
#' @importFrom pROC multiclass.roc
#' @export
#'
#' @examples
performance <- function(results, dataset, targetVars, ci=TRUE, conf_level = 0.95, bootstrap_samples=2000){
  pp <- predictIterations(results, targetVars)
  iter_aucs <- data.frame()
  iter_lci <- data.frame()
  iter_uci <- data.frame()
  
  if(ci){
    print("Computing AUCs with bootstrapped confidence intervals")
  } else{
    print("Computing AUCs")
  }
  
  pb <- txtProgressBar(min = 0, max = length(pp), style = 3)
  for(i in 1:length(pp)){
    setTxtProgressBar(pb, i)
    aucs <- c()
    lower_ci <- c()
    upper_ci <- c()
    for(target in targetVars){
      iter_dataset <- cbind(y = dataset[,target], pp[[i]][[target]])
      auc <- multiclass.roc(iter_dataset$y,iter_dataset[, 2:ncol(iter_dataset), drop=FALSE])$auc[[1]]
      aucs <- c(aucs, auc)
      # Bootstrapped CI
      if(ci){
        resample <- bootstraps(iter_dataset, strata = "y", times = bootstrap_samples)
        resample_aucs <- map_dbl(
          resample$splits,
          function(x) {
            dat <- as.data.frame(x)
            roc <- multiclass.roc(dat$y, dat[, 2:ncol(dat), drop=FALSE])
            roc$auc[[1]]
          }
        )
        alpha <- 1 - conf_level
        lower_bound <- quantile(resample_aucs, probs = alpha / 2)
        upper_bound <- quantile(resample_aucs, probs = 1 - alpha / 2)
        
        lower_ci <- c(lower_ci, lower_bound)
        upper_ci <- c(upper_ci, upper_bound)
      }
    }
    iter_aucs <- rbind(iter_aucs, aucs)
    iter_lci <- rbind(iter_lci, lower_ci)
    iter_uci <- rbind(iter_uci, upper_ci)
  }
  colnames(iter_aucs) <- targetVars
  
  if(ci){
    colnames(iter_lci) <- targetVars
    colnames(iter_uci) <- targetVars
    return(list(aucs = iter_aucs, llci = iter_lci, ulci = iter_uci))
  }
  iter_aucs
}
