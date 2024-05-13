#' Calculate BN with randomly selected variables
#' Calculates the posteriors for a set of target variables given a random subset of observed variables
#' This function is used for benchmarking CAT with randomly selected inputs.
#'
#' @param bn gRain BN object
#' @param target array of strings. Names of target variables
#' @param query array of strings. Names of observed variables which the data will be entered.
#' @param data dataframe row. Columns are variables, values are states of variables.
#' @param evidence dataframe row. Columns are variables, values are states of variables.
#' @param numquery number of input variables that will be entered. If this value is smaller than the number of variables in query, a random subset will be selected.
#' @param allquery logical. Run by entering all query variables.
#'
#' @return list containing posterior = gRain posterior object, evidence = evidence data entered
#' @importFrom gRain querygrain
#' @importFrom gRain setEvidence
#' @export
#'
#' @examples
rand_calc <-
  function(bn,
           target,
           query,
           data,
           evidence = NULL,
           numquery = 3,
           allquery = FALSE) {
    # This function randomly selects numvar number of variables from query variables and runs the BN for a single row of data
    # When all is true, this functoin runs the BN by entering all query variables
    # 1- Enter Initial Evidence
    if (!is.null(evidence)) {
      bn <- setEvidence(bn, evidence = as.list(evidence))
      
      # !!!Remove any evidence variable from query if they are in query
    }
    
    if (!allquery) {
      if (length(query) < numquery) {
        numquery <- length(query)
        print("Warning number of query variables is greater than numquery.")
        print(paste("The model will be run for:", length(query), "variables."))
      }
      
      # Select a random sample
      entered_ev <- sample(query, numquery)
      # print(entered_ev)
      query_ev <- data[entered_ev]
    } else  {
      entered_ev <- query
      query_ev <- data[query]
    }
    
    bn <-
      setEvidence(bn, nodes = names(query_ev), states = as.character(unlist(query_ev)))
    list(posterior = querygrain(bn, nodes = c(target), type = "marginal"),
         evidence = data[entered_ev])
    
    
  }




#' Runs CAT iterations on BN
#'
#' Iteratively calculates the BN with most informative variables
#'
#' @param bn gRain BN object
#' @param target array of strings. Names of target variables
#' @param query array of strings. Names of observed variables which the data will be entered.
#' @param data dataframe row. Columns are variables, values are states of variables.
#' @param evidence list of keys: variables, values: variable states. Initial evidence for the BN.
#' @param num_iter integer. Maximum number of iterations to run CAT
#' @param entropy_check_thres number. Stopping threshold for target variable entropy
#' @param mi_check_thres number. Stopping threshold for mutual information of query variables.
#'
#' @return list containing posterior = gRain posterior object, evidence = evidence data entered
#' @importFrom gRain setEvidence
#' @importFrom gRain querygrain
#' @export
#'
#' @examples
cat_iterations <-
  function(bn,
           target,
           query,
           data,
           evidence = NULL,
           num_iter = 3,
           entropy_check_thres = NULL,
           mi_check_thres = NULL,
           store_bns = FALSE) {
    # This function iteratively runs a single data.frame row by iteratively asking the vairables with maximum mutual information
    init_query <- query
    init_bn <- bn
    entered_ev <- c()
    bns <- vector(mode = "list", length = num_iter)
    posterior_targets <- vector(mode = "list", length = num_iter)
    counter <- 1
    nvars <- 1 # How many variables to enter data at each iteration
    
    # Enter prior evidence if present
    if (!is.null(evidence)) {
      bn <- setEvidence(bn, evidence = as.list(evidence))
      query <- setdiff(query, names(evidence))
    }
    
    
    for (i in 1:num_iter) {
      cat_result <- cat(
        bn,
        target,
        query,
        counter = counter,
        entropyThres = entropy_check_thres,
        miThres = mi_check_thres,
        maxIter = num_iter
      )
      
      if (cat_result[["stop"]]) {
        #print("Stopped")
        #print(cat_result$reason)
        post <- list(posterior = querygrain(bn, nodes = c(target), type = "marginal"),
                     evidence = data[entered_ev])
        post[["iter_posteriors"]] <- posterior_targets
        if (store_bns)
          post[["bns"]] <- bns
        return(post)
      }
      
      mi_vars <- cat_result[["mi"]]
      max_mi_vars <- names(mi_vars[1:nvars])
      if (mi_vars[1] == 0) {
        warning(paste("Maximum mutual information is 0 at iteration", i))
      }
      entered_ev <- unique(c(entered_ev, max_mi_vars))
      query_ev <- data[entered_ev]
      bn <-
        setEvidence(bn,
                    nodes = names(query_ev),
                    states = as.character(unlist(query_ev)))
      # Update query
      query <- setdiff(query, max_mi_vars)
      
      if (store_bns)
        bns[[i]] <- bn
      posterior <- querygrain(bn, nodes = c(target), type = "marginal")
      posterior_targets[[i]] <- posterior
      counter = counter + 1
      
    }
    post <- list(posterior = posterior, evidence = data[entered_ev])
    post[["iter_posteriors"]] <- posterior_targets
    if (store_bns)
      post[["bns"]] <- bns
    
    post
  }
