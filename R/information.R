#' Entropy of a probability distribution
#'
#' @param x A vector of numbers
#'
#' @return A number, entropy of the probability distribution
#' @export
#'
#' @examples
#' entropy(c(0.5,0.2,0.3))
#' entropy(c(0.5,0.5,0))

entropy <- function(x) {
  # to avoid log(0) = -Inf
  # 0 log(0) = 0
  x[x == 0] <- 1
  - sum(x * log(x))
}


#' Calculates entropy for BN nodes
#'
#' @param dist array. gRain posterior distribution
#' @param target array of target variable names
#'
#' @return array of numbers. Entropies.
#' @export
#'
#' @examples
#' library(gRain)
#' data(chest_cpt)
#' chest_bn <- grain(compileCPT(chest_cpt))
#' post <- querygrain(chest_bn, "dysp")
#' node_entropy(post)
dist_entropy <- function(dist, target = NULL) {
  marginal <- is.list(dist)
  if (marginal) {
    dist_variables <- names(dist)
  } else{
    dist_variables <- names(dimnames(dist))
  }
  
  if (is.null(target)) {
    target <- dist_variables
  }
  if (!all(target %in% dist_variables)) {
    stop(paste("Variable", target[!(target %in% dist_variables)], "is not present in the distribution"))
  }
  
  entropy <- 0
  
  if (marginal) {
    dist_entropy <- unlist(lapply(dist, entropy))
  } else {
    # Joint distribution
    dist_entropy <- entropy(dist)
    names(dist_entropy) <- paste(target, collapse = "_")
  }
  
  dist_entropy
}

#' Calculates conditional entropy
#'
#' @param dist array. gRain joint posterior distribution
#' @param target array of target variable names
#' @param conditional array of conditional variable names
#'
#' @return number. Conditional entropy
#' @importFrom gRbase tabExpand
#' @importFrom gRbase tabMarg
#' @export
#'
#' @examples
#' library(gRain)
#' data(chest_cpt)
#' chest_bn <- grain(compileCPT(chest_cpt))
#' target <- c("lung","asia")
#' conditional <- c("bronc","dysp")
#' cond_entropy(target, conditional, chest_bn)
cond_entropy <- function(dist,
                         target = NULL,
                         conditional = NULL) {
  dist_variables <- names(dimnames(dist))
  if (is.list(dist)) {
    stop("dist should be a posterior gRain distribution.")
  }
  if (is.null(target) & is.null(conditional)) {
    stop("Either target or conditional variable(s) need to be defined")
  }
  if (is.null(target)) {
    target <- setdiff(dist_variables, conditional)
  }
  if (is.null(conditional)) {
    conditional <- setdiff(dist_variables, target)
  }
  if (length(conditional) == 0) {
    stop("No conditional variables were defined.")
  }
  if (length(target) == 0) {
    stop("No target variables were defined.")
  }
  
  if (!setequal(unique(c(conditional, target)), unique(dist_variables))) {
    stop(
      "Variables in the posterior distirbution are not the same as those are in conditional and target."
    )
  }
  
  
  marg_cond <- tabMarg(dist, conditional)
  # Match dimensions
  marg_cond <- tabExpand(marg_cond, dimnames(dist))
  x <- dist * (log(dist) - log(marg_cond))
  x[dist == 0] <- 0
  - sum(x)
  
}

#' Mutual information
#'
#' Computes mutual information between two (sets of) variables using a gRain posterior joint distribution.
#'
#' @param dist array. gRain joint posterior distribution
#' @param var1 string or array. Variable names
#' @param var2 string or array. Variable names
#'
#' @return number. Mutual information
#' @importFrom gRbase tabMarg
#' @export
#'
#' @examples
mi <- function(dist, var1, var2 = NULL) {
  dist_variables <- names(dimnames(dist))
  cond_entropy <- cond_entropy(dist, var1, var2)
  target_marginal <- tabMarg(dist, var1)
  entropy(target_marginal) - cond_entropy
}




#' Mutual Informations Between a Target and a Set of Input Nodes from a BN
#'
#' This function returns the mutual informations between a set of input nodes and a target node.
#'
#' @param bn a gRain BN object
#' @param target a string, name of the target node
#' @param inputs optional, a vector, names of the input nodes. If NULL all nodes in the BN are used
#'
#' @return a vector of mutual information between the target and each of inputs
#' @importFrom gRain querygrain
#' @export
#'
#' @examples
#' library(gRain)
#' data(chest_cpt)
#' chest_bn <- grain(compileCPT(chest_cpt))
#' mi_bn(chest_bn, "lung")
#' mi_bn(chest_bn, "smoke", c("asia","dysp", "bronc"))
#'
mi_bn <- function(bn, target, inputs = NULL) {
  if (is.null(inputs)) {
    inputs <- bn$universe$nodes
  }
  
  mutInfs <- rep(NA, length(inputs))
  names(mutInfs) <- inputs
  
  # Mutual information of target nodes are NA
  mutInfs[target] <- NA
  # Remove target nodes from input nodes
  if (any(inputs %in% target))
    warning(
      "input nodes include target nodes. The mutual information of these nodes will be recorded as NA"
    )
  nodes <- inputs[which(!inputs %in% target)]
  
  # Mutual information of nodes with evidence are 0
  if (!is.null(bn$evidence$nodes)) {
    nodes <- nodes[!nodes %in% bn$evidence$nodes]
    # Previously we set it zero but sometimes mutual information of unobserved values are -1E10 which is smaller than zero
    # This caused the some non unobserved variables to never entered
    mutInfs[bn$evidence$nodes] <- -1E3
  }
  
  for (node in nodes) {
    dist <- querygrain(bn, c(target, node), type = "joint")
    mutInfs[node] <- mi(dist, target, node)
  }
  mutInfs
}

#' Sum of mutual information for a set of target nodes
#'
#' Plajner and Vomlel (2017)suggested using the sum of mutual information for marginal entropy (SMI) of each latent construct as an information measure for multiple latent constructs in CAT.
#' M. Plajner and J. Vomlel, “Student skill models in adaptive testing,” J. Mach. Learn. Res., vol. 52, no. 2016, pp. 403–414, 2016.
#'
#'
#' @param bn a gRain BN object
#' @param targets a vector, names of the target nodes
#' @param inputs optional, a vector, names of the input nodes. If NULL all nodes in the BN are used
#'
#' @return a vector of sum of mutual information between the targets and each of inputs
#' @importFrom gRain querygrain
#' @export
#'
#' @examples
#' library(gRain)
#' data(chest_cpt)
#' chest_bn <- grain(compileCPT(chest_cpt))
#' smi(chest_bn, c("lung","dysp"))
#'
smi <- function(bn, targets, inputs = NULL) {
  querygrain(bn, )
  smi <- 0
  for (target in targets) {
    smi <- smi + mi_bn(bn, target, inputs)
  }
  smi
}


#' CAT for BN
#'
#' Check various stopping rules. Returns marginals if stopping rule is satisfied, returns mutual informations if not
#'
#' @param targets a vector, names of the target nodes
#' @param bn a gRain BN object
#' @param inputs optional, a vector, names of the input nodes. If NULL all nodes in the BN are used
#' @param counter optional, a number, counter for keeping track of number of iterations
#' @param entropyThres optional, threshold value for the stopping rule of maximum entropy values in the target nodes
#' @param miThres optional, threshold value for the stopping rule of maximum mutual information values in the target nodes
#' @param maxIter optional, maximum number iterations for the iterations stopping rule
#'
#' @return list including posterior marginals of target variables if stopping rule is satisfied, ordered vector of input nodes according to their mutual information if stopping rule is not satisfied
#' @importFrom gRain querygrain
#' @export
#'
#' @examples
#' library(gRain)
#' data(chest_cpt)
#' chest_bn <- grain(compileCPT(chest_cpt))
#' cat(chest_bn,c("dysp","lung"),entropyThres = .1,counter=1)

cat <- function(bn,
                targets,
                inputs = NULL,
                counter = NULL,
                entropyThres = NULL,
                miThres = NULL,
                maxIter = NULL) {
  margPosterior <- querygrain(bn, nodes = targets, type = "marginal")
  maxTargetEntropy <- max(dist_entropy(margPosterior, targets))
  
  if ((!is.null(counter)) &&
      (!is.null(maxIter)) && counter > maxIter) {
    return(list(
      stop = TRUE,
      reason = paste("maximum number of iterations reached:", maxIter),
      posterior = margPosterior
    ))
  }
  
  if (((!is.null(counter)) &&
       is.null(maxIter)) ||
      (is.null(counter) &&
       (!is.null(maxIter))))
    warning(
      "You have defined one of counter or maxIter. In order to use maximum number of iterations stopping rule you need to define both"
    )
  
  # Stopping Rule: Entropy of Target Variables less than threshold
  if ((!is.null(entropyThres)) &&
      maxTargetEntropy <= entropyThres) {
    return(list(
      stop = TRUE,
      reason = paste("entropy of all target variables are less than", entropyThres),
      posterior = margPosterior
    ))
  }
  
  mi <- smi(bn, targets, inputs)
  # Sorting removes NA values (target variable) from inputs
  mi <- sort(mi, decreasing = TRUE)
  
  # Stopping Rule: No input left
  if (length(mi) == 0) {
    return(list(
      stop = TRUE,
      reason = "All inputs entered",
      posterior = margPosterior
    ))
  }
  
  # Stopping Rule: Mutual Information of Input Variables less than threshold
  if ((!is.null(miThres)) &&  max(mi) <= miThres) {
    return(list(
      stop = TRUE,
      reason = paste("mutual information of all inputs are less than", miThres),
      posterior = margPosterior
    ))
  }
  
  res <- list(stop = FALSE, mi = mi)
  
  # Update counter if it exists
  if (!is.null(counter)) {
    res[["counter"]] <- counter + 1
  }
  
  return(res)
}



#' KL Divergence of Two Vectors of Discrete Probability Distributions
#'
#' @param p array of numbers. A probability distribution.
#' @param q array of numbers. A probability distribution.
#'
#' @return number. Kullback-Leible divergence.
#' @export
#'
#' @examples kld(c(0.6,0.4),c(0.7,0.3))
kld <- function(p, q) {
  kl <- sum(p * log(p / q))
  if (is.nan(kl)) {
    kl <- 0
  }
  kl
}

#' Hellinger Distance of Two Vectors of Discrete Probability Distributions
#'
#' @param p array of numbers. A probability distribution.
#' @param q array of numbers. A probability distribution.
#'
#' @return number. Hellinger distance
#' @export
#'
#' @examples held(c(0.6,0.4),c(0.7,0.3))
held <- function(p, q) {
  # Hellinger Distance of Two Vectors of Discrete Probability Distributions
  hel <- sqrt(1 - sum(sqrt(p * q)))
  if (is.nan(hel)) {
    hel <- 0
  }
  hel
}
