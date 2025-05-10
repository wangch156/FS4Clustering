#' EM-test for mixture distributions
#' 
#' @description `em_test` performs the homogeneity test of several mixture models, including
#' mixtures of normal, Poisson and negative binomial (NB) distributions.
#' 
#' @param x A numerical vector. Data used for EM-test.
#' @param G An integer larger than or equals to 2. Number of groups in the mixture model.
#' @param dist A character string specifying the distribution used in the model, must be one of \code{"normal"}, \code{"Poisson"} or \code{"NB"} (default). See \sQuote{Details}.
#' @param alpha.init A numerical vector or matrix. If it is a vector, then its length should equal to \code{G}; If it is a matrix, then it should have exactly \code{G} columns and each row represents an initial value. Every element should be strictly positive and its rows will be automatically normalized. Note that this parameter is optional, the uniform mixing proportions will be automatically used as one initial value and will guarantee the validity of the test.
#' @param K An integer specifying number of iterations.
#' @param C A number specifying the penalty parameter. Usually, \code{1e-3}~\code{1e-5} will all work well. By default, \code{C=1e-4}.  See \sQuote{Details}.
#' @param labels A factor or vector representing pre-clustering labels used as prior for the EM-test. If this is \code{NULL}, no prior will be used.
#' @param prior.weight Numerical, specifying weight for the prior (is specified). If \code{labels} is \code{NULL}, this parameter will not have an effect on the test.
#' @param earlystop Numerical, specifying threshold of absolute changes in penalized log-likelihood (See \sQuote{Details}) for early-stopping. Empirically, if the data is far from homogeneously distributed, the penalized log-likelihood will increase rapidly in the early stage of iterations. Later iterations only increase the penalized likelihood marginally. Therefore, specifying an early-stopping threshold can save some time while still guarantee the validity of p-values.
#' @param limit_dist A character string specifying the limiting distribution used to compute p-values, must be one of \code{"chisq"} or \code{"precise"}. If \code{"chisq"} (default), a chi-square distribution will be used to generate valid p-values; if \code{"precise"}, precise limiting distribution will be used to generate precise p-values, but this will cost much more time.
#' @param cut_max An integer specifying error for gradient approximation. The approximation error is about \eqn{(cut_max)^{-2}}.
#' @returns A vector containing p-value and EM-statistic for the test.
#' @examples
#' # Example data
#' x <- rnbinom(100,mu = 2,size = 1)
#' y <- c(rnbinom(50,mu = 8,size = 1),rnbinom(50,mu = 2,size = 1))
#' # EM-test
#' em_test(x,2)
#' em_test(x,3)
#' em_test(y,2)
#' em_test(y,3)
#' em_test(y,2, dist = "Normal")
#' em_test(y,2, dist = "Poisson")
#' em_test(y,2, dist = "ZIP")
#' @details
#' This function performs EM-test to screen clustering-informative features as introduced in <\href{https://arxiv.org/abs/2306.12671}{arXiv:2306.12671}>.
#' The EM-test for normal mixtures was originally proposed by Li, Chen and Marriott (2009) <\href{https://www.jstor.org/stable/27798833}{doi:10.1093/biomet/asp011}>.
#' EM-test assumes that data are generated from a mixture distribution 
#' \deqn{\phi(x;\xi,\alpha) = \sum_{g=1}^{G}\alpha_g f(x;\theta_g),}
#' where \eqn{\xi = (\theta_1,...,\theta_G)}.
#' 
#' Largely speaking, we want to test whether all \eqn{\theta_g}'s are the same. (For more concrete definitions, see the two paper above.)
#' To achieve this, we use the EM algorithm to update the penalized log-likelihood \eqn{pl_n(\xi,\alpha)} under the heterogeneous model (not all \eqn{\theta_g}'s are the same), where
#' \deqn{pl_n(\xi,\alpha) = \sum_{i=1}^n \log \phi(x_i;\xi,\alpha) + C (\sum_{g=1}^G \log \alpha_g + G\log G).}
#' The hyper-parameter \eqn{C} was called the penalty parameter.
#' 
#' Currently, this function supports cases where \eqn{f(x;\theta)} is density of normal, Poisson or negative binomial distribution.
#' @export
#' 
em_test <- function(x,G,dist = "NB",alpha.init = NULL,K = 100,C = 1e-3,
                    labels = NULL,prior.weight = 0.05,earlystop = NA,limit_dist = "chisq", 
                    cut_max = 50){
  if (G<=1){
    stop("Number of components must be larger than 1.")
  }
  
  # prior
  if (is.null(labels)){
    use_prior <- F
    prior.weight <- 0
  }else {
    use_prior <- T
  }
  
  # limiting distribution
  if (limit_dist == "chisq"){
    use_ecdf <- F
  } else if (limit_dist == "precise"){
    use_ecdf <- T
  } else {
    stop("limit_dist must be 'chisq' or 'precise'.")
  }
  
  # tidy initial value for mixing proportions
  if (is.null(alpha.init)){
    alpha.init <- rep(1/G,G)
  } else{
    alpha.init <- rbind(alpha.init,rep(1/G,G))
  }
  # check length of alpha.init
  if (is.null(dim(alpha.init))){
    # alpha.ini is a vector
    if (length(alpha.init) != G){
      stop("Number of mixing proportions should equal to G.")
    }
    if (sum(alpha.init>0) < G){
      stop("All elements in alpha.ini should be strictly large than zero.")
    }
    alpha.init <- alpha.init/sum(alpha.init)
    alpha.init <- matrix(alpha.init[1:(G-1)],nrow = 1)
  } else{
    # alpha.ini is a matrix
    if (ncol(alpha.init) != G){
      stop("Number of mixing proportions should equal to G.")
    }
    if (sum(alpha.init>0) < G){
      stop("All elements in alpha.ini should be strictly large than zero.")
    }
    alpha.init <- t(apply(alpha.init,1, function(x) x/sum(x)))
    alpha.init <- alpha.init[,1:(G-1)]
  }
  
  # tidy labels
  labels <- factor(labels)
  levels(labels) <- 1:nlevels(labels)
  
  if (use_prior == F){
    labels <- factor(sample(1:G,length(x),replace = T))
  }
  # Check if prior labels have the same number of clusters as G
  if (use_prior == T & G != nlevels(labels)){
    stop("Classes of labels should be the same as G.")
  }
  
  if (dist == "NB"){
    em.stat(x,alpha.init,k0 = K,C = C,labels = labels,group.num = G,
            prior.weight = prior.weight,earlystop = earlystop, is_ecdf = use_ecdf, 
            cut_max = cut_max)
  } else if (dist == "Normal"){
    if (use_ecdf){
      stop("Not implemented yet! Please set limit_dist = 'chisq'")
    } else{
      em.stat.normal(x,alpha.init,k0 = K,C = C,labels = labels,group.num = G,
                     prior.weight = prior.weight,earlystop = earlystop, is_ecdf = use_ecdf)
    }
  } else if (dist == "Poisson"){
    if (use_ecdf){
      stop("Not implemented yet! Please set limit_dist = 'chisq'")
    } else{
      em.stat.poisson(x,alpha.init,k0 = K,C = C,labels = labels,group.num = G)
    }
  } else if (dist == "ZIP"){
    if (use_ecdf){
      stop("Not implemented yet! Please set limit_dist = 'chisq'")
    } else{
      em.stat.zero_inflated_poisson(x,alpha.init,k0 = K,C = C,labels = labels,group.num = G)
    }
  }
}