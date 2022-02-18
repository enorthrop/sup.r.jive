#sesJIVE and sesJIVE.predict functions


#' Sparse Exponential Family Supervised JIVE (sesJIVE)
#'
#' Given multi-source data and an outcome, sesJIVE can
#' simultaneously identify shared (joint) and source-specific (individual)
#' underlying structure while building a prediction model for an outcome
#' using these structures. These two components are weighted to compromise between
#' explaining variation in the multi-source data and in the outcome, and the method
#' can enforce sparsity in the results if specified. Data and the outcome
#' can follow a normal, Bernoulli, or Poisson distribution.
#'
#' @param X A list of two or more linked data matrices. Each matrix must
#' have the same number of columns which is assumed to be common.
#' @param Y A numeric outcome expressed as a vector with length equal
#' to the number of columns in each view of \code{X}.
#' @param rankJ An integer specifying the joint rank of the data.
#' If \code{rankJ=NULL}, ranks will be determined by the \code{method} option.
#' @param rankA A vector specifying the individual ranks of the data.
#' If \code{rankA=NULL}, ranks will be determined by the \code{method} option.
#' @param wts A value or vector of values between 0 and 1. If \code{wts}
#' is a single value, \code{X} will be weighted by \code{wts} and \code{Y}
#' will be weighted by \code{1-wts}. if \code{wts} is a vector, 5-fold
#' CV will pick the \code{wts} that minimizes the test deviance.
#' @param max.iter The maximum number of iterations for each instance
#' of the sJIVE algorithm.
#' @param threshold The threshold used to determine convergence of the algorithm.
#' @param family.x A vector of length equal to the number of \code{X} data matrices
#' with each element specifying which exponential family the data follows. Options are
#' "gaussian", "binomial", or "poisson". Default is that all \code{X} matrices are gaussian
#' @param family.y A string specifying which exponential family the outcome follows.
#' Options are "gaussian", "binomial", or "poisson". Default is "gaussian".
#' @param numCores The number of cores to use when determining the optimal lambda.
#' Default is 1.
#' @param show.error A boolean indicating whether or not to display the weighted
#' log-likelihood after each iteration. Default is FALSE
#' @param sparse A boolean indication whether or not to enforce sparsity in the loadings.
#' See description below for more information.
#' @param lambda A value or vector indicating what values of lambda to consider. If
#' a vector of values, the optimal lambda will be chosen based on the method
#' selected in "method.lambda".
#' @param intercept A boolean indicating whether or not there should be an
#' intercept term in the results.
#' @param method.lambda Either "CV" or "BIC" specifying if the optimal lambda
#' should be selected by cross-validation or by the Bayesian Information Criterion.
#' Default is "CV".
#' @param method.wts A string with which rank selection method to use.
#' Possible options are "permute" which uses JIVE's permutation method,
#' or "CV" which uses 5-fold forward CV to determine ranks.
#' @param initial Either "uninformative", "svd", or "jive" indicating how to
#' generate the initial values for the algorithm. See description for more details.
#'
#' @return A fitted model of class sesJIVE
#' @export
sesJIVE <- function(X, Y, rankJ = 1, rankA=rep(1,length(X)),wts=NULL, max.iter=1000,
                    threshold = 0.001, family.x = rep("gaussian", length(X)),
                    family.y="gaussian", numCores=1, show.error=F, sparse=F,
                    lambda=seq(0.15,0.35,0.02), intercept=T,
                    method.lambda="CV", method.wts="permute",
                    initial="uninformative"){
  #Notes:
  # Added method.wts=c("permute", "CV")
  # Removed orthogonal=NULL, sparse.svd=F, and var.none=NULL, show.lambda=F
  print("Hello World")
}




#' Prediction for sesJIVE
#'
#' @param sesJIVE.fit Fitted sesJIVE model
#' @param newdata New dataset
#' @param threshold threshold for convergence
#' @param max.iter max iterations
#' @param show.error show error during iterations
#' @param show.message show confirmatinon when method converges
#'
#' @return
#' @export
sesJIVE.predict <- function(sesJIVE.fit, newdata, threshold = 0.00001,
                            max.iter=2000, show.error=F,
                            show.message=T){
  #deleted irls_iter_1
  print("Goodbye World")
}


####################### Helper Functions ##########################
## insert here
