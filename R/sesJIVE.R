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
#' @param center.scale A boolean for whether or not to center and scale.
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
#
#' @details The method requires the data to be centered and scaled. This
#' can be done prior to running the method or by specifying center.scale=T.
#' The rank of the joint and individual components,
#' the weight between the data and the outcome, and the lambda value for
#' sparsity can be pre-specified
#' or adaptively selected within the function. The method will print the
#' ranks, the weight, the lambda value, and the number of iterations needed
#' to reach convergence.
#'
#' The sesJIVE algorithm uses an iterative reweighted least squares (IRLS)
#' approach to solve for the parameters. The parameter estimates are
#' initialized by the \code{initial} option in the function. "uninformative"
#' will use random values (via the \code{rnorm} function) to initialize the
#' starting values. "svd" will take the singular value decomposition (SVD) of the
#' concatenated X matrix to initialize the joint components, and will take the SVD
#' of each individual X matrix to initialize the individual components. Lastly,
#' "jive" will run Lock et al.'s Joint and Variation Explained (JIVE) (2013) method
#' and use the model fit to initialize the parameters.
#'
#' sesJIVE extends JIVE and sJIVE to allow for different data distributions and
#' sparsity. It decomposes multi-source data into low-rank,
#' orthogonal joint and individual components in a generalized framework that allows
#' each X dataset to follow any exponential family distribution. Each component is broken down
#' into the loadings, or left eigenvectors, and the scores, the product of the
#' eigenvalues and the right eigenvectors. The number of eigenvectors is equal to
#' the rank of the component, and the scores are used to predict \code{y}. Sparsity is enforced
#' on the loadings using a LASSO penalty (Tibshirani, 1996), but the fitted score matrices
#' do not have any penalization.
#'
#' @return \code{sesJIVE} returns an object of class "sesJIVE". The function \code{summary}
#' (i.e. \code{\link{summary.sesJIVE}}) can be used to summarize the model results, including a
#' variance table and testing the significance of the joint and individual components.
#'
#' An object of class "sesJIVE" is a list containing the following components.
#'  \item{S_J}{A matrix capturing the joint scores of the data.}
#'  \item{S_I}{A list containing matrices that capture the individual scores of the data.}
#'  \item{U_I}{A list containing matrices that capture the joint loadings of the data.}
#'  \item{W_I}{A list containing matrices that capture the individual loadings of the data.}
#'  \item{theta1}{A vector that captures the effect of the joint scores on the outcome.}
#'  \item{theta2}{A list containing vectors that capture the effect of the individual scores on the outcome.}
#'  \item{fittedY}{The fitted Y values.}
#'  \item{error}{The error value at which the model converged.}
#'  \item{all.error}{The error value at each iteration.}
#'  \item{iterations}{The number of iterations needed to reach convergence.}
#'  \item{rankJ}{The rank of the joint structure.}
#'  \item{rankA}{The rank of the individual structure.}
#'  \item{eta}{The weight between the data and the outcome.}
#'  \item{data}{A list containing the centered and scaled data sets, if applicable.}
#'
#' @export
#'
#' @seealso \code{\link{predict.sesJIVE}}  \code{\link{summary.sesJIVE}}
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
  result <- "Hello World"
  class(result) <- "sesJIVE"
  return(result)
}




#' Prediction for sesJIVE
#'
#' Predicted values based on the an sesJIVE model.
#'
#' @param object An object of class "sesJIVE", usually a fitted sesJIVE model.
#' @param newdata A list of matrices representing the new X datasets.
#' @param threshold threshold for convergence
#' @param max.iter max iterations
#' @param show.error show error during iterations
#' @param show.message show confirmation when method converges
#'
#' @details \code{predict.sesJIVE} calculates predicted values for \code{newdata}
#' based on the fitted model. The function first calculates the joint and
#' individual score matrices for \code{newdata}. Note that the fitted model's
#' loadings and coefficients are treated as known and will not get re-calculated.
#' Once the new score matrices are obtained, the linear prediction model will be
#' evaluated using the new scores as the data matrix.
#'
#' @return A list of the following components is returned:
#'  \item{Ypred}{The fitted Y values.}
#'  \item{S_J}{A matrix capturing the joint scores of newdata.}
#'  \item{S_I}{A list containing matrices that capture the individual scores of newdata.}
#'  \item{iterations}{The number of iterations needed to reach convergence.}
#'  \item{error}{The error value at which the model converged.}
#' @export
predict.sesJIVE<- function(object, newdata, threshold = 0.00001,
                            max.iter=2000, show.error=F,
                            show.message=T){
  #deleted irls_iter_1
  print("Goodbye World")
}



#' Summarizing sesJIVE Model Fits
#'
#' Summary methods for an sesJIVE model of class "sesJIVE".
#'
#' @param object An object of class "sesJIVE", usually a fitted sesJIVE model.
#'
#' @details Both the \code{print} and the \code{summary} functions
#' give summary results for a fitted sesJIVE model.
#'
#' For the \code{summary} function, amount of variance explained
#' is expressed in terms of the standardized Frobenious
#' norm. Partial R-squared values are calculated for the
#' joint and individual components. If rank=1, a z-statistic
#' is calculated to determine the p-value. If rank>1, an F-statistic
#' is calculated.
#'
#' For the \code{print} function, the coeffecients are simply
#' printouts of theta1 and theta2 from the sesJIVE model.
#'
#' @return Summary measures
#' @export
summary.sesJIVE <- function(object, ...){
  print("Summarizing My World")
}


#' @describeIn summary.sesJIVE
#'
#' @param x a fitted sesJIVE model.
#' @param ... further arguments passed to or from other methods.
print.sesJIVE <- function(x, ...) {
  print("Printing My World")
}



####################### Helper Functions ##########################
## insert here
