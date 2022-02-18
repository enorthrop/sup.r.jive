#Combined predict function for JIVE.pred, sJIVE, and sesJIVE




#' Prediction using a Supervised JIVE model
#'
#' Computes updated scores and predicts an outcome for new data given
#' a fitted supervised joint and individual variation explained (JIVE)
#' model of class sJIVE, JIVE.pred, or sesJIVE.
#'
#' @param sup.JIVE.fit Fitted model. Can be a fitted JIVE.predict,
#' sJIVE, or sesJIVE model.
#' @param newdata a list of datasets containing the same predictors
#' as used for the fitted model with which to predict.
#' @param threshold threshold used to define convergence
#' @param max.iter maximum number of iterations
#'
#' @return A list with the following components
#' @export
sup.JIVE.predict <- function(sup.JIVE.fit, newdata,
                             threshold = 0.001, max.iter=2000){
  print("So Long World")
}
