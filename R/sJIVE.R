#sJIVE and sJIVE.predict functions

#' sJIVE Converge
#'
#' @param X A list of p_i by n matrices containing centered and scaled data
#' @param Y A length n vector of a centered and scaled outcome
#' @param eta The weight
#' @param max.iter The max interations
#' @param threshold The threshold
#' @param show.error Show error after each iteration
#' @param rankJ Joint rank
#' @param rankA Individual rank
#' @param show.message Show messages
#' @param reduce.dim Reduce dimensions to improve computation time
#' @param center.scale Center and scale X and Y
#'
#' @return A list of the fitted model
#' @export
#'
#' @examples
#' X <- 5
#' Y <- 3
#' sJIVE.converge(X,Y)
sJIVE.converge <- function(X, Y, eta=NULL, max.iter=1000, threshold = 0.001,
                           show.error =F, rankJ=NULL, rankA=NULL,
                           show.message=T, reduce.dim=T, center.scale=T){
 X + Y
}

