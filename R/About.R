

#' sup.r.jive: Supervised JIVE methods
#'
#' The sup.r.jive package performs supervised joint and individual
#' explained (JIVE) methods including JIVE.predict, supervised JIVE (sJIVE),
#' and sparse exponential family JIVE (sesJIVE). The package also includes
#' three plotting functions for visualizing the variance attributed to
#' each source and its effect with the outcome.
#'
#' @author Elise F. Palzer, Sandra E. Safo, and Eric F. Lock
#'
#' @references
#' Kaplan, A, and EF Lock. 2017. "Prediction With Dimension
#' Reduction of Multiple Molecular Data Sources for Patient  Survival."
#' Cancer Informatics 16:1-11.
#'
#' Lock, EF, KA Hoadley, JS Marron, and AB Nobel. 2013.
#' “Joint and Individual Variation Explained (JIVE) for Integrated
#' Analysis of Multiple Data Types.” The Annals of
#' Applied Statistics 7 (1): 523–42.
#'
#' Palzer, EF, C Wendt, R Bowler, CP Hersh, SE Safo, and EF Lock.
#' 2021. "sJIVE: Supervised Joint and Individual Variation Explained."
#' Pre-print on arXiv.
#'
#' @docType package
#' @name sup.r.jive-package
NULL


#' Multi-source simulated dataset with continuous outcome.
#'
#' A multi-source dataset with a normally distributed outcome,
#'  saved as a list. The first two elements of the list are
#'  the X data matrices, each with 30 predictors and 20 people.
#'  The third element of the list is a continuous outcome vector
#'  Y. Dataset was simulated using same constructed described in
#'  sJIVE by Palzer et al. (2022) with true rankJ=1, rankI=c(1,1),
#'  a large error in X (90\%) and a small error in Y (10\%).
#'
#' @format A list
#' \describe{
#'   \item{X}{A list of length 2. Each element in the list is a
#'    centered and scaled data matrix
#'   with 30 predictors (rows) and 20 people (columns)}
#'   \item{Y}{A centered and scaled outcome vector of length 20}
#'   ...
#' }
#' @references Palzer E.F., Wendt C., Bowler R., Hersh C., Safo S.E.,
#' Lock E.F. sJIVE: Supervised joint and individual variation explained.
#' \url{https://arxiv.org/pdf/2102.13278.pdf}
#'
"SimData.norm"


#' Multi-source simulated dataset with a binary outcome.
#'
#' A multi-source dataset with a binary outcome,
#'  saved as a list. The first two elements of the list are
#'  the X data matrices, each with 30 predictors and 20 people.
#'  The third element of the list is a binary outcome vector
#'  Y. Dataset was simulated using same constructed with
#'  true rankJ=1, rankI=c(1,1),
#'  a large error in X (90\%) and a large signal in Y.
#'  Only the first 15 predictors in each dataset predict Y.
#'
#' @format A list with 3 elements
#' \describe{
#'   \item{X}{A list of length 2. Each element in the list is a
#'    centered and scaled data matrix
#'   with 30 predictors (rows) and 20 people (columns)}
#'   \item{Y}{binary outcome vector of length 20}
#'   ...
#' }
"SimData.bin"



