#JIVE.pred and JIVE.pred.predict() functions


#' Perform JIVE-predict
#'
#' Given multi-source data and an outcome, JIVE-predict will
#' identify shared (joint) and source-specific (individual)
#' underlying structure using Joint and Individual Variation Explained by
#' Lock et al. (2013). These structures will then be used to construct
#' a generalized linear model (GLM) for the outcome.
#'
#' @param X A list of two or more linked data matrices. Each matrix must
#' have the same number of columns which is assumed to be common.
#' @param Y An outcome expressed as a vector with length equal
#' to the number of columns in each view of \code{X}.
#' @param rankJ An integer specifying the joint rank of the data.
#' If \code{rankJ=NULL}, ranks will be determined by the \code{method} option.
#' @param rankA vector specifying the individual ranks of the data.
#' If \code{rankA=NULL}, ranks will be determined by the \code{method} option.
#' @param family A string specifying the type of prediction model to fit. Options
#' are "gaussian", "binomial", and "poisson". Model is fit using GLM.
#' @param center A boolean indicating whether or not X should be centered. Default is FALSE.
#' @param scale A boolean indicating whether or not X should be scaled. Default is FALSE.
#' @param orthIndiv A boolean indicating whether or not the algorithm should enforce
#' orthogonality between individual structures. Default is FALSE.
#' @param method A string with the method to use for rank selections.
#' Possible options are "given", "perm", and "bic". Default is "perm". If ranks are
#' specified, the method "given" will be used. See below for more details.
#' @param maxiter The maximum number of iterations allowed in the JIVE method.
#' @param showProgress A boolean indicating whether or not to give output showing
#' the progress of the algorithm.
#'
#' @details The method requires the data to be centered and scaled. This
#' can be done prior to running the method or by specifying center=T and scale=T.
#' The rank of the joint and individual components can be pre-specified
#' or adaptively selected within the function using JIVE's permutation approach
#' (\code{method="perm"}) or based on the Bayesian Information Criterion (\code{method="bic"}).
#' When \code{method="given"}, JIVE will use the user-specified ranks. If TRUE, \code{showProgress}
#' will print out updates about the number of iterations the JIVE algorithm is taking and the
#' progress of the rank selection method, if applicable.
#'
#' JIVE-predict extends JIVE to allow for supervision.
#' The function first runs JIVE to decompose multi-source data into low-rank,
#' orthogonal joint and individual components. Each component can be further broken down
#' into the loadings, or left eigenvectors, and the scores, the product of the
#' eigenvalues and the right eigenvectors. The number of eigenvectors is equal to
#' the rank of the component. JIVE-predict takes the scores from the fitted JIVE
#' model and uses them in a GLM function to predict \code{y}.
#'
#'
#' @return Returns an object of class JIVEpred. The function \code{summary}
#' (i.e. \code{\link{summary.JIVEpred}}) can be used to summarize the model results, including a
#' variance table and testing the significance of the joint and individual components.
#'
#' An object of class "JIVEpred" is a list containing the following components.
#'  \item{jive.fit}{A list of class "jive" containing the JIVE output. See the
#'  \code{jive} function from the r.jive R package for details about this object.}
#'  \item{mod.fit}{A list of class "glm" containing the GLM output. See the
#'  \code{glm} function from the stats R package for details about this object.}
#'  \item{data.matrix}{The data matrix that was used when fitting the GLM. The first column
#'  is the centered and scaled outcome, if applicable, and the remaining columns are the joint
#'  and individual score matrices.}
#'  \item{family}{A string stating which family was used when fitting the GLM}
#' @export
#'
#' @seealso \code{\link{predict.JIVEpred}}  \code{\link{summary.JIVEpred}}
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20), matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' train.fit <- JIVE.pred(X=train.x,Y=train.y,rankJ=1,rankI=c(1,1))
JIVE.pred <- function(X, Y, family="gaussian",
                      rankJ=NULL, rankA=NULL,
                      center=F, scale=F, orthIndiv=F,
                      method="perm", maxiter=1000,
                      showProgress=T){
  #family=c("gaussian", "binomial", "poisson")

  #Fit JIVE
  if(is.null(rankJ)|is.null(rankA)){
    time1 <- Sys.time()
    fit <- r.jive::jive(X, center=center, scale=scale,orthIndiv = orthIndiv,
                        method=method, maxiter=maxiter,showProgress=showProgress)
    time2 <- Sys.time()
  }else{
    time1 <- Sys.time()
    fit <- r.jive::jive(X, center=center, scale=scale,orthIndiv = orthIndiv, rankJ = rankJ,
                        rankA = rankA, method = "given",  maxiter=maxiter,
                        showProgress=showProgress)
    time2 <- Sys.time()
  }


  #Create Data matrix
  k <- length(X)
  jive.mat <- matrix(Y, ncol=1)
  nms <- "Y"
  #Add Joint Scores
  if(fit$rankJ >0){
    joint <- NULL
    for(i in 1:k){joint <- rbind(joint, fit$joint[[i]])}
    svd.temp <- svd(joint)
    t1 <- svd.temp$d[1:fit$rankJ] * t(svd.temp$v[,1:fit$rankJ])
    jive.mat <- cbind(jive.mat, t(t1))
    nms <- c(nms, paste0("J",1:fit$rankJ))
  }
  #Add Indiv Scores
  for(i in 1:k){
    if(fit$rankA[[i]]>0){
      svd.temp <- svd(fit$individual[[i]])
      t1 <- svd.temp$d[1:fit$rankA[i]] * t(svd.temp$v[,1:fit$rankA[i]])
      jive.mat <- cbind(jive.mat, t(t1))
      nms <- c(nms, paste0("I",i,"_", 1:fit$rankA[i]))
    }
  }
  jive.mat <- as.data.frame(jive.mat)
  names(jive.mat) <- nms

  # Fit Model
  if(family=="gaussian"){
    mod.fit <- stats::lm(Y ~ ., data=jive.mat)
  }else if(family=="binomial"){
    mod.fit <- stats::glm(Y ~ ., data=jive.mat, family = stats::binomial)
  }else if(family == "poisson"){
    mod.fit <- stats::glm(Y ~ ., data=jive.mat, family = stats::poisson)
  }else{
    stop(paste0("Model ", family, " doesn't exist"))
  }

  result <- list(jive.fit=fit,
                 mod.fit=mod.fit,
                 data.matrix=jive.mat,
                 family=family)
  class(result) <- "JIVEpred"


  return(result)

}


#' Prediction for JIVE-predict
#'
#' Predicted values based on the an JIVE-predict model.
#'
#' @param object An object of class "JIVEpred", usually a fitted JIVE-predict model.
#' @param newdata A list of matrices representing the new X datasets.
#' @param ... further arguments passed to or from other methods.
#'
#' @details \code{predict.JIVEpred} calculates predicted values for \code{newdata}
#' based on the fitted model. The function first calculates the joint and
#' individual score matrices for \code{newdata} treating the loading matrices from
#' the initial JIVE model as fixed.
#' Once the new score matrices are obtained, the GLM prediction will be
#' evaluated using the new scores as the data matrix.
#'
#' @return A list of the following components is returned:
#'  \item{Ypred}{The fitted Y values.}
#'  \item{Ynat}{The fitted values of Y's natural parameter space.}
#'  \item{joint.scores}{A matrix capturing the the joint scores of newdata.}
#'  \item{indiv.scores}{A list containing matrices that capture the individual scores of newdata.}
#'  \item{error}{The error value at which the model converged.}
#' @export
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20), matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' test.x <- list(matrix(rnorm(600), ncol=40),matrix(rnorm(400), ncol=40))
#' train.fit <- JIVE.pred(X=train.x,Y=train.y,rankJ=1,rankA=c(1,1))
#' test.fit <- predict(train.fit, test.x)
predict.JIVEpred <- function(object, newdata, ...){
  center=F
  scale=F
  n <- c()
  for (i in 1:length(newdata)) {
    n[i] <- nrow(newdata[[i]]) * ncol(newdata[[i]])
  }

  for (i in 1:length(newdata)) {
    if (center) {
      centerValues <- apply(newdata[[i]], 1, mean, na.rm = T)
      newdata[[i]] <- newdata[[i]] - matrix(rep(centerValues,
                                          ncol(newdata[[i]])), nrow = nrow(newdata[[i]]))
    }
    if (scale) {
      scaleValues <- norm(newdata[[i]], type = "f") * sqrt(sum(n))
      newdata[[i]] <- newdata[[i]]/scaleValues
    }
  }

  k <- length(newdata)
  jive.fit <- object$jive.fit
  mod.fit <- object$mod.fit
  fit_test2 <- r.jive::jive.predict(newdata, jive.fit)
  jive.mat <- cbind(matrix(rep(0, ncol(newdata[[1]])), ncol=1),
                    t(as.matrix(fit_test2$joint.scores)))
  for(i in 1:k){
    jive.mat <- cbind(jive.mat, t(as.matrix(fit_test2$indiv.scores[[i]])))
  }
  jive.mat <- as.data.frame(jive.mat)
  names(jive.mat) <-  names(object$data.matrix)
  y.pred <- stats::predict(mod.fit, newdata = jive.mat, type = "response")
  if(object$family != "gaussian"){
    y.pred2 <- round(y.pred)
  }else{y.pred2 <- y.pred}


  return(list(Ypred = y.pred2,
              Ynat = y.pred,
              joint.scores=fit_test2$joint.scores,
              indiv.scores = fit_test2$indiv.scores,
              error = fit_test2$errors))
}


#' Summarizing JIVE-predict Model Fits
#'
#' Summary methods for an JIVE-predict model of class "JIVEpred".
#'
#' @param object An object of class "JIVEpred", usually a fitted JIVE-predict model.
#'
#' @details Both the \code{print} and the \code{summary} functions
#' give summary results for a fitted JIVE-predict model.
#'
#' For the \code{summary} function, amount of variance explained
#' is expressed in terms of the standardized Frobenious
#' norm. Partial R-squared values are calculated for the
#' joint and individual components. If rank=1, a z-statistic
#' is calculated to determine the p-value. If rank>1, an F-statistic
#' is calculated.
#'
#' For the \code{print} function, the output displays the \code{print.glm}
#' result for the fitted GLM model.
#'
#' @return Summary measures
#' @export
summary.JIVEpred <- function(object, ...){
  k <- length(object$jive.fit$data)
  tbl_ranks <- data.frame(Source = c("Joint", paste0("Data", 1:k)),
                          Rank = c(object$jive.fit$rankJ, object$jive.fit$rankA))

  #cat("\n $Variance \n")
  var.table <- NULL
  for (i in 1:k) {
    j <- object$jive.fit$joint[[i]]
    a <- object$jive.fit$individual[[i]]
    ssj <- norm(j, type="f")^2
    ssi <- norm(a, type="f")^2
    sse <- norm(object$jive.fit$data[[i]] -j-a, type="f")^2
    sst <- norm(object$jive.fit$data[[i]], type="f")^2
    var.table <- cbind(var.table, round(c(ssj/sst, ssi/sst, sse/sst),4))
  }
  r_j <- object$jive.fit$rankJ
  j <- as.matrix(object$data.matrix[,c(2:(r_j+1))], ncol=r_j) %*% object$mod.fit$coefficients[c(2:(r_j+1))]
  ThetaS <- as.matrix(object$data.matrix[,-c(1:(r_j+1))]) %*% object$mod.fit$coefficients[-c(1:(r_j+1))]
  ssj <- norm(j, type="f")^2
  ssi <- norm(ThetaS, type="f")^2
  ypred <- j + ThetaS
  sse <- norm(as.matrix(object$data.matrix$Y-ypred), type="f")^2
  sst <- norm(as.matrix(object$data.matrix$Y), type="f")^2
  var.table <- cbind(c("Joint", "Indiv", "Error"),
                     var.table, round(c(ssj/sst, ssi/sst, sse/sst),4))
  var.table <- as.data.frame(var.table)
  names(var.table) <- c("Component", paste0("X", 1:k), "Y")
  if(object$family != "guassian"){
    var.table <- var.table[,-ncol(var.table)]
  }

  if(object$family == "guassian"){
  #cat("\n $pred.model \n")
  j <- as.matrix(object$data.matrix[,c(2:(r_j+1))], ncol=r_j) %*% object$mod.fit$coefficients[c(2:(r_j+1))]
  a <- list(); a2 <- 0
  sse <- sum((object$data.matrix$Y - j)^2)
  ssr <- sum((mean(object$data.matrix$Y) - j)^2)
  ssnum <-  sum((mean(j) - j)^2)
  coefs <- object$mod.fit$coefficients[2]
  for(i in 1:k){
    obs <- which(grepl(paste0("I",i),names(object$mod.fit$coefficients)))
    a[[i]] <- as.matrix(object$data.matrix[,obs]) %*% object$mod.fit$coefficients[obs]
    a2 <- a2 + a[[i]]
    sse <- c(sse, sum((object$data.matrix$Y - a[[i]])^2))
    ssr <- c(ssr, sum((mean(object$data.matrix$Y) - a[[i]])^2))
    ssnum <- c(ssnum, sum((mean(a[[i]]) - a[[i]])^2))
    coefs <- c(coefs, object$mod.fit$coefficients[obs[1]])
  }

  sse_full <- sum((object$data.matrix$Y - (j+a2))^2)
  sse_partial <- sum((object$data.matrix$Y - a2)^2)
  for(i in 1:k){
    temp <- j
    for(ii in 1:k){
      if(i != ii){temp <- temp + a[[i]]}
    }
    sse_partial <- c(sse_partial, sum((object$data.matrix$Y - temp)^2))
  }
  r_squared <- (sse_partial - sse_full)/sse_partial
  ranks <- c(object$jive.fit$rankJ, object$jive.fit$rankA)

  b <- which(ranks>1)
  n <- length(object$data.matrix$Y)
  if(length(b)>0){
    msr <- ssr[b] / (ranks[b]-1)
    mse <- sse[b] / (n-ranks[b])
    fstat <- msr/mse; pval <- NULL
    for(j in 1:length(b)){
      pval <- c(pval, 1-stats::pf(abs(fstat[j]), df1=ranks[b[j]]-1,
                                  df2=n-ranks[b[j]], lower.tail = T) )
    }}
  bb <- which(ranks==1)
  if(length(bb)>0){
    se <- sqrt((1/(n-1) * sse[bb])  / ssnum[bb] )
    z.stat <- coefs[bb]/se
    pval2 <- 2*(1-stats::pnorm(abs(z.stat), lower.tail = T))
  }
  pvalfinal <- ranks*1
  if(length(b)>0) pvalfinal[which(ranks>1)] <- pval
  if(length(bb)>0) pvalfinal[which(ranks==1)] <- pval2

  tbl <- data.frame(Component=c("Joint", paste0("Indiv", 1:k)),
                    Rank = ranks,
                    Partial_R2=r_squared,
                    Pvalue=pvalfinal)
  }else{
    tbl <- stats::anova(object$mod.fit)
  }

  return(list(ranks=tbl_ranks,
              variance=var.table,
              pred.model=tbl,
              Model=summary(object$mod.fit)))
}


#' @describeIn summary.JIVEpred
#'
#' @param x a fitted JIVEpred model.
#' @param ... further arguments passed to or from other methods.
print.JIVEpred <- function(x, ...) {
  k <- length(x$jive.fit$data)
  tbl_ranks <- data.frame(Source = c("Joint", paste0("Data", 1:k)),
                          Rank = c(x$jive.fit$rankJ, x$jive.fit$rankA))

  cat("Ranks: \n")
  print(tbl_ranks)
  cat("\n Model Fit: \n")
  print(x$mod.fit)
}

