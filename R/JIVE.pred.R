#JIVE.pred and JIVE.pred.predict() functions


#' Perform JIVE.Predict
#'
#' @param X A list of two or more linked data matrices. Each matrix must
#' have the same number of columns which is assumed to be common.
#' @param Y A numeric outcome expressed as a vector with length equal
#' to the number of columns in each view of \code{X}.
#' @param rankJ An integer specifying the joint rank of the data.
#' If \code{rankJ=NULL}, ranks will be determined by the \code{method} option.
#' @param rankI A vector specifying the individual ranks of the data.
#' If \code{rankA=NULL}, ranks will be determined by the \code{method} option.
#' @param family A string specifying the type of prediction model to fit. Options
#' are "gaussian", "binomial", and "poisson". Model is fit using GLM.
#' @param center A boolean indicating whether or not the data should be centered. Default is FALSE.
#' @param scale A boolean indicating whether or not the data should be scaled. Default is FALSE.
#' @param orthIndiv A boolean indicating whether or not the algorithm should enforce
#' orthogonality between individual structures. Default is FALSE.
#' @param method A string with the method to use for rank selections.
#' Possible options are "given", "perm", and "bic". Default is "perm". If ranks are
#' specified, the method "given" will be used.
#' @param maxiter The maximum number of iterations in the JIVE method
#' @param showProgress A boolean indicating whether or not to give output showing
#' the progress of the algorithm.
#'
#' @return Returns a fitted JIVE.predict model
#' @export
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20), matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' train.fit <- JIVE.pred(X=train.x,Y=train.y,rankJ=1,rankI=c(1,1))
JIVE.pred <- function(X, Y, family="gaussian",
                      rankJ=NULL, rankI=NULL,
                      center=F, scale=F, orthIndiv=F,
                      method="perm", maxiter=1000,
                      showProgress=T){
  #family=c("gaussian", "binomial", "poisson")

  #Fit JIVE
  if(is.null(rankJ)|is.null(rankI)){
    time1 <- Sys.time()
    fit <- r.jive::jive(X, center=center, scale=scale,orthIndiv = orthIndiv,
                        method=method, maxiter=maxiter,showProgress=showProgress)
    time2 <- Sys.time()
  }else{
    time1 <- Sys.time()
    fit <- r.jive::jive(X, center=center, scale=scale,orthIndiv = orthIndiv, rankJ = rankJ,
                        rankA = rankI, method = "given",  maxiter=maxiter,
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


  return(list(jive.fit=fit,
              mod.fit=mod.fit,
              data.matrix=jive.mat,
              family=family))

}



#' Prediction with JIVE.predict model
#'
#' @param JIVE.pred.fit Fitted JIVE.predict model
#' @param newdata A list of two or more linked data matrices. Each matrix must
#' have the same number of columns which is assumed to be common.
#'
#' @return Predictions for Y
#' @export
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20), matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' test.x <- list(matrix(rnorm(600), ncol=40),matrix(rnorm(400), ncol=40))
#' train.fit <- JIVE.pred(X=train.x,Y=train.y,rankJ=1,rankI=c(1,1))
#' test.fit <- JIVE.pred.predict(train.fit, test.x)
JIVE.pred.predict <- function(JIVE.pred.fit, newdata){
  k <- length(newdata)
  jive.fit <- JIVE.pred.fit$jive.fit
  mod.fit <- JIVE.pred.fit$mod.fit
  fit_test2 <- r.jive::jive.predict(newdata, jive.fit)
  jive.mat <- cbind(matrix(rep(0, ncol(newdata[[1]])), ncol=1),
                    t(as.matrix(fit_test2$joint.scores)))
  for(i in 1:k){
    jive.mat <- cbind(jive.mat, t(as.matrix(fit_test2$indiv.scores[[i]])))
  }
  jive.mat <- as.data.frame(jive.mat)
  names(jive.mat) <-  names(JIVE.pred.fit$data.matrix)
  y.pred <- stats::predict(mod.fit, newdata = jive.mat, type = "response")
  if(JIVE.pred.fit$family != "gaussian"){
    y.pred2 <- round(y.pred)
  }else{y.pred2 <- y.pred}


  return(list(Ypred = y.pred2,
              Ynat = y.pred,
              joint.scores=fit_test2$joint.scores,
              indiv.scores = fit_test2$indiv.scores,
              error = fit_test2$errors))
}
