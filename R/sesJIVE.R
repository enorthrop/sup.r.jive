#sesJIVE and sesJIVE.predict functions
#Using sesJIVE_Functions_V16.R Updated on 4/27/2022

`%dopar%` <- foreach::`%dopar%`

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
#' @param show.lambda A boolean indicating if an intermediate table should be printed
#' that shows the predictive performance of each candidate lambda value.
#' @param irls_iter The number of iterations the IRLS algorithm should run to update
#' each parameter in the algorithm. If NULL, will run until convergence.
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
#' @seealso \code{\link{predict.sesJIVE}}  \code{\link{summary.sesJIVE}}
#' @export
sesJIVE <- function(X, Y, rankJ = 1, rankA=rep(1,length(X)),wts=NULL, max.iter=1000,
                      threshold = 0.001, family.x = rep("gaussian", length(X)),
                      family.y="gaussian", numCores=1, show.error=F, sparse=F,
                      lambda= c(0.0001, 0.001, 0.005, 0.01, 0.025, 0.05, 1),
                      orthogonal=NULL, intercept=T,
                      show.lambda=F, method.lambda="CV", var.none=NULL,
                      initial="uninformative", irls_iter=1){
    ############################################################################
    #X is a list of 2 or more datasets, each with dimensions p_i by n
    #Y is continuous vector length n
    #wts is a tuning parameter between 0 and 1. When eta=NULL, a gridsearch
    #   is conducted to tune eta. You can specify a value of eta to use,
    #   or supply a vector of eta values for esJIVE to consider.
    #rankJ is a value for the low-rank of the joint component
    #rankA is a vector of the ranks for each X dataset.
    #method.lambda = CV or BIC
    ############################################################################

    k <- length(X)
    n <- ncol(X[[1]])
    if(length(Y) != n){stop("Number of columns differ between datasets")}

    if(is.null(orthogonal)){ #If Sparse model, don't enforce orthogonality
      orthogonal <- (sparse == F)
    }

    if(is.null(wts)){
      wt.vec=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
    }else{
      wt.vec=wts
    }
    lambda.dat <- NULL

    # Choose the best Weights
    if(length(wt.vec)>1){
      #get cv folds
      n <- length(Y)
      fold <- list()
      cutoff <- round(stats::quantile(1:n, c(.2,.4,.6,.8)))
      fold[[1]] <- 1:cutoff[1]
      fold[[2]] <- (cutoff[1]+1):cutoff[2]
      fold[[3]] <- (cutoff[2]+1):cutoff[3]
      fold[[4]] <- (cutoff[3]+1):cutoff[4]
      fold[[5]] <- (cutoff[4]+1):n

      cat("Choosing Tuning Parameter: Weights \n")
      doParallel::registerDoParallel(cores=numCores)
      test.best <- foreach::foreach(e=wt.vec, .combine=rbind) %dopar% {
        find.wts(e=e, YY=Y, XX=X,
                 max.iters=max.iter,
                 folds = fold,
                 family.xx = family.x, initials = initial,
                 family.yy = family.y, intercepts=intercept,
                 rankJJ=rankJ, rankAA=rankA)
      }
      #if(show.lambda){print(test.best)}
      doParallel::registerDoParallel(cores=1)
      best.wt <- which(test.best[,2] == min(test.best[,2], na.rm=T))
      best.wt <- best.wt[1]
      cat(paste0("Using wts= ", test.best[best.wt,1], "\n"))
    }else{
      test.best <- t(as.matrix(rep(wt.vec,2)))
      best.wt <- 1
    }

    # Choose the best lambda
    if(sparse & length(lambda)>1){
      if(method.lambda=="CV"){
        #get cv folds
        n <- length(Y)
        fold <- list()
        cutoff <- round(stats::quantile(1:n, c(.2,.4,.6,.8)))
        fold[[1]] <- 1:cutoff[1]
        fold[[2]] <- (cutoff[1]+1):cutoff[2]
        fold[[3]] <- (cutoff[2]+1):cutoff[3]
        fold[[4]] <- (cutoff[3]+1):cutoff[4]
        fold[[5]] <- (cutoff[4]+1):n

        cat("Choosing Tuning Parameter: Lambda \n")

        #no sparsity -- for variance estimate
        if(is.null(var.none)){ #find.lambda function doesn't work with var0 yet...
          test.zero <- find.lambda(lambda=0, YY=Y, XX=X,
                                   max.iters=max.iter,
                                   folds = fold,initials = initial,
                                   weights=c(rep(test.best[best.wt,1],length(X)), 1-test.best[best.wt,1]),
                                   family.xx = family.x, intercepts=intercept,
                                   family.yy = family.y, ortho=orthogonal,
                                   rankJJ=rankJ, rankAA=rankA)

          var.none <- test.zero[(length(test.zero)-k):length(test.zero) ]
        }
        #with sparsity -- using variance estimate
        test.best2 <- foreach::foreach(lambda=lambda, .combine=rbind) %dopar% {
          #for(i in 1:length(lambda)){
          find.lambda(lambda=lambda, YY=Y, XX=X,
                      max.iters=max.iter,
                      folds = fold,initials = initial,
                      weights=c(rep(test.best[best.wt,1],length(X)), 1-test.best[best.wt,1]),
                      family.xx = family.x, intercepts=intercept,
                      family.yy = family.y, ortho=orthogonal,
                      rankJJ=rankJ, rankAA=rankA, var0= var.none)
          #print(paste0("Made it past ", i))
          #}
        }
        doParallel::registerDoParallel(cores=1)
        if(show.lambda){print(test.best2)}
        lambda.dat <- test.best2
        temp <- which(test.best2[,2] == min(test.best2[,2], na.rm=T))
        #temp2 <- which(test.best[,2] < test.best[temp,2]+test.best[temp,3])
        best.lambda <- max(test.best2[temp,1])
        cat(paste0("Using lambda= ", best.lambda, "\n"))

      }else{
        cat("Choosing Tuning Parameter: Lambda \n")
        #no sparsity -- for variance estimate
        if(is.null(var.none)){
          test.zero <- find.lambda(lambda=0, YY=Y, XX=X,
                                   max.iters=max.iter,initials = initial,
                                   folds = fold, method=method.lambda,
                                   weights=c(rep(test.best[best.wt,1],length(X)), 1-test.best[best.wt,1]),
                                   family.xx = family.x, intercepts=intercept,
                                   family.yy = family.y, ortho=orthogonal,
                                   rankJJ=rankJ, rankAA=rankA, var0= var.none)

          var.none <- test.zero[(length(test.zero)-k):length(test.zero) ]
        }
        #with sparsity -- using variance estimate
        doParallel::registerDoParallel(cores=numCores)
        test.best2 <- foreach::foreach(lambda=lambda, .combine=rbind) %dopar% {
          #for(i in 1:length(lambda)){
          find.lambda(lambda=lambda, YY=Y, XX=X,
                      max.iters=max.iter,
                      method=method.lambda,initials = initial,
                      weights=c(rep(test.best[best.wt,1],length(X)), 1-test.best[best.wt,1]),
                      family.xx = family.x, intercepts=intercept,
                      family.yy = family.y, ortho=orthogonal,
                      rankJJ=rankJ, rankAA=rankA, var0= var.none)
          #print(paste0("Made it past ", i))
          #}
        }
        doParallel::registerDoParallel(cores=1)
        if(show.lambda){print(test.best2)}
        lambda.dat <- test.best2
        temp <- which(test.best2[,2] == min(test.best2[,2], na.rm=T))
        best.lambda <- max(test.best2[temp,1])
        cat(paste0("Using lambda= ", best.lambda, "\n"))
      }
    }else{
      best.lambda <- lambda
      if(sparse){var.none <- rep(1,k+1)}
    }

    cat("Estimating Loadings and Scores \n")
    test.best <- sesJIVE.converge(X, Y,
                                  max.iter=max.iter, threshold = threshold,
                                  family.x = family.x,
                                  family.y = family.y,
                                  rankJ=rankJ, rankA=rankA,
                                  weights=c(rep(test.best[best.wt,1],length(X)), 1-test.best[best.wt,1]),
                                  show.message=T, show.error=show.error, var00=var.none,
                                  irls_iter=irls_iter, intercept=intercept, sparse=sparse,
                                  lambda=best.lambda, orthogonal=orthogonal,
                                  initial=initial)
    if(sparse){
      cat("Re-estimating Scores \n")
      test.best.pred <- stats::predict(test.best, X, show.error=show.error,
                                        train=T)

      x.mat <- cbind(matrix(Y, ncol=1), t(test.best.pred$Sj))
      for(i in 1:k){ x.mat <- cbind(x.mat, t(test.best.pred$Si[[i]]))}
      x.mat <- as.data.frame(x.mat)
      names(x.mat) <- paste0("V", 1:ncol(x.mat))
      fit <- stats::glm(V1 ~ ., data=x.mat, family = family.y)

      #Combine results into final results
      test.final <- test.best
      test.final$S_J <- test.best.pred$Sj
      test.final$S_I <- test.best.pred$Si
      test.final$theta1 <- stats::coef(fit)[2:(1+rankJ)]
      bad.obs <- which(is.na(test.final$theta1))
      if(length(bad.obs)>0){test.final$theta1[bad.obs] <- 0}
      obs <- 2+rankJ
      for(i in 1:k){
        test.final$theta2[[i]] <- stats::coef(fit)[obs:(obs-1+rankA[i])]
        bad.obs <- which(is.na(test.final$theta2[[i]]))
        if(length(bad.obs)>0){test.final$theta2[[i]][bad.obs] <- 0}
        obs <- obs+rankA[i]
      }
      #Re-calculate natX and natY
      int <- t(as.matrix(rep(1,ncol(test.final$natX[[1]]))))
      thetaS <- 0; X.tilde <- NULL
      for(i in 1:k){
        X.tilde <- rbind(X.tilde, X)
        t1 <- test.final$intercept[[i]] %*%  int
        t2 <- test.final$U_I[[i]] %*% test.final$S_J
        t3 <- test.final$W_I[[i]] %*% test.final$S_I[[i]]
        test.final$natX[[i]] <- t1 + t2 + t3
        thetaS <- thetaS + test.final$theta2[[i]] %*% test.final$S_I[[i]]
      }
      test.final$natY <- test.final$intercept[[k+1]] %*% int +
        test.final$theta1 %*% test.final$S_J + thetaS
      #Update Error
      error.new <- get_loglik(test.final, X, Y, family.x, family.y)
      test.final$error <- error.new$log_lik
    }else{
      test.final <- test.best
    }

    if(sparse & length(lambda)>1){
      if(method.lambda=="CV"){
        lambda.dat <- as.data.frame(lambda.dat[,1:3])
        row.names(lambda.dat) <- c()
        names(lambda.dat) <- c("lambda", "Test MSE", "SE")
      }else if(method.lambda=="AIC"){
        lambda.dat <- as.data.frame(lambda.dat[,1:2])
        row.names(lambda.dat) <- c()
        names(lambda.dat) <- c("lambda", "AIC")
      }else if(method.lambda=="BIC"){
        lambda.dat <- as.data.frame(lambda.dat[,1:2])
        row.names(lambda.dat) <- c()
        names(lambda.dat) <- c("lambda", "BIC")

      }
      test.final$lambda.dat <- lambda.dat
    }

    if(length(var.none)>1){var.none = var.none[length(var.none)]}
    dev.resid <- get_deviance(Y, test.final$natY, family.y, var0=var.none)

    test.final$deviance <- dev.resid

    return(test.final)

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
#' @param irls_iter the number of iterations of the IRLS algorithm that should
#' be performed for each parameter update within the algorithm. If NULL, will
#' run until convergence.
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
                              max.iter=2000, irls_iter=1, show.error=F,
                              show.message=T, train=F, ...){
    ##############################################
    # object is the output from sesJIVE
    # newdata is list with the same predictors and
    #     number of datasets as used in sJIVE.fit
    ##############################################
    U.norm <- object$U.norm
    W.norm <- object$W.norm

    if(U.norm == 0){
      sparse <- FALSE
    }else{
      sparse <- TRUE
    }

    if(object$rankJ==0 & sum(object$rankA)==0){
      return(list(Ypred = 0,
                  Sj = 0,
                  Si = 0,
                  iteration = 0,
                  error = NA))
    }

    #Initalize values
    weights <- object$weights
    k <- length(newdata)
    n <- ncol(newdata[[1]])
    W <- object$W_I
    U <- object$U_I
    mu <- object$intercept
    int <- t(as.matrix(rep(1,n)))
    rankJ <- ncol(as.matrix(U[[1]]))
    if(train){
      Sj <- object$S_J
    }else{
      Sj <- matrix(rep(0,rankJ*n), ncol = n)
    }

    if(sparse){
      for(i in 1:k){
        U[[i]] <- U.norm * U[[i]]
        W[[i]] <- W.norm[[i]] * W[[i]]
      }
    }

    obs <- rankA <- Si <- list(); temp <- 0; X.tilde <- NULL
    for(i in 1:k){
      max.obs <- max(temp)
      temp <- (max.obs+1):(max.obs+nrow(newdata[[i]]))
      obs[[i]] <- temp

      X.tilde <- rbind(X.tilde, newdata[[i]])

      rankA[[i]] <- ncol(as.matrix(W[[i]]))
      if(train){
        Si[[i]] <- object$S_I[[i]]
      }else{
        Si[[i]] <- matrix(rep(0, rankA[[i]]*n), ncol=n)
      }
    }

    fam.list <- list()
    for(i in 1:k){
      if(object$family.x[i]=="gaussian"){
        fam.list[[i]] <- stats::gaussian()
      }else if(object$family.x[i]=="binomial"){
        fam.list[[i]] <- stats::binomial()
      }else if(object$family.x[i]=="poisson"){
        fam.list[[i]] <- stats::poisson()
      }else{
        print(paste0(object$family.x[i], " Distribution Does Not Exist"))
        stop()
      }
    }
    if(object$family.y=="gaussian"){
      famY <- stats::gaussian()
    }else if(object$family.y=="binomial"){
      famY <- stats::binomial()
    }else if(object$family.y=="poisson"){
      famY <- stats::poisson()
    }else{
      print(paste0(object$family.y, " Distribution Does Not Exist"))
      stop()
    }

    #Get Error
    error.old <- sesJIVE.error(X.tilde, U, Sj, W, Si, k, mu, fam.list, obs, kk=k,
                               wt.vec=weights)

    #Set up IRLS
    irls.list <- list()
    for(i in 1:(k+1)){
      if(i==1){p <- "Sj"; dfs <- 1:k
      beta <- Sj
      #matrix(rep(0,rankJ*n), ncol=n)   #rnorm(rankJ*n,0,0.1), ncol=n)
      }else{p <- paste0("S",i-1); dfs <- i-1
      beta <- Si[[i-1]]
      #matrix(rep(0,rankA[[i-1]]*n), ncol=n) #rnorm(rankA[[i-1]]*n,0,0.1), ncol=n)
      }
      irls.list[[i]] <- list(param=p,
                             dfs = dfs,
                             eta.old=0,
                             eta.new=0,
                             beta.old=beta,
                             beta.new=beta,
                             dev.old = 0,
                             dev.new = 0,
                             iter=0,
                             dev.all = 0)
    }
    eta.temp <- NULL
    for(i in 1:k){
      t1 <- mu[[i]] %*% t(as.matrix(rep(1,ncol(as.matrix(Sj))))) +
        U[[i]] %*% matrix(Sj, ncol=n) + W[[i]] %*% Si[[i]]
      eta.temp <- rbind(eta.temp, t1)
    }
    irls.list[[k+2]] <- eta.temp
    fit <- k+2

    ############################ Loop ################
    for(iter in 1:max.iter){

      #Optimize Sj
      U.mat <- A <- mu.mat <- NULL
      for(i in 1:k){
        mu.mat <- rbind(mu.mat, as.matrix(mu[[i]]))
        U.mat <- rbind(U.mat, as.matrix(U[[i]]))
        A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
      }
      off <- mu.mat %*% int + A
      irls.list <- irls_func(irls.list,U.mat, #U is known
                             num_iter = irls_iter, list_num = 1, offsets=off,
                             outcome = X.tilde, thresholds=threshold, famlist=fam.list,
                             transpose=F,  ob=obs, score=sparse, predicting = T,
                             eta1 = irls.list[[k+2]], Xtilde=X.tilde,
                             wt.vec=object$weights)
      Sj <- irls.list[[1]]$beta.new


      for(i in 1:k){
        #Calculate Outcome and Predictor
        W_temp <- matrix(rep(0,rankA[[i]]*nrow(X.tilde)), ncol=rankA[[i]])
        W_temp[obs[[i]],] <- W[[i]]
        off <- mu.mat %*% int + U.mat %*% Sj

        #Optimize Si
        irls.list <- irls_func(irls.list,predictor=W_temp, #W is known
                               num_iter = irls_iter, list_num = i+1, outcome = X.tilde,
                               offsets=off, Xtilde = X.tilde,
                               thresholds=threshold, famlist=fam.list,
                               transpose=F, ob=obs, predicting = T,
                               eta1 = irls.list[[fit]], score=sparse,
                               wt.vec=object$weights)
        Si[[i]] <- irls.list[[i+1]]$beta.new
      }

      #Get Error
      #Figure out the error
      error.new <- sesJIVE.error(X.tilde, U, Sj, W, Si, k, mu,
                                 fam.list,obs, kk=k, wt.vec=weights)
      if(show.error){print(error.new$log_lik)}
      #Check for Convergence
      if(abs(error.old$log_lik - error.new$log_lik) < threshold){
        if(show.message){cat(paste0("Converged in ", iter, " iterations \n"))}
        break
      }else{
        error.old <- error.new
      }
    }

    Xnat <- list()
    for(i in 1:k){
      Xnat[[i]] <- object$intercept[[i]] %*% int + U[[i]] %*% Sj +
        W[[i]] %*% Si[[i]]
    }

    Ynat <- object$intercept[[k+1]] %*% int + object$theta1 %*% Sj
    for(i in 1:k){
      Ynat <- Ynat + object$theta2[[i]] %*% Si[[i]]
    }
    Ypred <- Yprob <- famY$linkinv(Ynat)
    if(object$family.y=="binomial"){
      Ypred <- round(Ypred,0)
    }else if(object$family.y=="poisson"){
      Ypred <- round(Ypred,0)
    }

    #If sparse, re-scale Scores/loadings
    if(sparse){
      Sj <- Sj * object$U.norm
      for(i in 1:k){
        Si[[i]] <- Si[[i]] * object$W.norm[[i]]
      }
    }


    return(list(Ypred = Ypred,
                Ynat = Ynat,
                Xnat = Xnat,
                Yprob = Yprob,
                Sj = Sj,
                Si = Si,
                iteration = iter,
                error = error.new$log_lik))
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

backtrack <- function(x, dx, f, df, alpha=0.01, beta=0.8, b.full,
                      b.sparse, t) {
  b.full <- as.matrix(b.full)
  b.sparse <- as.matrix(as.vector(b.sparse))
  t <- 0.1
  g <- df(x)
  u <- alpha*sum(g*dx)
  #kk <- 1
  repeat {
    if(x-t*dx>=0 & x-t*dx<=1){
      if (f(as.numeric(x - t*dx)) <= f(x) + t*u ) break
    }
    t <- beta*t
    #kk <- kk + 1
  }
  return(c(x - t*dx, t))
}


optim.error2 <- function(irlslist2, famlist2, kk2, ob2,
                        Xtilde2, wt.vec2, sparse, var0,
                        lambda2){
  Sj <- irlslist2[[1]]$beta.new
  U <- irlslist2[[2]]$beta.new
  mu.temp <- as.matrix(irlslist2[[length(irlslist2)-1]]$beta.new) %*% rep(1,ncol(Sj))
  WS <- NULL; thetaS <- 0
  penalty <- norm(Sj, type="F") + lambda2 * norm(matrix(U), type="O")
  for(i in 1:kk2){
    temp <- irlslist2[[i*2+2]]$beta.new %*% irlslist2[[i*2+1]]$beta.new
    WS <- rbind(WS, temp[-nrow(temp),])
    thetaS <- thetaS + temp[nrow(temp),]
    Si_pen <- norm(irlslist2[[i*2+1]]$beta.new, type="F")
    Wi_pen <- norm(irlslist2[[i*2+2]]$beta.new, type="O")
    penalty <- penalty + Si_pen +  lambda2 *Wi_pen
  }
  Y.pred <-mu.temp + U %*% Sj + rbind(WS,thetaS)

  #Calculate Likelihood
  data_ll2 <- NULL
  for(i in 1:(kk2+1)){
    X <- Xtilde2[ob2[[i]],]
    natX <- Y.pred[ob2[[i]],]
    Xfit <- famlist2[[i]]$linkinv(natX)
    n <- ncol(natX)
    if(is.null(n)){n <- 1}
    if(famlist2[[i]]$family=="gaussian"){
      if(is.null(var0)){
        sigma2 <- 1 #var(as.vector(X - Xfit))
      }else{ sigma2 <- var0[i]}
      ll <- -wt.vec2[i] * sum((X - Xfit)^2) / (2 * sigma2)
    }else if(famlist2[[i]]$family=="binomial"){
      ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
    }else if(famlist2[[i]]$family=="poisson"){
      fx <- log(factorial(X))
      fx[which(fx==Inf)] <- -55.22 + 4.32*X[which(fx==Inf)]
      ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
    }
    data_ll2 <- c(data_ll2, ll)
  }
  #############
  #Get Deviance
  data_dev <- NULL
  for(i in 1:(kk2+1)){
    X <- Xtilde2[ob2[[i]],]
    natX <- Y.pred[ob2[[i]],]
    Xfit <- famlist2[[i]]$linkinv(natX)
    dev_temp <- 0

    for(j in 1:length(ob2[[i]])){
      if(is.vector(X)){
        ytrue <- X
        ynat <- natX
      }else{
        ytrue <- as.vector(X[j,])
        ynat <- as.vector(natX[j,])
      }
      n <- length(ytrue)
      if(famlist2[[i]]$family=="gaussian"){
        if(is.null(var0)){
          sigma2 <- 1
        }else{ sigma2 <- var0}
        dev.resid <-  2*sum((ytrue - ynat)^2) / (2 * sigma2)
      }else if(famlist2[[i]]$family=="binomial"){
        p.hat <- exp(ynat) / (1 + exp(ynat))
        dev.resid <- -2*( ytrue *log(p.hat) + (1-ytrue)*log(1-p.hat))
        t1 <- which(dev.resid == Inf)
        t2 <- which(is.na(dev.resid))
        if(length(t1)>0){dev.resid[t1] <- 100000}
        if(length(t2)>0){dev.resid[t1] <- 0}
        dev.resid <- sum(dev.resid, na.rm = T)
      }else if(famlist2[[i]]$family=="poisson"){
        mu <- exp(ynat)
        zero.obs <- which(ytrue==0)
        if(length(zero.obs)>0){
          dev.resid <- 2*sum( ytrue[-zero.obs]*log(ytrue[-zero.obs]/mu[-zero.obs]) -
                                ytrue[-zero.obs] + mu[-zero.obs]) + 2*sum(mu[zero.obs])
        }else{dev.resid <- 2*sum( ytrue*log(ytrue/mu) - ytrue + mu)}
      }
      dev_temp <- dev_temp + dev.resid

    }
    data_dev <- c(data_dev, dev_temp)
  }

  #############
  if(sparse==F){penalty<-0}

  return(list(log_lik = sum(data_ll2),#-penalty,
              dev=sum(data_dev),
              data_lik = data_ll2))
}

get_loglik <- function(sesJIVE.fit, XX, YY, family.xx, family.yy){
  #Calculate Likelihood
  kk2 <- k <- length(XX)

  famlist <- list()
  family.x <- c(family.xx, family.yy)
  for(i in 1:(k+1)){
    if(family.x[i]=="gaussian"){
      famlist[[i]] <- stats::gaussian()
    }else if(family.x[i]=="binomial"){
      famlist[[i]] <- stats::binomial()
    }else if(family.x[i]=="poisson"){
      famlist[[i]] <- stats::poisson()
    }else{
      print(paste0(family.x[i], " Distribution Does Not Exist"))
      stop()
    }
  }
  wt.vec2 <- sesJIVE.fit$weights
  data_ll2 <- NULL
  for(i in 1:(kk2+1)){
    if(i <= kk2){
      X <- XX[[i]]
      natX <- sesJIVE.fit$natX[[i]]
    }else{
      X <- YY
      natX <- sesJIVE.fit$natY
    }
    Xfit <- famlist[[i]]$linkinv(natX)
    n <- ncol(natX)
    if(is.null(n)){n <- 1}
    if(famlist[[i]]$family=="gaussian"){
      ll <- -wt.vec2[i] * sum((X - Xfit)^2) / 2
    }else if(famlist[[i]]$family=="binomial"){
      ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
    }else if(famlist[[i]]$family=="poisson"){
      fx <- log(factorial(X))
      fx[which(fx==Inf)] <- -55.22 + 4.32*X[which(fx==Inf)]
      ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
    }
    data_ll2 <- c(data_ll2, ll)
  }
  #############
  return(list(log_lik = sum(data_ll2),
              data_lik = data_ll2))
}

get_deviance <- function(ytrue, ynat, family.yy, var0=NULL){
  ytrue <- as.vector(ytrue)
  ynat <- as.vector(ynat)
  n <- length(ytrue)
  if(family.yy=="gaussian"){
    if(is.null(var0)){
      sigma2 <- 1
    }else{ sigma2 <- var0}
    dev.resid <-  2*sum((ytrue - ynat)^2) / (2 * sigma2)
  }else if(family.yy=="binomial"){
    p.hat <- exp(ynat) / (1 + exp(ynat))
    dev.resid <- -2*( ytrue *log(p.hat) + (1-ytrue)*log(1-p.hat))
    t1 <- which(dev.resid == Inf)
    t2 <- which(is.na(dev.resid))
    if(length(t1)>0){dev.resid[t1] <- 100000}
    if(length(t2)>0){dev.resid[t1] <- 0}
    dev.resid <- sum(dev.resid, na.rm = T)
  }else if(family.yy=="poisson"){
    mu <- exp(ynat)
    zero.obs <- which(ytrue==0)
    if(length(zero.obs)>0){
      dev.resid <- 2*sum( ytrue[-zero.obs]*log(ytrue[-zero.obs]/mu[-zero.obs]) -
                            ytrue[-zero.obs] + mu[-zero.obs]) + 2*sum(mu[zero.obs])
    }else{dev.resid <- 2*sum( ytrue*log(ytrue/mu) - ytrue + mu)}
  }
  return(dev.resid)
}

irls_func <- function(irlslist, predictor, offsets, list_num, num_iter=1,
                      thresholds, famlist,
                      outcome, transpose=F, ob, predicting=F,
                      eta1=NULL, wt.vec, sparse=F, lambda=1, kk,
                      Xtilde, old.err=NULL,scale=F, score=F, var0,
                      full.obs=obs){
  dat <- irlslist[[list_num]]
  time1 <- Sys.time()
  time8 <- time7 <- time6 <- time5 <- time4 <- time3 <- time2 <- time1
  bls1 <- bls2 <- 0
  bls.t <- NULL

  converge <- F
  if(is.null(num_iter)){
    num_iter <- 1000
    converge <- T
  }
  #update.beta <- T
  for(i in 1:num_iter){
    if(i>1 & transpose){eta1 <- t(eta1)}
    beta1 <- NULL; dev1 <- 0
    W <- NULL; z <- NULL; rownum <- NULL
    #big.eta <- NULL; eta.old2 <- NULL
    for(j in dat$dfs){
      if(transpose){
        X <- t(as.matrix(predictor))
        nobs = nrow(X)
        Y <- t(outcome[ob[[j]],])
        off <-  t(offsets[ob[[j]],])
        if(length(ob[[j]])==1){
          Y <- as.matrix(as.numeric(Y), ncol=1)
          off <- as.matrix(as.numeric(off), ncol=1)
        }
        nvars = ncol(Y)
      }else{
        X <- as.matrix(predictor[ob[[j]],]) ##Check if dimensions check out when
        nobs = length(ob[[j]])               ### nobs=1 and rank>1
        Y <- outcome[ob[[j]],]
        off <- offsets[ob[[j]],]
        if(length(ob[[j]])==1){
          Y <- t(as.matrix(as.numeric(Y), ncol=1))
          off <- t(as.matrix(as.numeric(off), ncol=1))
          X <- t(as.matrix(as.numeric(X), ncol=1))
        }
        nvars = ncol(Y)
      }
      if(transpose==F){
        wt <- rep(wt.vec[j],nobs)
      }else{
        wt <- rep(1,nobs)
      }
      if(dat$iter==0 &
         ((predicting==F & list_num==(length(irlslist)-1)) | #sesJIVE=finding mu
          (predicting==T & list_num==1))){   #sesJIVE.predict=finding Sj
        if(famlist[[j]]$family =="binomial"){
          eta.old <-  famlist[[j]]$linkfun((Y + 0.5)/2)
        }else if(famlist[[j]]$family =="gaussian"){
          eta.old =  Y - off
        }else if(famlist[[j]]$family =="poisson"){
          eta.old = off + 0.1
        }
      }else{
        if(length(ob[[j]])==1){
          if(transpose){
            eta.old <- as.matrix(as.numeric(eta1[,ob[[j]]]), ncol=1)
          }else{
            eta.old <- t(as.matrix(as.numeric(eta1[ob[[j]],]), ncol=1))
          }
        }else{
          if(transpose){
            eta.old <- as.matrix(eta1[,ob[[j]]])
          }else{
            eta.old <- as.matrix(eta1[ob[[j]],])
          }
        }
      }
      z_k <- NULL; w_k <- NULL
      for(l in 1:ncol(Y)){
        etastart <- eta.old[,l]
        weights <- wt
        y <- Y[,l]
        if(famlist[[j]]$family=="binomial"){
          weights <- wt*0+1
          eval(famlist[[j]]$initialize)
        }
        if(famlist[[j]]$family == "poisson"){
          eval(famlist[[j]]$initialize)
        }

        mu = famlist[[j]]$linkinv(eta.old[,l])
        varg = famlist[[j]]$variance(mu)
        gprime = famlist[[j]]$mu.eta(eta.old[,l])
        dev.new = sum(famlist[[j]]$dev.resids(y, mu, wt))

        dev1 <- dev1 + dev.new
        z_k <- cbind(z_k, eta.old[,l] + (y - mu) / gprime  - off[,l])
        w_k <- cbind(w_k, as.matrix(as.vector(wt[1]*gprime^2 / varg)))
      }

      if(transpose){
        z  = cbind(z, z_k)
        W = cbind(W, w_k)
      }else{
        z  =rbind(z, z_k)
        W = rbind(W,  w_k)
      }
      rownum <- c(rownum, ob[[j]])
    }
    time2 <- Sys.time()
    #Final Calculation
    if(sparse==F | predicting==T){ #If sparse==F, do NOT put a penalty on anything
      beta.new <- NULL
      if(transpose){
        for(k in 1:ncol(W)){
          Xtran <- t(predictor)
          temp <- rep(0, ncol(Xtran))
          try(
            temp <- solve(crossprod(Xtran,W[,k]*Xtran),
                          crossprod(Xtran,W[,k]*z[,k]), tol=2*.Machine$double.eps),
            silent=T)
          beta.new <- cbind(beta.new, temp)
        }
        eta.new = Xtran %*% beta.new + t(offsets[rownum,])
        beta.new <- t(beta.new)
      }else{
        for(k in 1:ncol(outcome)){
          Xtran <- predictor[rownum,]
          temp <- dat$beta.new[,k] * 0.01
          try(
            temp <- solve(crossprod(Xtran,W[,k]*Xtran),
                          crossprod(Xtran,W[,k]*z[,k]), tol=2*.Machine$double.eps) ,
            silent=T)
          beta.new = cbind(beta.new, temp)
        }
        eta.new = Xtran %*% beta.new + offsets[rownum,]
      }
    }else if(sparse & score==F){ ###Sparse Loadings
      time3 <- Sys.time()
      y.temp <- NULL
      x.temp <- list()
      for(k in 1:ncol(W)){
        x.temp[[k]] <- sqrt(diag(W[,k])) %*% t(as.matrix(predictor)) #X=W^(1/2)*X
        y.temp <- rbind(y.temp, sqrt(diag(W[,k]))%*%z[,k])
      }
      x.temp2 <- Matrix::bdiag(x.temp)

      r <- ncol(x.temp2)
      #time1 <- Sys.time()
      if(sum(abs(y.temp))==0 | sum(abs(x.temp2))==0){
        beta.sparse <- rep(0,r)
      }else{
        fit.lasso <- glmnet::glmnet(as.matrix(x.temp2),y.temp,
                            family="gaussian",
                            alpha=1, nlambda=1, lambda=lambda, #*r/100,
                            intercept = FALSE, standardize = TRUE,
                            thresh=0.001)
        #time2 <- Sys.time()
        beta.sparse <- as.numeric(fit.lasso$beta)
      }
      time4 <- Sys.time()
      #Incase Rank>1, re-format matrix
      beta.new <- matrix(beta.sparse, ncol=dim(predictor)[1], byrow=F)

      #Part2: Backtracking line search
      fl <- function(x, y=Xtilde, z2=predictor, b.full=dat$beta.new,
                     b.sparse=beta.new, of=offsets[rownum,], wt.vec2=wt.vec,
                     famlist2=famlist, ob2=ob,var01=var0,
                     lambda2=lambda){
        b.full <- as.matrix(b.full)
        b.sparse <- as.matrix(b.sparse)
        Y.pred <-  of + t(x * t(b.full) + (1-x) * t(b.sparse)) %*% z2

        #Calculate Likelihood
        data_ll2 <- NULL
        numvar <- 0
        for(i in dat$dfs){
          X <- y[ob2[[i]],]
          natX <- Y.pred[(numvar+1):(numvar + length(ob2[[i]])),]
          numvar <- numvar + length(ob2[[i]])
          Xfit <- famlist2[[i]]$linkinv(natX)
          n <- ncol(natX)
          if(is.null(n)){n <- 1}
          if(famlist2[[i]]$family=="gaussian"){
            if(is.null(var01)){
              sigma2 <- 1 #var(as.vector(X - Xfit))
            }else{ sigma2 <- var01[i]}
            ll <- -wt.vec2[i] * sum((X - Xfit)^2) / (2 * sigma2)
          }else if(famlist2[[i]]$family=="binomial"){
            ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
          }else if(famlist2[[i]]$family=="poisson"){
            fx <- log(factorial(X))
            fx[which(fx==Inf)] <- -55.22 + 4.32*X[which(fx==Inf)]
            ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
          }
          data_ll2 <- c(data_ll2, -ll)
        }
        return(sum(data_ll2))
      }
      dfl <- function(x){
        if(x==1){h <- -0.001
        }else{h <- 0.001}
        e <- (fl(x+h)-fl(x))/h
        return(e)}
      old.t<-0; new.t <-0; iters <- 0; t.vec <- NULL
      bb <- 0.5; run.backtrack=T; temp <- c(.1,.1)
      #time3 <- Sys.time()
      if(dfl(0)>0){
        new.t=0
        run.backtrack = F
      }
      if(dfl(1)<0){
        new.t=1
        run.backtrack = F
      }
      if(run.backtrack){
        repeat{
          iters <- iters + 1
          old.t <- new.t
          if(abs(dfl(new.t))<0.0000001){break}
          if(new.t == 0 & dfl(new.t)>0){break}
          t.vec <- c(t.vec, new.t)
          temp<-backtrack(x=as.numeric(old.t), dx=sign(dfl(as.numeric(old.t))), f=fl, df=dfl,
                          b.full=dat$beta.new,
                          b.sparse=beta.sparse, beta=bb, t=temp[2])
          new.t <- temp[1]
          if(abs(old.t-new.t)<0.00001){break}
          if(new.t %in% t.vec){bb <- bb/3*2}
          if(iters>500){break}
        }
      }
      bls1 <- iters
      bls.t <- new.t
      new.t <- as.numeric(new.t)
      beta.final <- new.t*as.matrix(dat$beta.new) + (1-new.t) * as.matrix(beta.new)

      # Need to calculate eta.new and beta.new
      eta.new = t(predictor) %*% t(beta.final) + t(offsets[rownum,])
      beta.new <- beta.final
      time5 <- Sys.time()
    }else{ ###Sparse Scores
      time6 <- Sys.time()
      y.temp <- NULL
      x.temp <- list()
      for(k in 1:ncol(W)){
        x.temp[[k]] <- sqrt(diag(W[,k])) %*% as.matrix(predictor) #X=W^(1/2)*X
        y.temp <- rbind(y.temp, sqrt(diag(W[,k]))%*%z[,k])
      }
      x.temp2 <- Matrix::bdiag(x.temp)

      if(sum(abs(as.vector(x.temp2)))!=0){
        # if temp.x== all zeros, then what? force beta=0 and move on
        r <- ncol(x.temp2)
        if(sum(abs(y.temp))==0){
          beta.sparse <- rep(0,r)
        }else{
          # time5 <- Sys.time()
          fit.ridge <- glmnet::glmnet(as.matrix(x.temp2),y.temp,
                              family="gaussian",
                              alpha=0, nlambda=1, lambda=1,
                              intercept = FALSE, standardize = TRUE,
                              thresh=0.001)
          #thresh=10^(-4))
          #time6 <- Sys.time()
          beta.sparse <- as.numeric(fit.ridge$beta)
        }
        time7 <- Sys.time()
        #Incase Rank>1, re-format matrix
        beta.new <- t(matrix(beta.sparse, ncol=dim(predictor)[2], byrow=F))
        #Part2: Backtracking line search
        fs <- function(x, y=Xtilde, z2=predictor, b.full=dat$beta.new,
                       b.sparse=beta.new, of=offsets,wt.vec2=wt.vec,
                       famlist2=famlist,  ob2=ob, var01=var0){
          b.full <- as.matrix(b.full)
          b.sparse <- as.matrix(b.sparse)
          Y.pred <-  of + z2 %*% t(x * t(b.full) + (1-x) * t(b.sparse))

          #Calculate Likelihood
          data_ll2 <- NULL
          for(i in 1:(length(famlist2))){
            X <- y[ob2[[i]],]
            natX <- Y.pred[ob2[[i]],]
            Xfit <- famlist2[[i]]$linkinv(natX)
            n <- ncol(natX)
            if(is.null(n)){n <- 1}
            if(famlist2[[i]]$family=="gaussian"){
              if(is.null(var01)){
                sigma2 <- 1 #var(as.vector(X - Xfit))
              }else{ sigma2 <- var01[i]}
              ll <- -wt.vec2[i] * sum((X - Xfit)^2) / (2 * sigma2)
            }else if(famlist2[[i]]$family=="binomial"){
              ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
            }else if(famlist2[[i]]$family=="poisson"){
              fx <- log(factorial(X))
              fx[which(fx==Inf)] <- -55.22 + 4.32*X[which(fx==Inf)]
              ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
            }
            data_ll2 <- c(data_ll2, -ll)
          }
          return(sum(data_ll2))
        }
        dfs <- function(x){
          if(x==1){h <- -0.001
          }else{h <- 0.001}
          e <- (fs(x+h)-fs(x))/h
          return(e)}
        old.t<-1; new.t <-1; iters <- 0; t.vec <- NULL
        bb <- 0.5; run.backtrack=T; temp <- c(0.1,0.1)
        #time7 <- Sys.time()
        if(dfs(0)>0){
          new.t=0
          run.backtrack = F
        }
        if(dfs(1)<0){
          new.t=1
          run.backtrack = F
        }
        if(run.backtrack){
          repeat{
            iters <- iters + 1
            old.t <- new.t
            if(abs(dfs(new.t))<0.0001){break}
            t.vec <- c(t.vec, new.t)
            temp<-backtrack(x=as.numeric(old.t), dx=sign(dfs(as.numeric(old.t))), f=fs, df=dfs,
                            b.full=dat$beta.new,
                            b.sparse=beta.sparse, beta=bb, t=temp[2])
            new.t <- temp[1]
            if(abs(old.t-new.t)<0.00001){break}
            if(new.t %in% t.vec){bb <- bb/3*2}
            if(iters>500){break}
          }
        }
        bls2 <- iters
        bls.t <- new.t
        new.t <- as.numeric(new.t)
        beta.final <- new.t*as.matrix(dat$beta.new) + (1-new.t) * as.matrix(beta.new)
        #time8 <- Sys.time()
        # Need to calculate eta.new and beta.new
        eta.new = predictor %*% beta.final + offsets[rownum,]
        beta.new <- beta.final
        time8 <- Sys.time()

      }else{
        beta.new <- dat$beta.new * 0
        eta.new = predictor %*% beta.new + offsets[rownum,]
      }
    }

    if(is.null(bls.t)==F){
      dat$bls.t <- c(dat$bls.t, bls.t)
    }
    dat$dev.old <- dat$dev.new
    dat$dev.new <- dev1
    dat$dev.all <- c(dat$dev.all, dev1)
    temp <- dat$dev.new - dat$dev.old
    dat$iter <- dat$iter+1
    #dat$eta.old <- dat$eta.new
    dat$beta.old <- dat$beta.new
    #dat$eta.new <- eta1 <- eta.new
    dat$beta.new <- beta.new
    dat$times <- dat$times + c(as.numeric(as.POSIXct(time2)-as.POSIXct(time1), units="secs"),
                               bls1,
                               as.numeric(as.POSIXct(time4)-as.POSIXct(time3), units="secs"),
                               as.numeric(as.POSIXct(time5)-as.POSIXct(time4), units="secs"),
                               bls2,
                               as.numeric(as.POSIXct(time7)-as.POSIXct(time6), units="secs"),
                               as.numeric(as.POSIXct(time8)-as.POSIXct(time7), units="secs"),
                               as.numeric(as.POSIXct(time1)-as.POSIXct(Sys.time()), units="secs"))


    if(transpose){
      irlslist[[length(irlslist)]][rownum,] <- t(eta.new)
      eta1 <- irlslist[[length(irlslist)]]
    }else{
      irlslist[[length(irlslist)]][rownum,] <- eta.new
      eta1 <- irlslist[[length(irlslist)]]
    }
    irlslist[[list_num]] <- dat

    if(converge & (temp < thresholds)){break}
  }
  return(irlslist)

}


sesJIVE.converge <- function(X, Y, max.iter=2000, threshold = 0.0001,
                             family.x = rep("gaussian",length(X)),
                             family.y = "gaussian",
                             rankJ=1, rankA=rep(1,length(X)),
                             weights=rep(1,length(X)+1),
                             show.message=F, show.error=F,
                             irls_iter=1, intercept=T, var00=NULL,
                             sparse=F, lambda=1, orthogonal=NULL,
                             initial="uninformative"){
  timeA <- Sys.time()
  #NOTES:
  #  -irls_iter=NULL means that we run the IRLS algorithm until convergence
  #             for each parameter, then repeat until likelihood converges.
  #  -Family must be "gaussian", "poisson", or "binomial"
  #  -initial must be "uninformative", "svd", "JIVE"

  if(is.null(orthogonal)){ #If Sparse model, don't enforce orthogonality
    orthogonal <- (sparse == F)
  }

  ##Step 1: Run IRLS Algorithm ##
  diverged<-F
  k <- length(X)
  obs <- list(); temp <- 0
  for(i in 1:k){
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(X[[i]]))
    obs[[i]] <- temp
  }
  obs[[k+1]] <- max(temp)+1
  X.tilde <- NULL
  for(i in 1:k){X.tilde <- rbind(X.tilde, X[[i]]) }
  X.tilde <- rbind(X.tilde, Y)

  y <- nrow(X.tilde)
  n <- ncol(X.tilde)

  #Get initial natural parameter matrix
  Xnat.start <- NULL
  fam.list <- list()
  family.x <- c(family.x, family.y)
  for(i in 1:(k+1)){
    if(family.x[i]=="gaussian"){
      fam.list[[i]] <- stats::gaussian()
      Xnat.start <- rbind(Xnat.start, X.tilde[obs[[i]],])
    }else if(family.x[i]=="binomial"){
      fam.list[[i]] <- stats::binomial()
      Xnat.start <- rbind(Xnat.start,
                          fam.list[[i]]$linkfun((X.tilde[obs[[i]],]+.5)/2))
    }else if(family.x[i]=="poisson"){
      fam.list[[i]] <- stats::poisson()
      Xnat.start <- rbind(Xnat.start,
                          fam.list[[i]]$linkfun((X.tilde[obs[[i]],]+.5)))
    }else{
      print(paste0(family.x[i], " Distribution Does Not Exist"))
      stop()
    }
  }

  # Start with uninformative values
  irls.list <- list()
  for(i in 1:(2*(k+1))){
    if(i==1){p <- "Sj"; dfs <- 1:(k+1)
    e <- X.tilde*0
    beta <- matrix(stats::rnorm(rankJ*n,0,0.1), ncol=n)
    }else if(i==2){p<- "U"; dfs <- 1:(k+1)
    e <- t(X.tilde*0)
    beta <- matrix(rep(0,rankJ*nrow(X.tilde)), ncol=rankJ)
    }else if(i %% 2 == 1){p <- paste0("S",(i-1)/2); dfs <- c(1,2)#c((i-1)/2,(k+1))
    e <- matrix(rep(0.01, (length(obs[[(i-2)/2]])+1)*n), ncol=n)
    beta <- matrix(stats::rnorm(rankA[(i-1)/2]*n,0,1), ncol=n)
    }else{p <- paste0("W", (i-2)/2); dfs <- c((i-2)/2,(k+1))
    e <- t(matrix(rep(0, (length(obs[[(i-2)/2]])+1)*n), ncol=n))
    beta <- matrix(rep(0, rankA[(i-2)/2]*(length(obs[[(i-2)/2]])+1)), ncol=rankA[(i-2)/2])
    }
    irls.list[[i]] <- list(param=p,
                           dfs = dfs,
                           eta.old=e,
                           eta.new=e,
                           beta.old=beta,
                           beta.new=beta,
                           dev.old = 0,
                           dev.new = 0,
                           iter=0,
                           dev.all = 0,
                           times= rep(0,8),
                           bls.t=NULL)
  }
  m <- length(irls.list)+1
  fit <- m+1
  X_temp <- NULL
  for(i in 1:k){X_temp <- c(X_temp, fam.list[[i]]$linkfun(apply(X.tilde[obs[[i]],], 1, mean)))}
  X_temp <- c(X_temp, fam.list[[k+1]]$linkfun(mean(X.tilde[obs[[k+1]],])))
  irls.list[[m]] <- list(param="mu",
                         dfs = 1:(k+1),
                         eta.old = X.tilde*0, eta.new = X.tilde*0,
                         beta.old = X_temp,
                         beta.new = X_temp,
                         dev.old = 0, dev.new = 0,
                         iter = 0, dev.all=0, times=rep(0,8))
  if(intercept==F){
    irls.list[[m]]$beta.new <- as.matrix(irls.list[[m]]$beta.new * 0)
  }
  irls.list[[fit]] <- X.tilde * 0


  ####If informative values, insert them here
  if(initial=="svd"){
    if(intercept){Xnat.start <- Xnat.start - apply(Xnat.start,1,mean)}
    #joint
    if(rankJ > 0){
      X.svd <- svd(Xnat.start, nu=rankJ, nv=rankJ)
      irls.list[[2]]$beta.old <- irls.list[[2]]$beta.new <- as.matrix(X.svd$u)
      if(rankJ==1){
        irls.list[[1]]$beta.old <- irls.list[[1]]$beta.new <- as.matrix(X.svd$d[1] * t(X.svd$v))
      }else{
        irls.list[[1]]$beta.old <- irls.list[[1]]$beta.new <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) }
    }
    #Individual
    for(i in 1:k){
      X.tilde_i <- Xnat.start[c(obs[[i]],obs[[k+1]]),]
      yi <- nrow(X.tilde_i)
      xi <- (X.tilde_i - irls.list[[2]]$beta.new[c(obs[[i]],obs[[k+1]]),] %*% irls.list[[1]]$beta.new) #X-joint
      if(rankJ==0){
        vi <- diag(rep(1,n))
      }else{
        vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v)
      }
      if(rankA[i] > 0){
        X2.svd <- svd(xi %*% vi, nu=rankA[i], nv=rankA[i])
        irls.list[[(i+1)*2]]$beta.old <- irls.list[[(i+1)*2]]$beta.new <- as.matrix(X2.svd$u)
        if(rankA[i]==1){
          irls.list[[(i+1)*2-1]]$beta.old <- irls.list[[(i+1)*2-1]]$beta.new <- as.matrix(X2.svd$d[1] * t(X2.svd$v))
        }else{
          irls.list[[(i+1)*2-1]]$beta.old <- irls.list[[(i+1)*2-1]]$beta.new  <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) }
      }
    }

  }else if(initial=="jive"){
    if(intercept){Xnat.start <- Xnat.start - apply(Xnat.start,1,mean)}
    Xnat.list <- list()
    for(i in 1:k){Xnat.list[[i]] <- Xnat.start[obs[[i]],]}
    jive.fit <- r.jive::jive(Xnat.list, rankJ=rankJ, rankA = rankA,
                     center = F, scale = F, orthIndiv = F,
                     method="given", showProgress=F)
    #joint
    if(rankJ > 0){
      joint <-NULL
      for(i in 1:k){joint <- rbind(joint, jive.fit$joint[[i]])}
      X.svd <- svd(joint, nu=rankJ, nv=rankJ)
      irls.list[[2]]$beta.old <- irls.list[[2]]$beta.new <- rbind(as.matrix(X.svd$u), rep(0,rankJ))
      if(rankJ==1){
        irls.list[[1]]$beta.old <- irls.list[[1]]$beta.new <- as.matrix(X.svd$d[1] * t(X.svd$v))
      }else{
        irls.list[[1]]$beta.old <- irls.list[[1]]$beta.new <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) }
    }
    #Individual
    for(i in 1:k){
      if(rankA[i] > 0){
        X2.svd <- svd(jive.fit$individual[[i]], nu=rankA[i], nv=rankA[i])
        irls.list[[(i+1)*2]]$beta.old <- irls.list[[(i+1)*2]]$beta.new <- rbind(as.matrix(X2.svd$u), rep(0, rankA[i]))
        if(rankA[i]==1){
          irls.list[[(i+1)*2-1]]$beta.old <- irls.list[[(i+1)*2-1]]$beta.new <- as.matrix(X2.svd$d[1] * t(X2.svd$v))
        }else{
          irls.list[[(i+1)*2-1]]$beta.old <- irls.list[[(i+1)*2-1]]$beta.new  <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) }
      }
    }

  }

  if(sparse){
    scld <- T
  }else{scld <- F}
  int <- rep(1,n)
  error.old <- error <- 0
  evec <- err.dat <- NULL
  temp.err1 <- temp.err2 <- -Inf

  if((family.y != "gaussian") | ("poisson" %in% family.x) |
     ("binomial" %in% family.x)){
    score_iter <- NULL
  }else{
    score_iter <- irls_iter
  }

  timeB <- Sys.time()
  time.dat <- c(as.numeric(as.POSIXct(timeB)-as.POSIXct(timeA), units="secs"),
                rep(0, 10))
  for(iter in 1:max.iter){
    timeC <- Sys.time()
    ###Joint Effect
    WS <- NULL; thetaS <- 0
    for(i in 1:k){ #Calculate Outcome
      temp <-  irls.list[[2+2*i]]$beta.new %*% irls.list[[1+2*i]]$beta.new #W_i S_i
      WS <- rbind(WS, as.matrix(temp)[-nrow(temp),])
      thetaS <- thetaS + temp[nrow(temp),]
    }
    out <- X.tilde
    off <- (irls.list[[2]]$beta.new %*% irls.list[[1]]$beta.new) + rbind(WS,thetaS)


    #Optimize Intercept mu
    if(intercept & iter==1){
      irls.list <- irls_func(irls.list,predictor=t(as.matrix(int)),
                             num_iter = score_iter, list_num = m, offsets=off,
                             outcome = out, Xtilde=X.tilde,
                             transpose = T, ob=obs, thresholds=threshold, famlist=fam.list,
                             eta1=t(irls.list[[fit]]),
                             wt.vec=weights, var0=var00)
      #temp.err2 <- optim.error(irls.list, famlist2=fam.list, kk2=k, ob2=obs,
      #                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
      #                         var0=var00, lambda2=lambda)[[1]]
      #if(temp.err2-temp.err1< -10^-8){cat(paste0("iter:", iter, " mu: ", round(temp.err2-temp.err1,4), " \n"))}
      #temp.err1 <- temp.err2

    }
    timeD <- Sys.time()

    #Optimize U
    if(sparse & iter==1){
      #If sparse, set joint loadings to zero for sparsity
      irls.list[[2]]$beta.new <- irls.list[[2]]$beta.new * 0
      irls.list[[2]]$beta.old <- irls.list[[2]]$beta.old * 0
    }
    off <-irls.list[[m]]$beta.new %*% int + rbind(WS,thetaS)
    irls.list <- irls_func(irls.list,predictor=irls.list[[1]]$beta.new, #Sj is known
                           num_iter = score_iter, list_num = 2, offsets=off,
                           outcome = out, transpose=T, thresholds=threshold, famlist=fam.list,
                           ob=obs, eta1=t(irls.list[[fit]]),
                           wt.vec=weights, score=F,
                           sparse=sparse, lambda=lambda,
                           kk=k, Xtilde=X.tilde, old.err=temp.err1,
                           var0=var00)
    sum.U <- sum(abs(irls.list[[2]]$beta.new[-nrow(irls.list[[2]]$beta.new),]))
    if(as.numeric(sum.U)==0){
      irls.list[[2]]$beta.new[nrow(irls.list[[2]]$beta.new),] <- 0
    }
    #temp.err2 <- optim.error(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
    #                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
    #                         var0=var00, lambda2=lambda)[[1]]
    #if(temp.err2-temp.err1< -10^-8){cat(paste0("iter:", iter, " U: ", round(temp.err2-temp.err1,4), " \n"))}
    #temp.err1 <- temp.err2
    #print("Made U")
    timeE <- Sys.time()

    #Optimize Sj
    irls.list <- irls_func(irls.list,predictor=irls.list[[2]]$beta.new, #U is known
                           num_iter = score_iter, list_num = 1, offsets=off,
                           outcome = out, thresholds=threshold, famlist=fam.list,
                           transpose=F, ob=obs, eta1=irls.list[[fit]], Xtilde=X.tilde,
                           wt.vec=weights, scale=scld, score=T, sparse=sparse, var0=var00,
                           full.obs=obs)
    #temp.err2 <- optim.error(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
    #                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
    #                         var0=var00, lambda2=lambda)[[1]]
    #if(temp.err2-temp.err1< -10^-8){print(paste0("iter:", iter, " Sj: ", temp.err2-temp.err1))}
    #temp.err1 <- temp.err2
    #print("Made Sj")
    timeF <- Sys.time()
    timeG <- timeH <- list()
    ###Individual Effect
    for(i in 1:k){
      timeG[[i]] <- Sys.time()
      if(sparse & iter==1){
        irls.list[[(i+1)*2]]$beta.old <- irls.list[[(i+1)*2]]$beta.new <-irls.list[[(i+1)*2]]$beta.new * 0
      }
      #Calculate Outcome and Predictor
      thetaS <-0
      for(j in 1:k){
        if(j != i){
          temp <- irls.list[[2*j+2]]$beta.new %*% irls.list[[2*j+1]]$beta.new
          thetaS <- thetaS + temp[nrow(temp),]
        }
      }
      out <- X.tilde
      off <- irls.list[[m]]$beta.new %*% int + irls.list[[2]]$beta.new %*% irls.list[[1]]$beta.new
      off[nrow(off),] <- off[nrow(off),] + thetaS


      #Optimize Wi
      irls.list <- irls_func(irls.list,predictor=irls.list[[2*i+1]]$beta.new, #Si is known
                             num_iter = score_iter,
                             list_num = (2*i+2), outcome = out,
                             offsets=off,
                             transpose=T,
                             thresholds=threshold, famlist=fam.list,
                             ob=obs, eta1=t(irls.list[[fit]]),
                             wt.vec=weights, sparse=sparse, lambda=lambda,
                             kk=k, Xtilde=X.tilde, old.err=temp.err1,
                             var0=var00)
      sum.W <- sum(abs(irls.list[[2*i+2]]$beta.new[-nrow(irls.list[[2*i+2]]$beta.new),]))
      if(as.numeric(sum.W)==0){
        irls.list[[2*i+2]]$beta.new[nrow(irls.list[[2*i+2]]$beta.new),] <- 0
      }
      #temp.err2 <- optim.error(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
      #                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
      #                         var0=var00, lambda2=lambda)[[1]]
      #if(temp.err2-temp.err1< -10^-8){cat(paste0("iter:", iter, " W", i, ": ",round(temp.err2-temp.err1,4), " \n"))}
      #temp.err1 <- temp.err2
      #print("Made Wi")
      timeH[[i]] <- Sys.time()

      #Optimize Si
      keep.eta <- irls.list[[fit]]
      keep.obs <- c(obs[[i]], obs[[k+1]])
      obs.temp <- list(1:length(obs[[i]]), length(obs[[i]])+1)
      irls.list[[fit]] <- irls.list[[fit]][keep.obs,]
      irls.list <- irls_func(irls.list, predictor=irls.list[[2*i+2]]$beta.new,
                             num_iter = score_iter, list_num = (2*i+1), outcome = out[keep.obs,],
                             offsets=off[keep.obs,], Xtilde=X.tilde[keep.obs,],
                             thresholds=threshold, famlist=list(fam.list[[i]], fam.list[[k+1]]),
                             transpose=F, ob=obs.temp, eta1=irls.list[[fit]],
                             wt.vec=weights[c(i, k+1)], sparse=sparse, score=sparse,
                             var0=var00[c(i,k+1)], full.obs=obs)
      keep.eta[keep.obs,] <- irls.list[[fit]]
      irls.list[[fit]] <- keep.eta
      #temp.err2 <- optim.error(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
      #                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
      #                         var0=var00, lambda2=lambda)[[1]]
      #if(temp.err2-temp.err1< -10^-8){cat(paste0("iter:", iter, " S", i, ": ", round(temp.err2-temp.err1,4), " \n"))}
      #temp.err1 <- temp.err2
      #print("Made Si")
    }
    timeI <- Sys.time()

    time.dat[2:8] <- time.dat[2:8] +
      c(as.numeric(as.POSIXct(timeD)-as.POSIXct(timeC), units="secs"),
        as.numeric(as.POSIXct(timeE)-as.POSIXct(timeD), units="secs"),
        as.numeric(as.POSIXct(timeF)-as.POSIXct(timeE), units="secs"),
        #as.numeric(as.POSIXct(timeG)-as.POSIXct(timeF), units="secs"),
        as.numeric(as.POSIXct(timeH[[1]])-as.POSIXct(timeG[[1]]), units="secs"),
        as.numeric(as.POSIXct(timeG[[2]])-as.POSIXct(timeH[[1]]), units="secs"),
        as.numeric(as.POSIXct(timeH[[2]])-as.POSIXct(timeG[[2]]), units="secs"),
        as.numeric(as.POSIXct(timeI)-as.POSIXct(timeH[[2]]), units="secs"))


    #Figure out the error
    error <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                         Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                         var0=var00, lambda2=lambda)
    err.dat <- rbind(err.dat, error$data_lik)
    evec <- c(evec, error$log_lik)
    error <- error$log_lik
    if(show.error){print(paste0("Iter: ", iter, " Error: ", round(error,4)))}
    if(abs(error.old-error) < threshold){
      #if converged, then stop loop
      if(show.message){cat(paste0("Converged in ", iter, " iterations \n"))}
      break
    }else if(iter == max.iter){
      if(show.message){cat(paste0("Warning: ", iter, " iterations reached \n"))}
      break
    }else if(iter>10){
      elength <- length(evec)
      if(max(evec[(elength-8):elength] - evec[(elength-9):(elength-1)])<0 &
         diverged==F){
        diverged<-T
        if(show.message){cat(paste0("Warning: Diverging at ", iter, " iterations \n"))}
        #  break
      }
    }
    #If didn't converge, prep for another loop
    error.old <- error
  }

  ##Step 2: Save results ##
  U.norm <- W.norm <- 0
  muu <- irls.list[[m]]$beta.new
  Sj <- irls.list[[1]]$beta.new
  U <- irls.list[[2]]$beta.new
  W <- Si <- list()
  for(i in 1:k){
    Si[[i]] <- irls.list[[i*2+1]]$beta.new
    W[[i]] <- irls.list[[i*2+2]]$beta.new
  }

  if(orthogonal==F){
    ##New scaling (7/13/21)
    muu_new <- muu
    U.norm <- norm(U, type="F")
    if(U.norm != 0){
      U_new <- U / U.norm
    }else{
      U_new <- U
    }
    Sj_new <- Sj * U.norm
    W_new <- Si_new <- list()
    W.norm <- list()
    for(i in 1:k){
      W.norm[[i]] <- norm(W[[i]], type="F")
      if(W.norm[[i]] != 0){
        W_new[[i]] <- W[[i]] / norm(W[[i]], type="F")
      }else{
        W_new[[i]] <- W[[i]]
      }
      Si_new[[i]] <- Si[[i]] * norm(W[[i]], type="F")
    }

  }else{
    #Force Orthogonality
    W_new <- Si_new <- list()
    WS_old <- WS_new <- NULL
    thetaS_old <- thetaS_new <- rep(0,n)
    if(rankJ==0){
      vi <- diag(rep(1,n))
    }else{
      X.svd <- svd(U %*% Sj, nu=rankJ, nv=rankJ)
      vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v)
    }
    for(i in 1:k){
        temp.svd <- svd(W[[i]]%*% Si[[i]] %*% vi, nu=rankA[i], nv=rankA[i])
      W_new[[i]] <- temp.svd$u
      Si_new[[i]] <- diag(temp.svd$d[1:rankA[i]],ncol=rankA[i]) %*% t(temp.svd$v)
      WS_old <- rbind(WS_old, W[[i]]%*%Si[[i]])
      WS_new <- rbind(WS_new, W_new[[i]]%*%Si_new[[i]])

      thetaS_old <- thetaS_old + WS_old[nrow(WS_old),]
      thetaS_new <- thetaS_new + WS_new[nrow(WS_new),]
      WS_old <- WS_old[-nrow(WS_old),]
      WS_new <- WS_new[-nrow(WS_new),]
    }
    WS_old <- rbind(WS_old, thetaS_old)
    WS_new <- rbind(WS_new, thetaS_new)

    J_temp <- muu %*% int + U %*% Sj + WS_old - WS_new
    if(intercept){
      muu_new <- apply(J_temp, 1, mean)
    }else{
      muu_new <- muu
    }

      temp.svd <- svd(J_temp-as.matrix(muu_new) %*% int, nu=rankJ, nv=rankJ)

    U_new <- temp.svd$u
    Sj_new <- diag(temp.svd$d[1:rankJ],ncol=rankJ)%*% t(temp.svd$v)
  }

  #Scale so first value in U and W are always positive
  if(U_new[1,1]<0){
    U_new <- -1 * U_new
    Sj_new <- -1 * Sj_new
  }
  for(i in 1:k){
    if(W_new[[i]][1,1]<0){
      W_new[[i]] <- W_new[[i]] * -1
      Si_new[[i]] <- Si_new[[i]] * -1
    }
  }

  #Run times
  run.times <- NULL
  bls.values <- list()
  for(i in 1:(length(irls.list)-1)){
    run.times <- rbind(run.times, c(irls.list[[i]]$param, irls.list[[i]]$iter, irls.list[[i]]$times))
    bls.values[[i]] <- c(irls.list[[i]]$param, irls.list[[i]]$iter, irls.list[[i]]$bls.t)
  }
  run.times <- as.data.frame(run.times)
  names(run.times) <- c("Parameter", "Total IRLS Iters", "IRLS setup time",
                        "BLS1", "glmnet1 time", "BLS1 time",
                        "BLS2", "glmnet2 time", "BLS2 time", "total time")

  #Step 3: Export the results
  U_i <- W_i <- theta_2 <- muu_final <- list()
  natX <- list(); thetaS <- 0
  for(i in 1:k){
    U_i[[i]] <- as.matrix(U_new[obs[[i]],])
    muu_final[[i]] <- as.matrix(muu_new[obs[[i]]])
    W_i[[i]] <- as.matrix(W_new[[i]][-nrow(W_new[[i]]),])
    theta_2[[i]] <- W_new[[i]][nrow(W_new[[i]]),]
    natX[[i]] <- as.matrix(muu_new[obs[[i]]]) %*% int + U_i[[i]] %*% Sj_new + W_i[[i]] %*% Si_new[[i]]
    thetaS <- thetaS + theta_2[[i]] %*% Si_new[[i]]
  }
  theta_1 <-  U_new[obs[[k+1]],]
  natY <- muu_new[obs[[k+1]]] %*% int + theta_1 %*% Sj_new + thetaS
  muu_final[[k+1]] <- muu_new[obs[[k+1]]]

  if(sparse){
    lambdas <- lambda
  }else{
    lambdas <- NULL
  }
  timeJ <- Sys.time()
  time.dat[9] <- as.numeric(as.POSIXct(timeJ)-as.POSIXct(timeA), units="secs")

  result <- list(S_J=Sj_new, S_I=Si_new, U_I=U_i, W_I=W_i,
                 theta1=theta_1, theta2=theta_2, intercept=muu_final,
                 natX=natX, natY=natY,
                 error=error, all.error=evec,
                 iterations = iter, rankJ=rankJ, rankA=rankA,
                 family.x=family.x, family.y=family.y,
                 diverged=diverged, err.dat=err.dat,
                 weights=weights, U.norm=U.norm, W.norm=W.norm,
                 lambda=lambdas, run.time=run.times, time.dat=time.dat,
                 vars=var00, bls.values=bls.values)

  class(result) <- "sesJIVE"
  return(result)
}



find.wts <- function(e, YY, XX, max.iters,
                     folds,
                     family.xx,
                     family.yy, intercepts,
                     rankJJ, rankAA, initials){
  #err.test <- NA
  #for(e in wt.vec){
  err.fold <- NA
  for(i in 1:5){
    #Get train/test sets
    sub.train.x <- sub.test.x <- list()
    sub.train.y <- YY[-folds[[i]]]
    sub.test.y <- YY[folds[[i]]]
    temp.mat <- NULL
    for(j in 1:length(XX)){
      sub.train.x[[j]] <- XX[[j]][,-folds[[i]]]
      sub.test.x[[j]] <- XX[[j]][,folds[[i]]]
      temp.mat <- rbind(temp.mat, e*sub.train.x[[j]])
    }
    temp.mat <-  rbind(temp.mat,
                       (1-e) * sub.train.y)
    temp.norm <- norm(temp.mat, type="F")

    fit1 <- NULL
    attempt <- 0
    while( is.null(fit1) && attempt <= 3 ) {
      attempt <- attempt + 1
      try(
        fit1 <- sesJIVE.converge(sub.train.x, sub.train.y,
                                 max.iter=max.iters, threshold = temp.norm/50000,
                                 family.x = family.xx,
                                 family.y = family.yy,
                                 rankJ=rankJJ, rankA=rankAA,
                                 weights=c(rep(e,length(XX)), 1-e),
                                 show.message=F, show.error=F, initial = initials,
                                 irls_iter=attempt, intercept=intercepts)
      )
    }
    #Record Error for fold
    fit_test1 <- stats::predict(fit1, sub.test.x, show.message=F)
    fit.dev <- get_deviance(sub.test.y, fit_test1$Ynat,
                            family.yy=family.yy)
    err.fold <- c(err.fold, fit.dev)
  }

  #Record Test Error (using validation set)
  fit.dev <- mean(err.fold, na.rm = T)
  fit.dev2 <- sqrt(stats::var(err.fold, na.rm = T))
  return(c(e, fit.dev, fit.dev2))
}


find.lambda <- function(lambda, YY, XX, method="CV", max.iters,
                        folds, weights,
                        family.xx,
                        family.yy, intercepts,
                        rankJJ, rankAA, ortho,
                        var0=NULL, initials){

  if(method=="CV"){
    err.fold <- NA
    for(i in 1:5){
      #Get train/test sets
      sub.train.x <- sub.test.x <- list()
      sub.train.y <- YY[-folds[[i]]]
      sub.test.y <- YY[folds[[i]]]
      for(j in 1:length(XX)){
        sub.train.x[[j]] <- XX[[j]][,-folds[[i]]]
        sub.test.x[[j]] <- XX[[j]][,folds[[i]]]
      }
      fit1 <- NULL
      if(lambda>0){
        attempt <- 0
        while( is.null(fit1) && attempt <= 5 ) {
          attempt <- attempt + 1
          try(
            fit1 <- sesJIVE.converge(sub.train.x, sub.train.y,
                                     max.iter=max.iters, threshold = 0.001,
                                     family.x = family.xx,
                                     family.y = family.yy,
                                     rankJ=rankJJ, rankA=rankAA,
                                     weights=weights, lambda=lambda,
                                     sparse=T, orthogonal=ortho,initial = initials,
                                     show.message=F, show.error=F, var00=var0,
                                     irls_iter=attempt, intercept=intercepts)
          )
        }
      }else{
        attempt <- 0
        while( is.null(fit1) && attempt <= 5 ) {
          attempt <- attempt + 1
          try(
            fit1 <- sesJIVE.converge(sub.train.x, sub.train.y,
                                     max.iter=max.iters, threshold = 0.001,
                                     family.x = family.xx,
                                     family.y = family.yy,
                                     rankJ=rankJJ, rankA=rankAA,
                                     weights=weights,
                                     sparse=F,initial = initials,
                                     show.message=F, show.error=F,
                                     irls_iter=attempt, intercept=intercepts)
          )
        }
        #get var
        var0 <- NULL
        for(i in 1:length(sub.train.x)){
          if(family.xx[i]=="gaussian"){
            var0 <- c(var0, sum((sub.train.x[[i]]-fit1$natX[[i]])^2)/
                        (length(as.vector(fit1$natX[[i]]))-1))
          }else{ var0 <- c(var0, NA)}
        }
        if(family.yy=="gaussian"){
          var0 <- c(var0, sum((sub.train.y-fit1$natY)^2)/
                      (length(as.vector(fit1$natY))-1))
        }else{ var0 <- c(var0, NA)}

      }
      #Record Error for fold
      fit_test1 <- stats::predict(fit1, sub.test.x, show.message = F)
      fit.dev <- get_deviance(sub.test.y, fit_test1$Ynat,
                              family.yy=family.yy)
      err.fold <- c(err.fold, fit.dev)
    }

    #Record Test Error (using validation set)
    fit.dev <- mean(err.fold, na.rm = T)
    fit.se <- sqrt(stats::var(err.fold, na.rm = T)/5)
    lambda2 <- c(lambda, fit.dev, fit.se, var0)
    return(lambda2)

  }else if(method=="BIC"){
    if(lambda > 0){
      fit1 <- sesJIVE.converge(XX, YY,
                               max.iter=max.iters, threshold = 0.001,
                               family.x = family.xx,
                               family.y = family.yy,
                               rankJ=rankJJ, rankA=rankAA,
                               weights=weights, lambda=lambda,
                               sparse=T, orthogonal=ortho,initial = initials,
                               show.message=F, show.error=F, var00=var0,
                               irls_iter=1, intercept=intercepts)
    }else{
      fit1 <- sesJIVE.converge(XX, YY,
                               max.iter=max.iters, threshold = 0.001,
                               family.x = family.xx,
                               family.y = family.yy,
                               rankJ=rankJJ, rankA=rankAA,
                               weights=weights,
                               sparse=F,initial = initials,
                               show.message=F, show.error=F,
                               irls_iter=1, intercept=intercepts)
      #get var
      var0 <- NULL
      for(i in 1:length(XX)){
        if(family.xx[i]=="gaussian"){
          var0 <- c(var0, sum((XX[[i]]-fit1$natX[[i]])^2)/
                      (length(as.vector(fit1$natX[[i]]))-1))
        }else{ var0 <- c(var0, NA)}
      }
      if(family.yy=="gaussian"){
        var0 <- c(var0, sum((YY-fit1$natY)^2)/
                    (length(as.vector(fit1$natY))-1))
      }else{ var0 <- c(var0, NA)}
    }
    ll <- fit1$err.dat[nrow(fit1$err.dat),]
    n <- length(YY)
    ss <- n.total <- 0
    for(i in 1:length(XX)){
      ss <- ss + length(which(as.vector(fit1$U_I[[i]]!=0))) +
        length(which(as.vector(fit1$W_I[[i]]!=0))) +
        length(which(as.vector(fit1$S_I[[i]]!=0)))
      n.total <- n.total + length(as.vector(fit1$U_I[[i]])) +
        length(as.vector(fit1$W_I[[i]])) +
        length(as.vector(fit1$S_I[[i]]))
    }
    ss <- ss + sum(fit1$rankJ, fit1$rankA) + length(which(as.vector(fit1$S_J!=0)))
    n.total <- n.total + sum(fit1$rankJ, fit1$rankA) + length(as.vector(fit1$S_J))
    #Formula based on Meta-analysis based variable selection
    #        by Li et al (2014)
    bic.vec <- -2 * ll + ss * log(n.total)
    lambda2 <- c(lambda, sum(bic.vec), var0)
    return(lambda2)

  }else if(method=="AIC"){
    if(lambda > 0){
      fit1 <- sesJIVE.converge(XX, YY,
                               max.iter=max.iters, threshold = 0.001,
                               family.x = family.xx,
                               family.y = family.yy,
                               rankJ=rankJJ, rankA=rankAA,
                               weights=weights, lambda=lambda,
                               sparse=T, orthogonal=ortho,initial = initials,
                               show.message=F, show.error=F,var00=var0,
                               irls_iter=1, intercept=intercepts)
    }else{
      fit1 <- sesJIVE.converge(XX, YY,
                               max.iter=max.iters, threshold = 0.001,
                               family.x = family.xx,
                               family.y = family.yy,
                               rankJ=rankJJ, rankA=rankAA,
                               weights=weights,
                               sparse=F,initial = initials,
                               show.message=F, show.error=F,
                               irls_iter=1, intercept=intercepts)
      #get var
      var0 <- NULL
      for(i in 1:length(XX)){
        if(family.xx[i]=="gaussian"){
          var0 <- c(var0, sum((XX[[i]]-fit1$natX[[i]])^2)/
                      (length(as.vector(fit1$natX[[i]]))-1))
        }else{ var0 <- c(var0, NA)}
      }
      if(family.yy=="gaussian"){
        var0 <- c(var0, sum((YY-fit1$natY)^2)/
                    (length(as.vector(fit1$natY))-1))
      }else{ var0 <- c(var0, NA)}
    }

    ll <- fit1$err.dat[nrow(fit1$err.dat),]
    ss <- 0
    for(i in 1:length(XX)){
      ss <- ss + length(which(as.vector(fit1$U_I[[i]]!=0))) +
        length(which(as.vector(fit1$W_I[[i]]!=0))) +
        length(which(as.vector(fit1$S_I[[i]]!=0)))
    }
    ss <- ss + sum(fit1$rankJ, fit1$rankA) + length(which(as.vector(fit1$S_J!=0)))
    #Formula based on Meta-analysis based variable selection
    #        by Li et al (2014)
    bic.vec <- -2 * ll + ss * 2
    lambda2 <- c(lambda, sum(bic.vec), var0)
    return(lambda2)


  }
}


sesJIVE.error <- function(Xtilde, U, Sj, W, Si, k, muu, family.x, ob2, kk,
                          wt.vec){
  #Needed for predict.sesJIVE() function
  Y.pred <- NULL
  for(i in 1:k){
    intercept <- as.matrix(muu[[i]]) %*% t(as.matrix(rep(1,ncol(as.matrix(Sj)))))
    J <- as.matrix(U[[i]]) %*% as.matrix(Sj)
    A <- as.matrix(W[[i]]) %*% as.matrix(Si[[i]])
    Y.pred <- rbind(Y.pred, family.x[[i]]$linkinv(intercept + J + A))
  }

  #Calculate Likelihood
  data_ll2 <- NULL
  for(i in 1:(k)){
    X <- Xtilde[ob2[[i]],]
    natX <- Y.pred[ob2[[i]],]
    Xfit <- natX #family.x[[i]]$linkinv(natX)
    Xfit[which(as.numeric(Xfit)>10^12)] <- 10^12
    n <- ncol(natX)
    if(is.null(n)){n <- 1}
    if(family.x[[i]]$family=="gaussian"){
      ll <- -wt.vec[i] * sum((X - Xfit)^2)/2
    }else if(family.x[[i]]$family=="binomial"){
      ll <- wt.vec[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
    }else if(family.x[[i]]$family=="poisson"){
      fx <- as.numeric(log(factorial(X)))
      if(length(which(fx==Inf))>0){
        min.x <- min(X[which(fx==Inf)])
        for(j in 0:(min.x-1)){
          temp <- temp+log(X[which(fx==Inf)]-j)
        }
        fx[which(fx==Inf)] <- temp + log(factorial(X-min.x))
      }
      if(length(which(fx==Inf))>0){
        fx[which(fx==Inf)] <- -55.22 + 4.32*X[which(fx==Inf)]
      }
      ll <- wt.vec[i]*(sum(as.numeric( X*log(Xfit) - Xfit - fx),
                           rm.na=T))
    }
    data_ll2 <- c(data_ll2, ll)
  }
  #############

  return(list(log_lik = sum(data_ll2),
              data_lik = data_ll2))

}


