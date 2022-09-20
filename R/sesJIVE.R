
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
#' a vector of values, the optimal lambda will be chosen based on CV.
#' @param intercept A boolean indicating whether or not there should be an
#' intercept term in the results.
#' @param initial Either "uninformative", "svd", "jive", or "no sparsity" indicating how to
#' generate the initial values for the algorithm. See description for more details.
#' @param show.lambda A boolean indicating if an intermediate table should be printed
#' that shows the predictive performance of each candidate lambda value.
#'
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
                    lambda=NULL,  intercept=T, show.lambda=F,
                    initial="uninformative"){
  ############################################################################
  #X is a list of 2 or more datasets, each with dimensions p_i by n
  #Y is continuous vector length n
  #wts is a tuning parameter between 0 and 1. When wts=NULL, a gridsearch
  #   is conducted to tune eta. You can specify a value of eta to use,
  #   or supply a vector of eta values for esJIVE to consider.
  #rankJ is a value for the low-rank of the joint component
  #rankA is a vector of the ranks for each X dataset.
  ############################################################################

  k <- length(X)
  n <- ncol(X[[1]])
  if(length(Y) != n){stop("Number of columns differ between datasets")}

  #If Sparse model, don't enforce orthogonality
  orthogonal <- (sparse == F)


  if(is.null(wts)){
    wt.vec=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
  }else{
    wt.vec=wts
  }
  if(is.null(lambda)){
    lambda <- c(0, 0.0001, 0.001, 0.01, 0.1, 1)
  }
  lambda.dat <- NULL

  if(length(wt.vec)>1){
    keep.me.wt <- 0.5
  }else{
    keep.me.wt <- wt.vec
  }

  # Choose the best lambda
  if(sparse & length(lambda)>1){
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

    bad.range <- T
    max.lambda <- length(Y)+sum(unlist(lapply(X, function(y) length(as.vector(y)))))
    lambda.vec <- unique(c(0,lambda, 1))
    good.lambdas <- NULL
    while(bad.range==T){
      #with sparsity -- using variance estimate
      doParallel::registerDoParallel(cores=numCores)
      test.best2 <- foreach::foreach(lambda=lambda.vec, .combine=rbind) %dopar% {
        find.lambda(lambda=lambda, YY=Y, XX=X,
                    max.iters=max.iter,
                    folds = fold,initials = "uninformative",
                    weights=c(rep(keep.me.wt,length(X)), 1-keep.me.wt),
                    family.xx = family.x, intercepts=intercept,
                    family.yy = family.y,
                    rankJJ=rankJ, rankAA=rankA)
      }
      doParallel::registerDoParallel(cores=1)
      test.best2 <- as.data.frame(test.best2)
      names(test.best2) <-  c("lambda", "val", "bad_lambda_ind", "pct_sparsity")
      row.names(test.best2) <- c()
      if(length(is.na(test.best2$pct_sparsity))>0){
        test.best2$bad_lambda_ind[which(is.na(test.best2$pct_sparsity))] <- -2
      }
      if(show.lambda){print(test.best2)}

      if(length(which(abs(test.best2$bad_lambda_ind) < 0.5))>3){
        bad.range=F
      }else{
        if(length(which(abs(test.best2$bad_lambda_ind) < 0.5))>0){
          good.lambdas <- rbind(good.lambdas,
                                test.best2[which(abs(test.best2$bad_lambda_ind) < 0.5),])
        }
        too.small <- which(test.best2$bad_lambda_ind < -0.5)
        too.big <- which(test.best2$bad_lambda_ind > 0.5)
        start1 <- test.best2[max(too.small),1]
        stop1 <- test.best2[min(too.big),1]
        lambda.vec <- seq(start1, stop1, by=(stop1-start1)/9)#[2:7]
        cat(paste0("Re-Tuning Parameter with lambda=c(", toString(lambda.vec), ") \n"))
      }
    }
    test.best2 <- rbind(good.lambdas, test.best2)
    if(show.lambda){print(test.best2)}
    good.obs <- which(abs(test.best2$bad_lambda_ind) < 0.5)
    lambda.dat <- test.best2
    temp <- which(test.best2[good.obs,2] == min(test.best2[good.obs,2], na.rm=T))
    best.lambda <- max(test.best2[good.obs[temp],1])
    cat(paste0("Using lambda= ", best.lambda, "\n"))
  }else{
    best.lambda <- lambda
  }

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
    test.best <- NULL
    test.best <- foreach::foreach(e=wt.vec, .combine=rbind) %dopar% {
      try(find.wts(e=e, YY=Y, XX=X,
                   max.iters=max.iter,
                   folds = fold, sparse2=sparse,
                   lambda2=best.lambda,
                   family.xx = family.x, initials = initial,
                   family.yy = family.y, intercepts=intercept,
                   rankJJ=rankJ, rankAA=rankA), silent = T)
    }
    test.best <- as.data.frame(test.best)
    if(show.lambda){print(test.best)}
    doParallel::registerDoParallel(cores=1)
    best.wt <- which(test.best[,2] == min(test.best[,2], na.rm=T))
    best.wt <- best.wt[1]
    keep.me.wt <- as.numeric(as.character(test.best[best.wt,1]))
    cat(paste0("Using wts= ", keep.me.wt, "\n"))
  }else{
    test.best <- t(as.matrix(rep(wt.vec,2)))
    best.wt <- 1
    keep.me.wt <- as.numeric(as.character(wt.vec))
  }


  cat("Estimating Loadings and Scores \n")
  if(initial=="no sparsity"){
    test.first <- sesJIVE.converge(X, Y,
                                   max.iter=500, threshold = threshold,
                                   family.x = family.x,
                                   family.y = family.y,
                                   rankJ=rankJ, rankA=rankA,
                                   weights=c(rep(keep.me.wt,length(X)), 1-keep.me.wt),
                                   show.message=F, show.error=F,
                                   intercept=intercept, sparse=F,
                                   irls_iter=1, lambda=best.lambda,
                                   initial="uninformative",
                                   sesJIVE.fit=NULL)
  }else{test.first=NULL}

  test.best <- sesJIVE.converge(X, Y,
                                max.iter=max.iter, threshold = threshold,
                                family.x = family.x,
                                family.y = family.y,
                                rankJ=rankJ, rankA=rankA,
                                weights=c(rep(keep.me.wt,length(X)), 1-keep.me.wt),
                                show.message=T, show.error=show.error,
                                irls_iter=1, intercept=intercept, sparse=sparse,
                                lambda=best.lambda,
                                initial=initial,
                                sesJIVE.fit=test.first)

  if(test.best$bad.lambda==1 & sparse & length(lambda)>1){
    lambda.red <- as.data.frame(lambda.dat[,1:4])
    row.names(lambda.red) <- c()
    keep.num <- which(abs(lambda.red$bad_lambda_ind)<0.5 & lambda.red$lambda < best.lambda)
    if(length(keep.num)>0){
      new.best.lambda <- lambda.red$lambda[keep.num[which(lambda.red$val[keep.num] == min(lambda.red$val[keep.num], na.rm = T))]][1]
      cat(paste0("Re-estimating Loadings and Scores with lambda= ", new.best.lambda, "\n"))
      test.best <- sesJIVE.converge(X, Y,
                                    max.iter=max.iter, threshold = threshold,
                                    family.x = family.x,
                                    family.y = family.y,
                                    rankJ=rankJ, rankA=rankA,
                                    weights=c(rep(keep.me.wt,length(X)), 1-keep.me.wt),
                                    show.message=T, show.error=show.error,
                                    irls_iter=1, intercept=intercept, sparse=sparse,
                                    lambda=new.best.lambda,
                                    initial=initial,
                                    sesJIVE.fit=test.first)
    }
  }

  test.best$data$X <- X
  test.best$data$Y <- Y

  if(sparse){
    cat("Re-estimating Scores \n")
    test.best.pred <- stats::predict(test.best, X, show.error=show.error,
                              train=T)

    #Combine results into final results
    test.final <- test.best
    test.final$S_J <- test.best.pred$Sj
    test.final$S_I <- test.best.pred$Si
    test.final$pred.all.error <- test.best.pred$all.error

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
  }else{
    test.final <- test.best
  }

  if(sparse & length(lambda)>1){
    lambda.dat <- as.data.frame(lambda.dat[,1:4])
    row.names(lambda.dat) <- c()
    test.final$lambda.dat <- lambda.dat
  }

  dev.resid <- get_deviance(Y, test.final$natY, family.y)

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
#' @param train A boolean for whether or not the predict function is running for the
#'  training dataset. Default is FALSE.
#' @param ... Additional arguments
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
                           show.message=T, train=F, ...){
  ##############################################
  # object is the output from sesJIVE
  # newdata is list with the same predictors and
  #     number of datasets as used in sJIVE.fit
  ##############################################
  U.norm <- object$U.norm
  W.norm <- object$W.norm

  sparse<- FALSE

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
    theta1 <- object$theta1
  }else{
    Sj <- matrix(rep(0,rankJ*n), ncol = n)
    theta1=rep(0, rankJ)
  }

  if(sparse){
    for(i in 1:k){
      U[[i]] <- U.norm * U[[i]]
      W[[i]] <- W.norm[[i]] * W[[i]]
    }
  }

  obs <- rankA <- Si <- theta2 <- list(); temp <- 0; X.tilde <- NULL
  for(i in 1:k){
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(newdata[[i]]))
    obs[[i]] <- temp

    X.tilde <- rbind(X.tilde, newdata[[i]])

    rankA[[i]] <- ncol(as.matrix(W[[i]]))
    if(train){
      Si[[i]] <- object$S_I[[i]]
      theta2[[i]] <- object$theta2[[i]]
    }else{
      Si[[i]] <- matrix(rep(0,rankA[[i]]*n), ncol=n)
      theta2[[i]] <- rep(0, rankA[[i]])
    }
  }
  if(train){
    obs[[k+1]] = max(obs[[k]])+1
    X.tilde <- rbind(X.tilde, object$data$Y)
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
  if(train){ fam.list[[k+1]] = famY}

  #Get Error
  error.old <- sesJIVE.error(X.tilde, U, Sj, W, Si, k, mu, fam.list, ob2=obs, kk=k,
                             wt.vec=weights, train2 = train, theta1 = theta1,
                             theta2 = theta2)

  #Set up IRLS
  k2 <- ifelse(train, k+1, k)
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
  eta.temp <- NULL; WS <- 0
  for(i in 1:k){
    t1 <- mu[[i]] %*% t(as.matrix(rep(1,ncol(as.matrix(Sj))))) +
      U[[i]] %*% matrix(Sj, ncol=n) + W[[i]] %*% matrix(Si[[i]], ncol=n)
    WS <-  WS + matrix(theta2[[i]], ncol=length(theta2[[i]])) %*% matrix(Si[[i]], ncol=n)
    eta.temp <- rbind(eta.temp, t1)
  }
  if(train){
    t2 <- as.numeric(mu[[k+1]]) +  theta1 %*% matrix(Sj, ncol=n) + WS
    eta.temp <- rbind(eta.temp, t2)
  }
  irls.list[[k+2]] <- eta.temp
  fit <- k+2

  ############################ Loop ################
  temp.err1 <- sesJIVE.error(X.tilde, U, Sj, W, Si, k, mu,
                             fam.list,ob2=obs, kk=k, wt.vec=weights,
                             train2 = train, theta1 = theta1,
                             theta2 = theta2)$log_lik

  error.vec <- NULL

  if(("poisson" %in% fam.list) | ("binomial" %in% fam.list)){
    irls_iter <- NULL
  }else{
    irls_iter <- 1
  }

  for(iter in 1:max.iter){

    #Optimize Sj
    U.mat <- A <- mu.mat <- NULL; WS=0
    for(i in 1:k){
      mu.mat <- rbind(mu.mat, as.matrix(mu[[i]]))
      U.mat <- rbind(U.mat, as.matrix(U[[i]]))
      A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
      WS <- matrix(theta2[[i]], ncol=length(theta2[[i]]))  %*% as.matrix(Si[[i]])
    }
    if(train){
      mu.mat <- rbind(mu.mat, as.matrix(mu[[k+1]]))
      U.mat <- rbind(U.mat, as.matrix(theta1))
      A <- rbind(A, WS)
    }
    off <- mu.mat %*% int + A
    beta.old <- irls.list[[1]]$beta.new
    irls.list <- irls_func(irls.list,U.mat, #U is known
                           num_iter = irls_iter, list_num = 1, offsets=off,
                           outcome = X.tilde, thresholds=threshold, famlist=fam.list,
                           transpose=F,  ob=obs, predicting = T,
                           eta1 = irls.list[[k+2]], Xtilde=X.tilde,
                           wt.vec=object$weights)
    temp.err2 <- sesJIVE.error(X.tilde, U, irls.list[[1]]$beta.new, W, Si, k, mu,
                               fam.list,ob2=obs, kk=k, wt.vec=weights,
                               train2 = train, theta1 = theta1,
                               theta2 = theta2)$log_lik

    if(is.na(as.numeric(as.character(temp.err2)))){
      irls.list[[1]]$beta.new <- beta.old
    }else if(as.numeric(as.character(temp.err2))-temp.err1< -1){
      #cat(paste0("Warning: S", i, "wanted to diverge iter ", iter))
      irls.list[[1]]$beta.new <- beta.old
    }else{
      temp.err1 <- as.numeric(as.character(temp.err2))
      Sj <- irls.list[[1]]$beta.new
    }

    for(i in 1:k){
      #Calculate Outcome and Predictor
      W_temp <- matrix(rep(0,rankA[[i]]*nrow(X.tilde)), ncol=rankA[[i]])
      W_temp[obs[[i]],] <- W[[i]]
      off <- mu.mat %*% int + U.mat %*% Sj
      if(train){
        WS <- 0
        for(j in 1:k){
          if(j != i){
            WS <- WS + as.matrix(theta2[[j]]) %*% as.matrix(Si[[j]])
          }
        }
        W_temp[obs[[k+1]],] <- theta2[[i]]
        off[obs[[k+1]],] <- mu[[k+1]] + theta1 %*% Sj + WS
      }

      #Optimize Si
      beta.old <- irls.list[[i+1]]$beta.new
      irls.list <- irls_func(irls.list,predictor=W_temp, #W is known
                             num_iter = irls_iter, list_num = i+1, outcome = X.tilde,
                             offsets=off, Xtilde = X.tilde,
                             thresholds=threshold, famlist=fam.list,
                             transpose=F, ob=obs, predicting = T,
                             eta1 = irls.list[[fit]],
                             wt.vec=object$weights)
      Si.temp <- Si
      Si.temp[[i]] <- irls.list[[i+1]]$beta.new
      temp.err2 <- sesJIVE.error(X.tilde, U, Sj, W, Si.temp, k, mu,
                                 fam.list,ob2=obs, kk=k, wt.vec=weights,
                                 train2 = train, theta1 = theta1,
                                 theta2 = theta2)$log_lik
      if(is.na(as.numeric(as.character(temp.err2)))){
        irls.list[[i+1]]$beta.new <- beta.old
      }else if(as.numeric(as.character(temp.err2))-temp.err1< -1){
        #cat(paste0("Warning: S", i, "wanted to diverge iter ", iter))
        irls.list[[i+1]]$beta.new <- beta.old
      }else{
        temp.err1 <- as.numeric(as.character(temp.err2))
        Si <- Si.temp
      }

    }

    #Get Error
    #Figure out the error
    error.new <- sesJIVE.error(X.tilde, U, Sj, W, Si, k, mu,
                               fam.list,ob2=obs, kk=k, wt.vec=weights,
                               train2 = train, theta1 = theta1,
                               theta2 = theta2)
    if(show.error){print(error.new$log_lik)}
    error.vec <- c(error.vec, error.new$log_lik)
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
              error = error.new$log_lik,
              all.error=error.vec))
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
  k <- length(object$data$X)
  tbl_ranks <- data.frame(Source = c("Joint", paste0("Data", 1:k)),
                          Rank = c(object$rankJ, object$rankA))


  #cat("\n $Variance \n")
  var.table <- NULL; ThetaS <- 0
  for (i in 1:k) {
    j <- object$U_I[[i]] %*% object$S_J
    a <- object$W_I[[i]] %*% object$S_I[[i]]
    ssj <- norm(j, type="f")^2
    ssi <- norm(a, type="f")^2
    sse <- norm(object$data$X[[i]] -j-a, type="f")^2
    sst <- norm(object$data$X[[i]], type="f")^2
    var.table <- cbind(var.table, round(c(ssj/sst, ssi/sst, sse/sst),4))
    ThetaS <- ThetaS + object$theta2[[i]] %*% object$S_I[[i]]
  }
  j <- object$theta1 %*% object$S_J
  ssj <- norm(j, type="f")^2
  ssi <- norm(ThetaS, type="f")^2
  sse <- norm(as.matrix(object$data$Y-j-ThetaS), type="f")^2
  sst <- norm(as.matrix(object$data$Y), type="f")^2
  var.table <- cbind(c("Joint", "Indiv", "Error"),
                     var.table, round(c(ssj/sst, ssi/sst, sse/sst),4))
  var.table <- as.data.frame(var.table)
  names(var.table) <- c("Component", paste0("X", 1:k), "Y")

  #cat("\n $pred.model \n")
  j <- object$theta1 %*% object$S_J
  a <- list(); a2 <- 0
  sse <- sum((object$data$Y - j)^2)
  ssr <- sum((mean(object$data$Y) - j)^2)
  ssnum <-  sum((mean(j) - j)^2)
  coefs <- object$theta1[1]
  for(i in 1:k){
    a[[i]] <- object$theta2[[i]] %*% object$S_I[[i]]
    a2 <- a2 + a[[i]]
    sse <- c(sse, sum((object$data$Y - a[[i]])^2))
    ssr <- c(ssr, sum((mean(object$data$Y) - a[[i]])^2))
    ssnum <- c(ssnum, sum((mean(a[[i]]) - a[[i]])^2))
    coefs <- c(coefs, object$theta2[[i]][1])
  }

  sse_full <- sum((object$data$Y - (j+a2))^2)
  sse_partial <- sum((object$data$Y - a2)^2)
  for(i in 1:k){
    temp <- j
    for(ii in 1:k){
      if(ii != i){temp <- temp + a[[ii]]}
    }
    sse_partial <- c(sse_partial, sum((object$data$Y - temp)^2))
  }
  r_squared <- (sse_partial - sse_full)/sse_partial
  ranks <- c(object$rankJ, object$rankA)

  b <- which(ranks>1)
  n <- length(object$data$Y)
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

  return(list(wts=object$weights, ranks=tbl_ranks,
              variance=var.table,
              pred.model=tbl))
}


#' @describeIn summary.sesJIVE
#'
#' @param x a fitted sesJIVE model.
#' @param ... further arguments passed to or from other methods.
print.sesJIVE <- function(x, ...) {
  k <- length(x$data$X)
  tbl_ranks <- data.frame(Source = c("Joint", paste0("Data", 1:k)),
                          Rank = c(x$rankJ, x$rankA))
  tbl_coef <- NULL
  if(x$rankJ>0){
    for(i in 1:x$rankJ){
      new.col=c(paste0("Joint_",i), x$theta1[i])
      tbl_coef <- cbind(tbl_coef, new.col)
    }}
  for(j in 1:k){
    if(x$rankA[j]>0){
      for(i in 1:x$rankA[j]){
        new.col=c(paste0("Individual",j,"_", i), x$theta2[[j]][i])
        tbl_coef <- cbind(tbl_coef, new.col)
      }}
  }

  cat("weights:", x$weights, "\n")
  cat("Ranks: \n")
  for(i in 1:nrow(tbl_ranks)){
    cat("   ",  unlist(tbl_ranks[i,]), "\n")
  }
  cat("Coefficients: \n")
  for(i in 1:ncol(tbl_coef)){
    cat("   ",  unlist(tbl_coef[,i]), "\n")
  }
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
                         Xtilde2, wt.vec2, sparse,
                         lambda2){
  Sj <- irlslist2[[1]]$beta.new
  U <- irlslist2[[2]]$beta.new
  mu.temp <- as.matrix(irlslist2[[length(irlslist2)-1]]$beta.new) %*% rep(1,ncol(Sj))
  mu.temp2 <- mu.temp[-nrow(mu.temp),]
  WS <- NULL; thetaS <- 0
  penalty <- prod(dim(mu.temp)) * norm(Sj, type="F")^2 +
    prod(dim(mu.temp2))*lambda2 * norm(matrix(matrix(U)[-nrow(matrix(U)),]), type="O")
  for(i in 1:kk2){
    Wi <- irlslist2[[i*2+2]]$beta.new
    Si <- irlslist2[[i*2+1]]$beta.new
    temp <- Wi %*% Si
    WS <- rbind(WS, temp[-nrow(temp),])
    thetaS <- thetaS + temp[nrow(temp),]
    Si_pen <- norm(matrix(as.numeric(Si), ncol=ncol(Si)), type="F")^2
    Wi_pen <- norm(matrix(as.numeric(Wi[-nrow(Wi),]), ncol=ncol(Wi)), type="O")
    penalty <- penalty + prod(dim(temp)) * Si_pen +
      prod(dim(temp))* lambda2 *Wi_pen
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
      ll <- -wt.vec2[i] * sum((X - Xfit)^2) / 2
    }else if(famlist2[[i]]$family=="binomial"){
      ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
    }else if(famlist2[[i]]$family=="poisson"){
      fx <- log(factorial(X))
      high.obs <- which(fx==Inf)
      if(length(high.obs)>0){
        for(j in high.obs){
          temp <- log(factorial(170))
          for(m in 171:X[j]){temp <- temp + log(m)}
          fx[j] <- temp
        }
      }
      ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
    }
    data_ll2 <- c(data_ll2, ll)
  }

  #############
  if(sparse==F){penalty<-0
  }else{data_ll2 <- data_ll2}

  return(list(log_lik = sum(data_ll2)-penalty,
              data_lik = data_ll2))
}

get_deviance <- function(ytrue, ynat, family.yy){
  ytrue <- as.vector(ytrue)
  ynat <- as.vector(ynat)
  n <- length(ytrue)
  if(family.yy=="gaussian"){
    dev.resid <-  2*sum((ytrue - ynat)^2) / 2
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
                      thresholds=0.0001, famlist,
                      outcome, transpose=F, ob, predicting=F,
                      eta1=NULL, wt.vec, sparse=F, lambda=1, kk=k,
                      Xtilde, score=F,
                      full.obs){
  dat <- irlslist[[list_num]]
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
      wt <- rep(wt.vec[j],nobs)
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
        z_k <- cbind(z_k, eta.old[,l] - off[,l] + (y - mu) / gprime )
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
    #Final Calculation
    if(sparse==F){ #If sparse==F, do NOT put a penalty on anything
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
      #part 1: save the unpenalized coefficients for the outcome
      Xtran <- t(predictor)
      temp <- rep(0, ncol(Xtran))
      k <- ncol(W)
      try(
        temp <- solve(crossprod(Xtran,W[,k]*Xtran),
                      crossprod(Xtran,W[,k]*z[,k]), tol=2*.Machine$double.eps),
        silent=T)
      beta.y.coef <- temp

      #Part 2: find the penalized loadings
      y.temp <- NULL
      x.temp <- list()
      for(k in 1:(ncol(W)-1)){
        x.temp[[k]] <- sqrt(diag(W[,k])) %*% t(as.matrix(predictor)) #X.temp=W^(1/2)*X
        y.temp <- rbind(y.temp, sqrt(diag(W[,k]))%*%z[,k])
      }
      x.temp2 <- Matrix::bdiag(x.temp)

      r <- ncol(x.temp2)
      if(sum(abs(y.temp))==0 | sum(abs(x.temp2))==0){
        beta.sparse <- rep(0,r)
      }else{
        fit.lasso <- glmnet::glmnet(as.matrix(x.temp2),y.temp,
                                    family="gaussian",
                                    alpha=1, nlambda=1, lambda=lambda,
                                    intercept = FALSE, standardize = FALSE,
                                    thresh=0.001)
        beta.sparse <- as.numeric(fit.lasso$beta)
      }
      #Incase Rank>1, re-format matrix
      beta.new <- matrix(beta.sparse, ncol=dim(predictor)[1], byrow=F)

      #Part2: Backtracking line search
      fl <- function(x, y=Xtilde[-nrow(Xtilde),], z2=predictor, b.full=dat$beta.new[-nrow(dat$beta.new),],
                     b.sparse=beta.new, of=offsets[rownum[-length(rownum)],], wt.vec2=wt.vec,
                     famlist2=famlist, ob2=ob,
                     lambda2=lambda){
        b.full <- as.matrix(b.full)
        b.sparse <- as.matrix(b.sparse)
        b <- t(x * t(b.full) + (1-x) * t(b.sparse))
        Y.pred <-  of + b %*% z2

        penalty <- prod(dim(Y.pred))*lambda2 * norm(matrix(b), type="O")
        #Calculate Likelihood
        data_ll2 <- NULL
        numvar <- 0
        for(i in dat$dfs[-length(dat$dfs)]){
          X <- y[ob2[[i]],]
          natX <- Y.pred[(numvar+1):(numvar + length(ob2[[i]])),]
          numvar <- numvar + length(ob2[[i]])
          Xfit <- famlist2[[i]]$linkinv(natX)
          n <- ncol(natX)
          if(is.null(n)){n <- 1}
          if(famlist2[[i]]$family=="gaussian"){
            ll <- -wt.vec2[i] * sum((X - Xfit)^2) / 2
          }else if(famlist2[[i]]$family=="binomial"){
            ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
          }else if(famlist2[[i]]$family=="poisson"){
            fx <- log(factorial(X))
            high.obs <- which(fx==Inf)
            if(length(high.obs)>0){
              for(j in high.obs){
                temp <- log(factorial(170))
                for(m in 171:X[j]){temp <- temp + log(m)}
                fx[j] <- temp
              }
            }
            ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
          }
          data_ll2 <- c(data_ll2, -ll)
        }
        return(sum(data_ll2)+penalty)
      }
      dfl <- function(x){
        if(x==1){h <- -0.001
        }else{h <- 0.001}
        e <- (fl(x+h)-fl(x))/h
        return(e)}
      old.t<-0; new.t <-0; iters <- 0; t.vec <- NULL
      bb <- 0.5; run.backtrack=T; temp <- c(.1,.1)
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
          if(new.t == 0 & df(new.t)>0){break}
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
      new.t <- as.numeric(new.t)
      beta.final <- new.t*as.matrix(dat$beta.new[-nrow(dat$beta.new),]) +
        (1-new.t) * as.matrix(beta.new)

      #add back in unpenalized y coefficient
      beta.final <- rbind(beta.final, beta.y.coef)

      # Need to calculate eta.new and beta.new
      eta.new = t(predictor) %*% t(beta.final) + t(offsets[rownum,])
      beta.new <- beta.final
    }else{ ###Sparse Scores
      y.temp <- NULL
      x.temp <- list()
      for(k in 1:ncol(W)){
        x.temp[[k]] <- sqrt(diag(W[,k])) %*% as.matrix(predictor) #X=W^(1/2)*X
        y.temp <- rbind(y.temp, sqrt(diag(W[,k]))%*%z[,k])
      }
      x.temp2 <- Matrix::bdiag(x.temp)

      if(sum(abs(as.vector(x.temp2)))!=0){
        # if temp.x== all zeros, then force beta=0 and move on
        r <- ncol(x.temp2)
        if(sum(abs(y.temp))==0){
          beta.sparse <- rep(0,r)
        }else{
          fit.ridge <- glmnet::glmnet(as.matrix(x.temp2),y.temp,
                                      family="gaussian",
                                      alpha=0, nlambda=1, lambda=1,
                                      intercept = FALSE, standardize = FALSE,
                                      thresh=0.001)
          beta.sparse <- as.numeric(fit.ridge$beta)
        }
        #Incase Rank>1, re-format matrix
        beta.new <- t(matrix(beta.sparse, ncol=dim(predictor)[2], byrow=F))
        #Part2: Backtracking line search
        fs <- function(x, y=Xtilde, z2=predictor, b.full=dat$beta.new,
                       b.sparse=beta.new, of=offsets,wt.vec2=wt.vec,
                       famlist2=famlist,  ob2=ob){
          b.full <- as.matrix(b.full)
          b.sparse <- as.matrix(b.sparse)
          Y.pred <-  of + z2 %*% t(x * t(b.full) + (1-x) * t(b.sparse))

          penalty <- prod(dim(Y.pred)) * norm(matrix(t(x * t(b.full) + (1-x) * t(b.sparse))), type="F")^2
          #Calculate Likelihood
          data_ll2 <- NULL
          for(i in 1:(length(famlist2))){
            X <- y[ob2[[i]],]
            natX <- Y.pred[ob2[[i]],]
            Xfit <- famlist2[[i]]$linkinv(natX)
            n <- ncol(natX)
            if(is.null(n)){n <- 1}
            if(famlist2[[i]]$family=="gaussian"){
              ll <- -wt.vec2[i] * sum((X - Xfit)^2) / 2
            }else if(famlist2[[i]]$family=="binomial"){
              ll <- wt.vec2[i]*(sum( X*log(Xfit) + (1-X)*log(1-Xfit)))
            }else if(famlist2[[i]]$family=="poisson"){
              fx <- log(factorial(X))
              high.obs <- which(fx==Inf)
              if(length(high.obs)>0){
                for(j in high.obs){
                  temp <- log(factorial(170))
                  for(m in 171:X[j]){temp <- temp + log(m)}
                  fx[j] <- temp
                }
              }
              ll <- wt.vec2[i]*(sum( X*log(Xfit) - Xfit - fx))
            }
            data_ll2 <- c(data_ll2, -ll)
          }
          return(sum(data_ll2)+penalty)
        }
        dfs <- function(x){
          if(x==1){h <- -0.001
          }else{h <- 0.001}
          e <- (fs(x+h)-fs(x))/h
          return(e)}
        old.t<-1; new.t <-1; iters <- 0; t.vec <- NULL
        bb <- 0.5; run.backtrack=T; temp <- c(0.1,0.1)
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
        new.t <- as.numeric(new.t)
        beta.final <- new.t*as.matrix(dat$beta.new) + (1-new.t) * as.matrix(beta.new)

        # Need to calculate eta.new and beta.new
        eta.new = predictor %*% beta.final + offsets[rownum,]
        beta.new <- beta.final

      }else{
        beta.new <- dat$beta.new * 0
        eta.new = predictor %*% beta.new + offsets[rownum,]
      }
    }

    dat$dev.old <- dat$dev.new
    dat$dev.new <- dev1
    dat$dev.all <- c(dat$dev.all, dev1)
    dev.change <- as.numeric(as.character(dat$dev.new)) -  as.numeric(as.character(dat$dev.old))
    dat$iter <- dat$iter+1
    dat$beta.old <- dat$beta.new
    dat$beta.new <- beta.new

    if(transpose){
      irlslist[[length(irlslist)]][rownum,] <- t(eta.new)
      eta1 <- irlslist[[length(irlslist)]]
    }else{
      irlslist[[length(irlslist)]][rownum,] <- eta.new
      eta1 <- irlslist[[length(irlslist)]]
    }
    irlslist[[list_num]] <- dat

    if(converge & (dev.change < thresholds)){break}
  }
  return(irlslist)

}


sesJIVE.converge <- function(X, Y, max.iter=2000, threshold = 0.0001,
                             family.x = rep("gaussian",length(X)),
                             family.y = "gaussian",
                             rankJ=1, rankA=rep(1,length(X)),
                             weights=rep(1,length(X)+1),
                             show.message=F, show.error=F,stop.lambda=0,
                             intercept=T, irls_iter=1,
                             sparse=F, lambda=1,
                             initial="uninformative", sesJIVE.fit=NULL){
  #NOTES:
  #  -Family must be "gaussian", "poisson", or "binomial"
  #  -initial must be "uninformative", "svd", "JIVE"
  #set.seed(061821)

  #If Sparse model, don't enforce orthogonality
  orthogonal <- (sparse == F)


  ##Step 1: Run IRLS Algorithm ##
  diverged<-F
  sm.lambda <- lg.lambda <- 0
  min.pct.sparsity = 0.4
  pct.sparsity <- NULL
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
                           dev.all = 0)
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
                         iter = 0, dev.all=0)
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

  }else if(initial=="no sparsity" & is.null(sesJIVE.fit)==F){
    #joint
    if(rankJ>0){
      irls.list[[1]]$beta.old <- irls.list[[1]]$beta.new <- sesJIVE.fit$S_J
    }
    #Individual
    temp.U <- NULL
    for(i in 1:k){
      if(rankJ>0){
        temp.U <- rbind(temp.U, matrix(sesJIVE.fit$U_I[[i]]))
      }
      if(rankA[i] > 0){
        irls.list[[(i+1)*2]]$beta.old <- irls.list[[(i+1)*2]]$beta.new <- rbind(as.matrix(sesJIVE.fit$W_I[[i]]),
                                                                                sesJIVE.fit$theta2[[i]])
        irls.list[[(i+1)*2-1]]$beta.old <- irls.list[[(i+1)*2-1]]$beta.new  <- as.matrix(sesJIVE.fit$S_I[[i]])
      }
    }
    if(rankJ>0){ irls.list[[2]]$beta.old <- irls.list[[2]]$beta.new <- rbind(temp.U, sesJIVE.fit$theta1) }
  }

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

  for(iter in 1:max.iter){
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
                             wt.vec=weights)
    }

    #Optimize U
    temp.err1 <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                              Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                              lambda2=lambda)[[1]]
    evec <- c(evec, temp.err1)
    beta.old <- irls.list[[2]]$beta.new
    if(sparse & iter==1){
      #If sparse, set joint loadings to zero for sparsity
      irls.list[[2]]$beta.new <- irls.list[[2]]$beta.new * 0
      irls.list[[2]]$beta.old <- irls.list[[2]]$beta.old * 0
    }
    off <- matrix(irls.list[[m]]$beta.new) %*% int + rbind(WS,thetaS)
    irls.list <- irls_func(irls.list,predictor=irls.list[[1]]$beta.new, #Sj is known
                           num_iter = score_iter, list_num = 2, offsets=off,
                           outcome = out, transpose=T, thresholds=threshold, famlist=fam.list,
                           ob=obs, eta1=t(irls.list[[fit]]),
                           wt.vec=weights, score=F,
                           sparse=sparse, lambda=lambda,
                           kk=k, Xtilde=X.tilde)
    sum.U <- sum(abs(irls.list[[2]]$beta.new[-nrow(irls.list[[2]]$beta.new),]))
    if(as.numeric(sum.U)==0){
      irls.list[[2]]$beta.new[nrow(irls.list[[2]]$beta.new),] <- 0
    }
    temp.err2 <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                              Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                              lambda2=lambda)[[1]]
    if(is.na(as.numeric(as.character(temp.err2)))){
      irls.list[[2]]$beta.new <- beta.old
    }else if(as.numeric(as.character(temp.err2))-temp.err1< -1 & show.message){
      cat(paste0("Warning: U wanted to diverge iter ", iter))
      irls.list[[2]]$beta.new <- beta.old
    }else{
      temp.err1 <- as.numeric(as.character(temp.err2))
    }
    evec <- c(evec, temp.err1)

    #Optimize Sj
    beta.old <- irls.list[[1]]$beta.new
    irls.list <- irls_func(irls.list,predictor=irls.list[[2]]$beta.new, #U is known
                           num_iter = score_iter, list_num = 1, offsets=off,
                           outcome = out, thresholds=threshold, famlist=fam.list,
                           transpose=F, ob=obs, eta1=irls.list[[fit]], Xtilde=X.tilde,
                           wt.vec=weights, score=T, sparse=sparse,
                           full.obs=obs)
    temp.err2 <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                              Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                              lambda2=lambda)[[1]]
    if(is.na(as.numeric(as.character(temp.err2)))){
      irls.list[[1]]$beta.new <- beta.old
    }else if(as.numeric(as.character(temp.err2))-temp.err1< -1 & show.message){
      cat(paste0("Warning: Sj wanted to diverge iter ", iter))
      irls.list[[1]]$beta.new <- beta.old
    }else{
      temp.err1 <- as.numeric(as.character(temp.err2))
    }
    evec <- c(evec, temp.err1)

    ###Individual Effect
    for(i in 1:k){

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
      off <- matrix(irls.list[[m]]$beta.new)  %*% int + irls.list[[2]]$beta.new %*% irls.list[[1]]$beta.new
      off[nrow(off),] <- off[nrow(off),] + thetaS


      #Optimize Wi
      beta.old <- irls.list[[2*i+2]]$beta.new
      irls.list <- irls_func(irls.list,predictor=irls.list[[2*i+1]]$beta.new, #Si is known
                             num_iter = score_iter,
                             list_num = (2*i+2), outcome = out,
                             offsets=off,
                             transpose=T,
                             thresholds=threshold, famlist=fam.list,
                             ob=obs, eta1=t(irls.list[[fit]]),
                             wt.vec=weights, sparse=sparse, lambda=lambda,
                             kk=k, Xtilde=X.tilde)
      sum.W <- sum(abs(irls.list[[2*i+2]]$beta.new[-nrow(irls.list[[2*i+2]]$beta.new),]))
      if(is.na(as.numeric(sum.W))==F){
        if(as.numeric(sum.W)==0){
          irls.list[[2*i+2]]$beta.new[nrow(irls.list[[2*i+2]]$beta.new),] <- 0
        }}
      temp.err2 <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                                Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                                lambda2=lambda)[[1]]
      #print(temp.err2)
      if(is.na(as.numeric(as.character(temp.err2)))){
        irls.list[[2*i+2]]$beta.new <- beta.old
      }else if(as.numeric(as.character(temp.err2))-temp.err1< -1 & show.message){
        cat(paste0("Warning: W", i, "wanted to diverge iter ", iter))
        irls.list[[2*i+2]]$beta.new <- beta.old
      }else{
        temp.err1 <- as.numeric(as.character(temp.err2))
      }
      evec <- c(evec, temp.err1)

      #Optimize Si
      keep.eta <- irls.list[[fit]]
      keep.obs <- c(obs[[i]], obs[[k+1]])
      obs.temp <- list(1:length(obs[[i]]), length(obs[[i]])+1)
      irls.list[[fit]] <- irls.list[[fit]][keep.obs,]
      beta.old <- irls.list[[2*i+1]]$beta.new
      irls.list <- irls_func(irls.list, predictor=irls.list[[2*i+2]]$beta.new,
                             num_iter = score_iter, list_num = (2*i+1), outcome = out[keep.obs,],
                             offsets=off[keep.obs,], Xtilde=X.tilde[keep.obs,],
                             thresholds=threshold, famlist=list(fam.list[[i]], fam.list[[k+1]]),
                             transpose=F, ob=obs.temp, eta1=irls.list[[fit]],
                             wt.vec=weights[c(i, k+1)], sparse=sparse, score=sparse,
                             full.obs=obs)
      keep.eta[keep.obs,] <- irls.list[[fit]]
      irls.list[[fit]] <- keep.eta
      temp.err2 <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                                Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                                lambda2=lambda)[[1]]
      if(is.na(as.numeric(as.character(temp.err2)))){
        irls.list[[2*i+1]]$beta.new <- beta.old
      }else if(as.numeric(as.character(temp.err2))-temp.err1< -1 & show.message){
        cat(paste0("Warning: S", i, "wanted to diverge iter ", iter))
        irls.list[[2*i+1]]$beta.new <- beta.old
      }else{
        temp.err1 <- as.numeric(as.character(temp.err2))
      }
      evec <- c(evec, temp.err1)
    }

    #Record Sparsity
    if(sparse){
      pct.sparsity <- NULL
      for(i in 2*(1:(k+1))){
        #Across all joint/indiv loadings
        pct.sparsity <- c(pct.sparsity, as.vector(unlist(irls.list[[i]]$beta.new)))
      }
      pct.sparsity <- length(which(abs(pct.sparsity)<0.00000001)) / length(pct.sparsity)
      if(is.na(pct.sparsity)){pct.sparsity=0}
      if(pct.sparsity >= min.pct.sparsity & pct.sparsity != 1){
        sm.lambda <- lg.lambda <- 0
      }
      if(pct.sparsity < min.pct.sparsity){
        sm.lambda <- sm.lambda + 1
        if(stop.lambda>0 & sm.lambda>5){
          sm.lambda <- T; lg.lambda<-F
          if(show.message){cat(paste0("Warning: Lambda=", lambda, " didn't induce enough sparsity \n"))}
          break()}
      }else if(pct.sparsity == 1){
        lg.lambda <- lg.lambda + 1
        if(stop.lambda>0 & lg.lambda>5){
          sm.lambda <- F; lg.lambda<-T
          if(show.message){cat(paste0("Warning: Lambda=", lambda, " caused intercept-only model \n"))}
          break()
        }
      }
    }

    #Figure out the error
    error <- optim.error2(irlslist2 = irls.list, famlist2=fam.list, kk2=k, ob2=obs,
                          Xtilde2=X.tilde, wt.vec2=weights, sparse=sparse,
                          lambda2=lambda)
    err.dat <- rbind(err.dat, error$data_lik)
    evec <- c(evec, error$log_lik, "new iter")
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
    }
    #If didn't converge, prep for another loop
    error.old <- error
  }

  if(sparse){
    if(pct.sparsity == 1){
      sm.lambda <- F
      lg.lambda<-T
      if(show.message){cat(paste0("Warning: Lambda=", lambda, " caused intercept-only model \n"))}
    }else if(pct.sparsity < min.pct.sparsity){
      sm.lambda <- T
      lg.lambda<-F
      if(show.message){cat(paste0("Warning: Lambda=", lambda, " didn't induce enough sparsity \n"))}
    }else{
      sm.lambda <- F
      lg.lambda<-F
    }
  }

  ##Step 2: Save results ##
  U.norm <- 0; W.norm <- list()
  muu <- irls.list[[m]]$beta.new
  Sj <- irls.list[[1]]$beta.new
  U <- irls.list[[2]]$beta.new
  W <- Si <- list()
  for(i in 1:k){
    Si[[i]] <- irls.list[[i*2+1]]$beta.new
    W[[i]] <- irls.list[[i*2+2]]$beta.new
  }

  if(orthogonal==F){
    muu_new <- muu
    U.norm <- norm(U, type="F")
    if(U.norm != 0){
      U_new <- U / U.norm
      Sj_new <- Sj * U.norm
    }else{
      U_new <- U
      Sj_new <- Sj
    }
    W_new <- Si_new <- list()
    W.norm <- list()
    for(i in 1:k){
      W.norm[[i]] <- norm(W[[i]], type="F")
      if(W.norm[[i]] != 0){
        W_new[[i]] <- W[[i]] / norm(W[[i]], type="F")
        Si_new[[i]] <- Si[[i]] * norm(W[[i]], type="F")
      }else{
        W_new[[i]] <- W[[i]]
        Si_new[[i]] <- Si[[i]]
      }
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
    #W.norm[[i]] <- norm(as.matrix(W_i[[i]]), type="F")
  }
  theta_1 <-  U_new[obs[[k+1]],]
  natY <- muu_new[obs[[k+1]]] %*% int + theta_1 %*% Sj_new + thetaS
  muu_final[[k+1]] <- muu_new[obs[[k+1]]]

  if(sparse){
    lambdas <- lambda
  }else{
    lambdas <- NULL
  }

  bad.lambda=0
  if(sm.lambda){
    bad.lambda=-1
  }else if(lg.lambda){
    bad.lambda=1
  }

  output <- list(S_J=Sj_new, S_I=Si_new, U_I=U_i, W_I=W_i,
       theta1=theta_1, theta2=theta_2, intercept=muu_final,
       natX=natX, natY=natY,
       error=error, all.error=evec,
       iterations = iter, rankJ=rankJ, rankA=rankA,
       family.x=family.x, family.y=family.y,
       weights=weights, U.norm=U.norm, W.norm=W.norm,
       lambda=lambdas,
       bad.lambda=bad.lambda, pct.sparsity=pct.sparsity)
  class(output) <- "sesJIVE"
  return(output)
}



find.wts <- function(e, YY, XX, max.iters,
                     folds, sparse2,
                     lambda2,
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

    if("poisson" %in% family.xx){
      temp.scale <- 1000
    }else{
      temp.scale <- 50000
    }

    fit1 <- NULL
    attempt <- 0
    while( is.null(fit1) && attempt <= 3 ) {
      attempt <- attempt + 1
      try(
        fit1 <- sesJIVE.converge(sub.train.x, sub.train.y,
                                 max.iter=max.iters, threshold = temp.norm/temp.scale,
                                 family.x = family.xx,
                                 family.y = family.yy,
                                 sparse=sparse2, lambda=lambda2,
                                 rankJ=rankJJ, rankA=rankAA,
                                 weights=c(rep(e,length(XX)), 1-e),
                                 show.message=F, show.error=F, initial = initials,
                                 irls_iter=attempt, intercept=intercepts)
      )
    }
    if(is.null(fit1)){
      fit.dev <- NA
    }else{
    #Record Error for fold
    fit_test1 <- stats::predict(fit1, sub.test.x, show.message=F)
    fit.dev <- get_deviance(sub.test.y, fit_test1$Ynat,
                            family.yy=family.yy)
    }
    err.fold <- c(err.fold, fit.dev)
  }

  #Record Test Error (using validation set)
  fit.dev <- mean(as.numeric(as.character(err.fold)), na.rm = T)
  fit.dev2 <- sqrt(stats::var(as.numeric(as.character(err.fold)), na.rm = T))
  return(c(e, fit.dev, fit.dev2))
}


find.lambda <-  function(lambda, YY, XX, max.iters,
                         folds, weights,
                         family.xx,
                         family.yy, intercepts,
                         rankJJ, rankAA,
                         initials){

  err.fold <- bad.lamb <- sparsity <- NA
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
                                   family.y = family.yy, stop.lambda = 1,
                                   rankJ=rankJJ, rankA=rankAA,
                                   weights=weights, lambda=lambda,
                                   sparse=T, initial = initials,
                                   show.message=F, show.error=F,
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
    }
    #Record Error for fold
    fit_test1 <- stats::predict(fit1, sub.test.x, show.message = F)
    fit.dev <- get_deviance(sub.test.y, fit_test1$Ynat,
                            family.yy=family.yy)
    err.fold <- c(err.fold, fit.dev)
    bad.lamb <- c(bad.lamb, fit1$bad.lambda)
    sparsity <- c(sparsity, fit1$pct.sparsity)
  }

  #Record Test Error (using validation set)
  fit.dev <- mean(err.fold, na.rm = T)
  fit.se <- sqrt(stats::var(err.fold, na.rm = T)/5)
  bad.lambda <- mean(bad.lamb, na.rm = T)
  pct.sparsity <- mean(sparsity, na.rm = T)
  lambda2 <- c(lambda, fit.dev, bad.lambda, pct.sparsity)
  return(lambda2)
}


sesJIVE.error <- function(Xtilde, U, Sj, W, Si, k, muu, family.x, ob2, kk,
                          wt.vec, train2, theta1=NULL, theta2=NULL){
  #Needed for sesJIVE.predict
  Y.pred <- NULL
  for(i in 1:k){
    intercept <- as.matrix(muu[[i]]) %*% t(as.matrix(rep(1,ncol(as.matrix(Sj)))))
    J <- as.matrix(U[[i]]) %*% as.matrix(Sj)
    A <- as.matrix(W[[i]]) %*% as.matrix(Si[[i]])
    Y.pred <- rbind(Y.pred, family.x[[i]]$linkinv(intercept + J + A))
  }
  if(train2){
    temp <-muu[[k+1]] + matrix(theta1, ncol=length(theta1))  %*% as.matrix(Sj)
    for(i in 1:k){
      temp <- temp + matrix(theta2[[i]], ncol=length(theta2[[i]]))  %*% as.matrix(Si[[i]])
    }
    Y.pred <- rbind(Y.pred, temp)
  }

  #Calculate Likelihood
  k2 <- ifelse(train2, k+1, k)
  data_ll2 <- NULL
  for(i in 1:k2){
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
      fx <- log(factorial(X))
      high.obs <- which(fx==Inf)
      if(length(high.obs)>0){
        for(j in high.obs){
          temp <- log(factorial(170))
          for(m in 171:X[j]){temp <- temp + log(m)}
          fx[j] <- temp
        }
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


