#sJIVE and sJIVE.predict functions


#' Supervised JIVE (sJIVE)
#'
#' Given multi-source data and a continuous outcome, sJIVE can
#' simultaneously identify shared (joint) and source-specific (individual)
#' underlying structure while building a linear prediction model for an outcome
#' using these structures. These two components are weighted to compromise between
#' explaining variation in the multi-source data and in the outcome.
#'
#' @param X A list of two or more linked data matrices. Each matrix must
#' have the same number of columns, which is assumed to be common, but the
#' number of rows may differ.
#' @param Y A numeric outcome expressed as a vector with length equal
#' to the number of columns in each view of \code{X}.
#' @param rankJ An integer specifying the joint rank of the data.
#' If \code{rankJ=NULL}, ranks will be determined by the \code{method} option.
#' @param rankA A vector specifying the individual ranks of the data.
#' If \code{rankA=NULL}, ranks will be determined by the \code{method} option.
#' @param eta A value or vector of values greater than 0 and less than 1. If \code{eta}
#' is a single value, \code{X} will be weighted by \code{eta} and \code{Y}
#' will be weighted by \code{1-eta}. if \code{eta} is a vector, 5-fold
#' CV will pick the \code{eta} that minimizes the test MSE.
#' @param max.iter The maximum number of iterations for each instance
#' of the sJIVE algorithm.
#' @param threshold The threshold used to determine convergence of the algorithm.
#' @param method A string specifying which rank selection method to use.
#' Possible options are "permute" which uses JIVE's permutation method,
#' or "CV" which uses 5-fold forward CV to determine ranks.
#' @param center.scale A boolean indicating whether or not the
#' data should be centered and scaled.
#' @param reduce.dim A boolean indicating whether or not dimension
#' reduction should be used to increase computation efficiency.
#' @param numCores An integer specifying the number of cores to use when
#' estimating eta. Default is 1.
#'
#' @details The method requires the data to be centered and scaled. This
#' can be done prior to running the method or by specifying center.scale=T.
#' The rank of the joint and individual components as well as
#' the weight between the data and the outcome can be pre-specified
#' or adaptively selected within the function. The method will print the
#' ranks, the weight, and the number of iterations needed
#' to reach convergence.
#'
#' \code{sJIVE} extends \code{jive} to allow for simultaneous prediction
#' of a continuous outcome. It decomposes multi-source data into low-rank,
#' orthogonal joint and individual components. Each component is broken down
#' into the loadings, or left eigenvectors, and the scores, the product of the
#' eigenvalues and the right eigenvectors. The number of eigenvectors is equal to
#' the rank of the component, and the scores are used to predict \code{y}.
#'
#
#' @return \code{sJIVE} returns an object of class "sJIVE". The function \code{summary}
#' (i.e. \code{\link{summary.sJIVE}}) can be used to summarize the model results, including a
#' variance table and testing the significance of the joint and individual components.
#'
#' An object of class "sJIVE" is a list containing the following components.
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
#' @seealso \code{\link{predict.sJIVE}}  \code{\link{summary.sJIVE}}
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20), matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' train.fit <- sJIVE(X=train.x,Y=train.y,rankJ=1,rankA=c(1,1),eta=0.5)
#'
#' \dontrun{
#' train.fit <- sJIVE(X=X,Y=Y,rankJ=1,rankA=c(1,1),eta=0.5)
#' }
sJIVE <- function(X, Y, rankJ = NULL, rankA=NULL,eta=c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99), max.iter=1000,
                  threshold = 0.001,  method="permute",
                  center.scale=TRUE, reduce.dim = TRUE, numCores=1){
  origX <- X
  origY <- Y
  k <- length(X)
  n <- ncol(X[[1]])
  if(length(Y) != n){stop("Number of columns differ between datasets")}

  if(center.scale){
    Y <- as.numeric(scale(as.numeric(Y)))
    for(i in 1:k){
      for(j in 1:nrow(X[[i]])){
        X[[i]][j,] <- as.numeric(scale(as.numeric(X[[i]][j,])))
      }
    }
  }
  e.vec=eta

  if(is.null(rankJ) | is.null(rankA)){
    if(method=="permute"){
      cat("Estimating Ranks via permutation \n")
      temp <- r.jive::jive(X, Y, center=F, scale=F,orthIndiv = F, showProgress = F)
      rankJ <- temp$rankJ
      rankA <- temp$rankA
    }else if(method=="CV"){
      cat("Estimating Ranks via cross-validation \n")
      temp <- sJIVE.ranks(X,Y, eta=eta, max.iter = max.iter, center.scale=center.scale,
                          reduce.dim=reduce.dim)
      rankJ <- temp$rankJ
      rankA <- temp$rankA
    }else{
      errorCondition("Invalid method chosen")
    }
    cat(paste0("Using rankJ= ", rankJ, " and rankA= ", paste(rankA, collapse = " "), "\n"))
  }

  if(length(e.vec)>1){
    #get cv folds
    n <- length(Y)
    fold <- list()
    cutoff <- round(stats::quantile(1:n, c(.2,.4,.6,.8)))
    fold[[1]] <- 1:cutoff[1]
    fold[[2]] <- (cutoff[1]+1):cutoff[2]
    fold[[3]] <- (cutoff[2]+1):cutoff[3]
    fold[[4]] <- (cutoff[3]+1):cutoff[4]
    fold[[5]] <- (cutoff[4]+1):n

    cat("Choosing Tuning Parameter: eta \n")
    err.test <- NA
    doParallel::registerDoParallel(cores=numCores)
    test.best <- foreach::foreach(e=e.vec, .combine=rbind) %dopar% {
      sJIVE.eta(e=e, Y2=Y, X2=X, fold2=fold, max.iter2=max.iter,
                rankJ2=rankJ, rankA2=rankA,
                center.scale2 = center.scale,
                reduce.dim2=reduce.dim)
    }
    doParallel::registerDoParallel(cores=1)
    best.eta <- test.best[which(test.best[,2] == min(test.best[,2], na.rm=T)),1]
    best.eta <- best.eta[1]
    cat(paste0("Using eta= ", best.eta, "\n"))
    test.best <- sJIVE.converge(X, Y, max.iter = max.iter,
                                rankJ = rankJ, rankA = rankA, eta = best.eta,
                                threshold = threshold, center.scale=center.scale,
                                reduce.dim=reduce.dim)
  }else{
    test.best <- sJIVE.converge(X, Y, max.iter = max.iter,
                                rankJ = rankJ, rankA = rankA, eta = e.vec,
                                threshold = threshold, center.scale=center.scale,
                                reduce.dim=reduce.dim)
  }

  return(test.best)

}




#' Prediction for sJIVE
#'
#' Predicted values based on the an sJIVE model.
#'
#' @param object An object of class "sJIVE", usually a fitted sJIVE model.
#' @param newdata A list of matrices representing the new X datasets.
#' @param threshold The threshold used to determine convergence of the algorithm.
#' @param max.iter The maximum number of iterations for each instance
#' of the sJIVE algorithm.
#' @param ... further arguments passed to or from other methods.
#'
#' @details \code{predict.sJIVE} calculates predicted values for \code{newdata}
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
#'
#' @export
#'
#' @examples
#' train.x <- list(matrix(rnorm(300), ncol=20),matrix(rnorm(200), ncol=20))
#' train.y <- rnorm(20)
#' test.x <- list(matrix(rnorm(600), ncol=40),matrix(rnorm(400), ncol=40))
#' train.fit <- sJIVE(X=train.x,Y=train.y,rankJ=1,rankA=c(1,1),eta=0.5)
#' test.fit <- predict(train.fit, newdata = test.x)
predict.sJIVE <- function(object, newdata, threshold = 0.001, max.iter=2000, ...){

  if(object$rankJ==0 & sum(object$rankA)==0){
    return(list(Ypred = 0,
                Sj = 0,
                Si = 0,
                iteration = 0,
                error = NA))
  }

  #Initialize values
  k <- length(newdata)
  n <- ncol(newdata[[1]])
  W <- object$W_I
  U <- object$U_I
  rankJ <- ncol(as.matrix(U[[1]]))
  Sj <- matrix(rep(0,rankJ*n), ncol = n)

  obs <- rankA <- Si <- list(); temp <- 0; X.tilde <- NULL
  for(i in 1:k){
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(newdata[[i]]))
    obs[[i]] <- temp

    X.tilde <- rbind(X.tilde, newdata[[i]])

    rankA[[i]] <- ncol(as.matrix(W[[i]]))
    Si[[i]] <- matrix(rep(0, rankA[[i]]*n), ncol=n)
  }

  #Get Error
  error.old <- sJIVE.pred.err(X.tilde, U, Sj, W, Si, k)

  U.mat <- NULL; pred.J<-NULL
  for(i in 1:k){
    U.mat <- rbind(U.mat, as.matrix(U[[i]]))}
  if(object$rankJ>0){
    pred.J= solve(t(U.mat)%*%U.mat)%*%t(U.mat)}
  pred.A <- list()
  for(i in 1:k){
    if(object$rankA[i]>0)
      pred.A[[i]] <- solve(t(W[[i]])%*%W[[i]])%*%t(W[[i]])
  }

  for(iter in 1:max.iter){

    #Update Sj
    A <- NULL
    for(i in 1:k){
      A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
    }
    if(object$rankJ>0){
      Sj <-pred.J %*% (X.tilde - A)}

    #Update Si
    for(i in 1:k){
      if(object$rankA[i]>0)
        Si[[i]] <- pred.A[[i]] %*% (newdata[[i]] - U[[i]] %*% Sj)
    }

    #Get Error
    error.new <- sJIVE.pred.err(X.tilde, U, Sj, W, Si, k)

    #Check for Convergence
    if(abs(error.old - error.new) < threshold){
      break
    }else{
      error.old <- error.new
    }
  }

  Ypred <- object$theta1 %*% Sj
  for(i in 1:k){
    Ypred <- Ypred + object$theta2[[i]] %*% Si[[i]]
  }

  return(list(Ypred = Ypred,
              Sj = Sj,
              Si = Si,
              iterations = iter,
              error = error.new))
}



#' Summarizing sJIVE Model Fits
#'
#' Summary methods for an sJIVE model of class "sJIVE".
#'
#' @param object An object of class "sJIVE", usually a fitted sJIVE model.
#'
#' @details Both the \code{print} and the \code{summary} functions
#' give summary results for a fitted sJIVE model.
#'
#' For the \code{summary} function, amount of variance explained
#' is expressed in terms of the standardized Frobenious
#' norm. Partial R-squared values are calculated for the
#' joint and individual components. If rank=1, a z-statistic
#' is calculated to determine the p-value. If rank>1, an F-statistic
#' is calculated.
#'
#' For the \code{print} function, the coeffecients are simply
#' printouts of theta1 and theta2 from the sJIVE model.
#'
#' @return Summary measures
#' @export
summary.sJIVE <- function(object, ...){
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
  dat <- data.frame(Y=object$data$Y)
  lm1 <- lm(dat$Y ~ 0)
  rss.old <-sum(residuals(lm1)^2)

  #Joint
  dat <- cbind(dat, t(object$S_J))
  lmfit <- lm(Y ~ 0 + ., data=dat)
  rss.new <- sum(residuals(lmfit)^2)
  aov.dat <- data.frame(Component = "Joint",
                        Df = ncol(dat)-1,
                       `Sum Sq`=rss.old - rss.new)

  #Indiv
  for(i in 1:k){
    rss.old <- rss.new
    dat <- cbind(dat, t(object$S_I[[i]]))
    names(dat)[-1] <- paste0("V", 1:(ncol(dat)-1))
    lmfit <- lm(Y ~ 0 + ., data=dat)
    rss.new <- sum(residuals(lmfit)^2)
    new.row <- c(paste("Indiv", i), ncol(t(object$S_I[[i]])),
                      rss.old - rss.new)
    aov.dat <- rbind(aov.dat, new.row)
  }
  new.row <- c("Residuals", nrow(dat)-sum(as.numeric(aov.dat$Df)), rss.new)
  new.row <- c(new.row, as.numeric(new.row[3])/as.numeric(new.row[2]), "", "")

  aov.dat$`Mean Sq` = as.numeric(aov.dat$Sum.Sq) / as.numeric(aov.dat$Df)
  aov.dat$`F value` = as.numeric(aov.dat$`Mean Sq`) / as.numeric(new.row[4])
  aov.dat$`Pr(>F)` = pf(as.numeric(aov.dat$`F value`), as.numeric(aov.dat$Df),
                        as.numeric(new.row[2]))
  aov.dat <- rbind(aov.dat, new.row)
  for(i in 3:6){
    aov.dat[,i] <- round(as.numeric(aov.dat[,i]), 5)
  }

  return(list(eta=object$eta, ranks=tbl_ranks,
              variance=var.table,
              anova=aov.dat))
}


#' @describeIn summary.sJIVE
#'
#' @param x a fitted sJIVE model.
#' @param ... further arguments passed to or from other methods.
#' @export
print.sJIVE <- function(x, ...) {
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

  cat("eta:", x$eta, "\n")
  cat("Ranks: \n")
  for(i in 1:nrow(tbl_ranks)){
    cat("   ",  unlist(tbl_ranks[i,]), "\n")
  }
  cat("Coefficients: \n")
  for(i in 1:ncol(tbl_coef)){
    cat("   ",  unlist(tbl_coef[,i]), "\n")
  }
}


########################## Helper Functions: #######################
sJIVE.converge <- function(X, Y, eta=NULL, max.iter=1000, threshold = 0.001,
                           show.error =F, rankJ=NULL, rankA=NULL,
                           show.message=T, reduce.dim=T, center.scale=T){
  #X = list(X_1  X_2  ... X_k) with each row centered and scaled
  #Y is continuous vector centered and scaled
  #eta is between 0 and 1, when eta=NULL, no weighting is done
  #rank is prespecified rank of svd approximation
  origX <- X
  origY <- Y

  #Set Ranks
  if(is.null(rankJ) | is.null(rankA)){
    temp <- sJIVE.ranks(X,Y, eta=eta, max.iter = max.iter)
    rankJ <- temp$rankJ
    rankA <- temp$rankA
    cat(paste0("Using rankJ= ", rankJ, " and rankA= ", paste(rankA, collapse = " "), "\n"))
  }
  k <- length(X)

  Y <- as.vector(Y)
  if(center.scale){
    Y <- as.numeric(scale(as.numeric(Y)))
    for(i in 1:k){
      for(j in 1:nrow(X[[i]])){
        X[[i]][j,] <- as.numeric(scale(as.numeric(X[[i]][j,])))
      }
    }
  }

  k <- length(X)
  n <- ncol(X[[1]])
  svd.bigX <- list()
  for(i in 1:k){
    if(ncol(X[[i]]) != n){
      stop("Number of columns differ between datasets")
    }
    if(nrow(X[[i]])>n & reduce.dim){
      svd.bigX[[i]]<- svd(X[[i]], nu=n)
      if(svd.bigX[[i]]$u[1,1]<0){
        svd.bigX[[i]]$u <- svd.bigX[[i]]$u *-1
        svd.bigX[[i]]$v <- svd.bigX[[i]]$v *-1
      }
      X[[i]] <- diag(svd.bigX[[i]]$d) %*% t(svd.bigX[[i]]$v)
    }else{
      svd.bigX[[i]]<-NULL
    }
  }
  svd.bigX[[k+1]] <- 1

  #Step 1: Initialize values
  k <- length(X)
  obs <- list(); temp <- 0
  for(i in 1:k){
    max.obs <- max(temp)
    temp <- (max.obs+1):(max.obs+nrow(X[[i]]))
    obs[[i]] <- temp
  }
  X.tilde <- NULL
  if(is.null(eta)==T){
    for(i in 1:k){X.tilde <- rbind(X.tilde, X[[i]]) }
    X.tilde <- rbind(X.tilde, Y)
  }else{
    for(i in 1:k){X.tilde <- rbind(X.tilde, sqrt(eta) * X[[i]])}
    X.tilde <- rbind(X.tilde, sqrt(1-eta)* Y)
  }
  y <- nrow(X.tilde)
  n <- ncol(X.tilde)

  #Initialize U, theta1, and Sj
  if(rankJ == 0){
    X.svd <- svd(X.tilde, nu=1, nv=1)
    U.new <- U.old <- as.matrix(X.svd$u[-y,]) * 0
    theta1.new <- theta1.old <- t(as.matrix(X.svd$u[y,])) * 0
    Sj.new <- Sj.old <- as.matrix(X.svd$d[1] * t(X.svd$v)) * 0
  }else{
    X.svd <- svd(X.tilde, nu=rankJ, nv=rankJ)
    U.new <- U.old <- as.matrix(X.svd$u[-y,])
    theta1.new <- theta1.old <- t(as.matrix(X.svd$u[y,]))
    if(rankJ==1){Sj.new <- Sj.old <- as.matrix(X.svd$d[1] * t(X.svd$v))
    }else{Sj.new <- Sj.old <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) }
  }

  #Initialize W and S_i = 0
  W.old <- S.old <- theta2.old <- list()
  WS <- NULL
  for(i in 1:k){
    #Get Wi, Si, and theta2i
    X.tilde_i <- X.tilde[c(obs[[i]],y),]
    yi <- nrow(X.tilde_i)
    xi <- (X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new)
    if(rankJ==0){
      vi <- diag(rep(1,n))
    }else{
      vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v)
    }
    if(rankA[i] == 0){
      X2.svd <- svd(xi, nu=1, nv=1)
      W.old[[i]] <- as.matrix(X2.svd$u[-yi,]) *0
      theta2.old[[i]] <- t(as.matrix(X2.svd$u[yi,])) *0
      S.old[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v)) *0
    }else{
      X2.svd <- svd(xi %*% vi, nu=rankA[i], nv=rankA[i])
      W.old[[i]] <- as.matrix(X2.svd$u[-yi,])
      theta2.old[[i]] <- t(as.matrix(X2.svd$u[yi,]))
      if(rankA[i]==1){S.old[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v))
      }else{S.old[[i]] <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) }
    }
    WS <- rbind(WS, W.old[[i]] %*% S.old[[i]])
  }
  error.old <- optim.error(X.tilde, U.old, theta1.old, Sj.old, W.old, S.old, theta2.old, k, obs)
  if(show.error){cat(paste0("Iter: ", 0, "  Error: ", error.old, "\n"))}
  e.vec <- NULL
  S.new <- S.old
  theta2.new <- theta2.old
  W.new <- W.old
  thetaS.sum <- matrix(rep(0,n), ncol=n)

  #Step 2: Interatively optimize problem
  for(iter in 1:max.iter){
    #Get U, theta1, and Sj
    if(rankJ >0){
      X.svd <- svd(X.tilde-rbind(WS, thetaS.sum), nu=rankJ, nv=rankJ)
      neg <- sign(X.svd$u[1,1])
      U.new <- as.matrix(X.svd$u[-y,]) * neg
      theta1.new <- t(as.matrix(X.svd$u[y,])) *neg
      if(rankJ==1){Sj.new <- X.svd$d[1] * t(X.svd$v) * neg
      }else{Sj.new <- diag(X.svd$d[1:rankJ]) %*% t(X.svd$v) * neg }
    }

    if(show.error){
      e <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.old, S.old, theta2.old, k, obs)
      cat(paste0("Iter: ", iter, "  Error: ", e, "  Updated: U, theta1, Sj \n"))
    }

    #Get Wi, Si, and theta2i
    thetaS.sum <- 0; WS <- NULL
    for(i in 1:k){
      thetaS <-0
      for(j in 1:k){if(j != i){ thetaS <- thetaS + theta2.new[[j]] %*% S.new[[j]]}}
      X.tilde_i <- X.tilde[c(obs[[i]],y),]
      yi <- nrow(X.tilde_i)
      X.tilde_i[yi,] <- X.tilde_i[yi,] - thetaS

      if(rankA[i]>0){
        #Get Wi, Si, and theta2i
        xi <- (X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new)
        if(rankJ==0){
          vi <- diag(rep(1,n))
        }else{
          vi <- diag(rep(1,n)) -  X.svd$v %*% t(X.svd$v)
        }
        X2.svd <- svd(xi %*% vi, nu=rankA[i], nv=rankA[i])

        #X2.svd <- svd(X.tilde_i - rbind(as.matrix(U.new[obs[[i]],]), theta1.new) %*% Sj.new, nu=rankA[i], nv=rankA[i])
        neg2 <- sign(X2.svd$u[1,1])
        W.new[[i]] <- as.matrix(X2.svd$u[-yi,]) * neg2
        theta2.new[[i]] <- t(as.matrix(X2.svd$u[yi,])) * neg2
        if(rankA[i]==1){S.new[[i]] <- as.matrix(X2.svd$d[1] * t(X2.svd$v)) * neg2
        }else{S.new[[i]] <- diag(X2.svd$d[1:rankA[i]]) %*% t(X2.svd$v) * neg2}
      }

      #prep for next iteration
      thetaS.sum <- thetaS.sum + theta2.new[[i]] %*% S.new[[i]]
      WS <- rbind(WS, W.new[[i]] %*% S.new[[i]])

      if(show.error){
        e <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.new, S.new, theta2.new, k, obs)
        cat(paste0("Iter: ", iter, "  Error: ", e, "  Updated: Wi, theta2i, Si  i=", i, "\n"))
      }
    }

    #Figure out the error
    error <- optim.error(X.tilde, U.new, theta1.new, Sj.new, W.new, S.new, theta2.new, k, obs)
    if(abs(error.old-error) < threshold){
      #If converged, then stop loop
      if(show.message){cat(paste0("Converged in ", iter, " iterations \n"))}
      break
    }else if(iter == max.iter){
      if(show.message){cat(paste0("Warning: ", iter, " iterations reached \n"))}
      break
    }else{
      #If didn't converge, prep for another loop
      e.vec <- c(e.vec, error)
      U.old <- U.new
      W.old <- W.new
      theta2.old <- theta2.new
      theta1.old <- theta1.new
      Sj.old <- Sj.new
      S.old <- S.new
      error.old <- error
    }
  }

  #Scale so first value in U and W are always positive
  if(U.new[1,1]<0){
    U.new <- -1 * U.new
    theta1.new <- -1 * theta1.new
    Sj.new <- -1 * Sj.new
  }
  for(i in 1:k){
    if(W.new[[i]][1,1]<0){
      W.new[[i]] <- W.new[[i]] * -1
      theta_2[[i]] <- theta_2[[i]] * -1
      S.new[[i]] <- S.new(i) * -1
    }
  }


  #Step 3: Export the results
  U_i <- W <- theta_2 <- list()
  if(is.null(eta)){
    for(i in 1:k){
      U_i[[i]] <- U.new[obs[[i]],]
    }
    theta_1 <- theta1.new
    W <- W.new
    theta_2 <- theta2.new
  }else{
    theta_1 <- (1/sqrt(1-eta)) * theta1.new
    U.new<- (1/sqrt(eta)) * U.new
    for(i in 1:k){
      W[[i]] <- (1/sqrt(eta)) * W.new[[i]]
      theta_2[[i]] <- (1/sqrt(1-eta)) * theta2.new[[i]]
    }
    #ReScale to be norm 1
    for(j in 1:ncol(U.new)){
      U.norm <-norm(rbind(as.matrix(U.new[,j]), theta_1[j]),type="F")^2
      if(U.norm==0){
        U.new[,j] <- as.matrix(U.new[,j])
        theta_1[j] <- theta_1[j]
      }else{
        U.new[,j] <- as.matrix(U.new[,j])/sqrt(U.norm)
        theta_1[j] <- theta_1[j]/sqrt(U.norm)
      }
      Sj.new[j,] <- Sj.new[j,] * sqrt(U.norm)
    }
    for(i in 1:k){
      for(j in 1:ncol(W[[i]])){
        W.norm <-norm(rbind(as.matrix(W[[i]][,j]),theta_2[[i]][j]),type="F")^2
        if(W.norm==0){
          W[[i]][,j] <- as.matrix(W[[i]][,j])
          theta_2[[i]][j] <- theta_2[[i]][j]
        }else{
          W[[i]][,j] <- as.matrix(W[[i]][,j])/sqrt(W.norm)
          theta_2[[i]][j] <- theta_2[[i]][j]/sqrt(W.norm)
        }
        S.new[[i]][j,] <- S.new[[i]][j,] * sqrt(W.norm)
        U_i[[i]] <- U.new[obs[[i]],]
      }
    }

  }

  Ypred <- theta_1 %*% Sj.new
  for(i in 1:k){
    Ypred <- Ypred + theta_2[[i]] %*% S.new[[i]]
  }


  #Map X back to original space
  for(i in 1:k){
    if(is.null(svd.bigX[[i]]) == FALSE){
      U_i[[i]] <- svd.bigX[[i]]$u %*% U_i[[i]]
      W[[i]] <- svd.bigX[[i]]$u %*% W[[i]]
    }
  }

  result <- list(S_J=Sj.new, S_I=S.new, U_I=U_i, W_I=W,
       theta1=theta_1, theta2=theta_2, fittedY=Ypred,
       error=error, all.error=e.vec,
       iterations = iter, rankJ=rankJ, rankA=rankA, eta=eta,
       data=list(X=origX,Y=origY))

  class(result) <- "sJIVE"

  return(result)
}

sJIVE.ranks <- function(X, Y, eta=NULL, max.iter=1000, threshold = 0.01,
                        max.rank=100, center.scale=T,
                        reduce.dim=T){
  cat("Estimating joint and individual ranks via cross-validation... \n")
  k <- length(X)
  n <- ncol(X[[1]])

  #get cv folds
  fold <- list()
  cutoff <- round(stats::quantile(1:n, c(.2,.4,.6,.8)))
  fold[[1]] <- 1:cutoff[1]
  fold[[2]] <- (cutoff[1]+1):cutoff[2]
  fold[[3]] <- (cutoff[2]+1):cutoff[3]
  fold[[4]] <- (cutoff[3]+1):cutoff[4]
  fold[[5]] <- (cutoff[4]+1):n

  rankJ <- 0
  rankA <- rep(0,k)

  #initialize error
  error.old  <- NULL
  for(i in 1:5){
    #get X train, Y.train
    train.X <- test.X <- list()
    for(j in 1:k){
      train.X[[j]] <- X[[j]][,-fold[[i]]]
      test.X[[j]] <- X[[j]][,fold[[i]]]
    }
    train.Y <- Y[-fold[[i]]]
    test.Y <- Y[fold[[i]]]
    fit.old <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter,
                              rankJ = rankJ, rankA = rankA, show.message = F,
                              center.scale=center.scale,
                              reduce.dim=reduce.dim)
    new.data <- stats::predict(fit.old, test.X)
    error.old <- c(error.old, sum((new.data$Ypred-test.Y)^2) )
  }
  err.old <- mean(error.old)

  #iteravely add ranks
  for(iter in 1:1000){
    error.j <- NULL; error.a <- list()

    for(i in 1:5){
      #get X train, Y.train
      train.X <- test.X <- list()
      for(j in 1:k){
        train.X[[j]] <- X[[j]][,-fold[[i]]]
        test.X[[j]] <- X[[j]][,fold[[i]]]
      }
      train.Y <- Y[-fold[[i]]]
      test.Y <- Y[fold[[i]]]

      #Add rank to joint
      if(rankJ < max.rank){
        fit.j <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter,
                                rankJ = rankJ+1, rankA = rankA, show.message=F,
                                center.scale=center.scale,
                                reduce.dim=reduce.dim)
        new.data <- stats::predict(fit.j, test.X)
        error.j <- c(error.j, sum((new.data$Ypred-test.Y)^2) )
      }else{
        error.j <- c(error.j, 99999999)
      }

      #Add rank to individual
      for(j in 1:k){
        if(rankA[j] < max.rank){
          rankA.new <- rankA
          rankA.new[j] <- rankA.new[j]+1

          #Add rank to individual
          fit.a <- sJIVE.converge(train.X, train.Y, eta=eta, max.iter = max.iter,
                                  rankJ = rankJ, rankA = rankA.new, show.message=F,
                                  center.scale=center.scale,
                                  reduce.dim=reduce.dim)
          new.data <- stats::predict(fit.a, test.X)
          if(length(error.a)<j){error.a[[j]] <- NA}
          error.a[[j]] <- c(error.a[[j]], sum((new.data$Ypred-test.Y)^2) )
        }else{
          if(length(error.a)<j){error.a[[j]] <- NA}
          error.a[[j]] <- c(error.a[[j]], 99999999 )
        }
      }

    }

    #average over folds
    err.j <- mean(error.j, na.rm = T)
    err.a <- lapply(error.a, function(x) mean(x, na.rm=T))
    #cat(error.j)
    #cat(error.a[[1]])
    #cat(error.a[[2]])

    #Determine which rank to increase
    asd <- c(err.old - err.j, err.old - as.vector(unlist(err.a)))
    if(max(asd) < threshold){
      break
    }else{
      temp <- which(asd == max(asd))
      if(temp==1){
        #cat(asd)
        rankJ <- rankJ+1
        err.old <- err.j
      }else{
        rankA[temp-1] <- rankA[temp-1]+1
        err.old <- err.a[[temp-1]]
      }
    }
    #cat(paste0("Joint Rank: ", rankJ, "\n"))
    #cat(paste0("Individual Ranks: ", rankA, "\n"))
  }

  return(list(rankJ = rankJ,
              rankA = rankA,
              error = err.old,
              error.joint = err.j,
              error.individual = err.a))
}

optim.error <- function(X.tilde, U, theta1, Sj, W, Si, theta2, k, obs){
  WS <- NULL; thetaS <- 0
  for(i in 1:k){
    temp <- W[[i]] %*% Si[[i]]
    WS <- rbind(WS, temp)
    thetaS <- thetaS + theta2[[i]] %*% Si[[i]]
  }
  WS.new <- rbind(WS, thetaS)
  error  <- norm(X.tilde - rbind(U, theta1) %*% Sj - WS.new, type = "F")^2
  return(error)
}


sJIVE.pred.err <- function(X.tilde, U, Sj, W, Si, k){
  J <- A <- NULL
  for(i in 1:k){
    J <- rbind(J, as.matrix(U[[i]]) %*% as.matrix(Sj))
    A <- rbind(A, as.matrix(W[[i]]) %*% as.matrix(Si[[i]]))
  }
  temp <- X.tilde - J - A
  error <- norm(temp, type = "F")^2
  return(error)
}

sJIVE.eta <- function(e, Y2, X2, fold2, max.iter2,
                      rankJ2, rankA2,
                      center.scale2,
                      reduce.dim2){
    err.fold <- NA
    for(i in 1:5){
      #Get train/test sets
      sub.train.x <- sub.test.x <- list()
      sub.train.y <- Y2[-fold2[[i]]]
      sub.test.y <- Y2[fold2[[i]]]
      for(j in 1:length(X2)){
        sub.train.x[[j]] <- X2[[j]][,-fold2[[i]]]
        sub.test.x[[j]] <- X2[[j]][,fold2[[i]]]
      }
      temp.mat <- rbind(e * sub.train.x[[1]],
                        e * sub.train.x[[2]],
                        (1-e) * sub.train.y)
      temp.norm <- norm(temp.mat, type="F")
      fit1 <- sJIVE.converge(sub.train.x, sub.train.y, max.iter = max.iter2,
                             rankJ = rankJ2, rankA = rankA2, eta = e,
                             show.message=F, center.scale=center.scale2,
                             reduce.dim=reduce.dim2, threshold = temp.norm/50000)
      #Record Error for fold
      fit_test1 <- stats::predict(fit1, sub.test.x)
      if(sum(is.na(fit_test1$Ypred))==0){
        fit.mse <- sum((fit_test1$Ypred - sub.test.y)^2)/length(sub.test.y)
        err.fold <- c(err.fold, fit.mse)
      }else{
        err.fold <- c(err.fold, NA)
      }
    }

    #Record Test Error (using validation set)
    fit.mse <- mean(err.fold, na.rm = T)
    return(c(e, fit.mse))
}
