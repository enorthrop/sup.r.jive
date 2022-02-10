#Private helper functions for sJIVE

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
    new.data <- sJIVE.predict(fit.old, test.X)
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
        new.data <- sJIVE.predict(fit.j, test.X)
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
          new.data <- sJIVE.predict(fit.a, test.X)
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
