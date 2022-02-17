

############ sJIVE dataset
##Function
sim.data <- function(k, p, n, rankJ, rankI, prop.causal=NULL,
                     eigval.J=1, eigval.I=1, X.error=.1, Y.error=0.1, n.pred=NULL){
  ######################################################################
  #k=integer, number of datasets
  #p=vector length k, number of predictors for each of the k datasets
  #n=integer, number of observations
  #rankJ=integer, rank for joint decomposition
  #rankI=vector length k, rank of individual decomposition
  #prop.causal=vector length k, proportion of predictors that are causal
  #eigval.J=number, weight of joint signal
  #eigval.I=number, weight of individual signal
  #X.error=number between 0 and 1, proportion of X variance contributed by error
  #Y.error=number between 0 and 1, proportion of Y variance contributed by error
  ######################################################################

  if(is.null(prop.causal)){prop.causal <- rep(1,k)}
  #Simulate U and theta 1
  if(is.null(n.pred)){
    U <- matrix(runif(rankJ, min=0.5, max=1), ncol=rankJ) #theta1 signal
  }else{
    t1 <-runif(n.pred, min=0.5, max=1)
    U <- matrix(c(t1, rep(0,rankJ-n.pred)), ncol=rankJ) #theta1 signal
  }
  for(i in k:1){
    if(prop.causal[i]==1){
      t1 <- matrix(runif(p[i]*rankJ, min=0.5, max=1), ncol=rankJ) #U_i signal
      U <- rbind(t1, U)
    }else if(prop.causal[i]==0){
      t2 <- matrix(rep(0,p[i]*rankJ), ncol = rankJ)
      U <- rbind(t2, U)
    }else{
      pc <- round(p[i]*prop.causal[i])
      t1 <- matrix(runif(pc*rankJ, min=0.5, max=1), ncol=rankJ) #U_i signal
      t2 <- matrix(rep(0,(p[i]-pc)*rankJ), ncol = rankJ)
      U <- rbind(t1, t2, U)
    }
  }
  U <- as.matrix(qr(U)[[1]]) #force orthogonality
  if(as.vector(U)[1] < 0){U <- -1 * U} #Added 3/30/20 for identifiability

  #Simulate S_j
  S_J <-eigval.J * as.matrix(qr(matrix(rnorm(rankJ*n, 0, 1), ncol=n))[[1]])


  #Simulate W_i
  W <- theta2 <- list()
  for(i in 1:k){
    if(prop.causal[i]==1){
      t1 <- matrix(runif(p[i]*rankI[i], min=0.5, max=1), ncol=rankI[i]) #W_i signal
      t2 <- NULL
    }else if(prop.causal[i]==0){
      t1 <- NULL #W_i signal
      t2 <- matrix(rep(0,p[i]*rankI[i]), ncol = rankI[i])
    }else{
      pc <- round(p[i]*prop.causal[i])
      t1 <- matrix(runif(pc*rankI[i], min=0.5, max=1), ncol=rankI[i]) #W_i signal
      t2 <- matrix(rep(0,(p[i]-pc)*rankI[i]), ncol = rankI[i])
    }
    if(is.null(n.pred)){
      t3 <- matrix(runif(rankI[i], min=0.5, max=1), ncol=rankI[i]) #theta2_i signal
    }else{
      tt1 <-runif(n.pred, min=0.5, max=1)
      t3 <- matrix(c(tt1, rep(0,rankI[i]-n.pred)), ncol=rankI[i]) #theta2_i signal
    }
    temp <- qr(rbind(t1, t2, t3))[[1]]
    if(as.vector(temp)[1] < 0){temp <- -1 * temp} #Added 3/30/20 for identifiability
    W[[i]] <- as.matrix(temp[-nrow(temp),])
    theta2[[i]] <- temp[nrow(temp),]
  }

  #Simulate S_i
  S_J.project <- t(S_J) %*% solve(S_J %*% t(S_J)) %*% S_J
  S_i <- list()
  for(i in 1:k){
    temp <-eigval.I * as.matrix(qr(matrix(rnorm(rankI[i]*n, 0, 1), ncol=n))[[1]])
    S_i[[i]] <- temp %*% (diag(rep(1,n)) - S_J.project)
  }

  #Form X matrices
  X <- list(); obs <- 0; error <- list()
  for(i in 1:k){
    error[[i]] <- X[[i]] <- matrix(rep(0, p[i]*n), ncol=n)
    for(j in 1:p[i]){
      temp <- var(as.vector(U[(obs+j),] %*% S_J + W[[i]][j,] %*% S_i[[i]]))
      error[[i]][j,] <- rnorm(n, 0, sqrt(temp*X.error/(1-X.error)))
      X[[i]][j,] <- U[(obs+j),] %*% S_J + W[[i]][j,] %*% S_i[[i]] + error[[i]][j,]
    }
    obs <- obs+p[i]
  }

  #Form Y matrix
  thetaS <- 0
  for(i in 1:k){thetaS <- thetaS + theta2[[i]] %*% S_i[[i]] }
  temp<-var(as.vector(U[nrow(U),] %*% S_J + thetaS))
  err <- matrix(rnorm(n,0, temp*Y.error/(1-Y.error)), ncol=n)
  Y <- U[nrow(U),] %*% S_J + thetaS + err

  #Get error of Unscaled matrix (Added 1/16/20)
  X.temp <- NULL
  for(i in 1:k){X.temp <- rbind(X.temp, X[[i]])}
  obs <- list(); temp <- 0
  obs[[1]] <- 1:p[1]
  for(i in 1:(k-1)){
    temp <- temp + p[i]
    obs[[i+1]] <- (temp + 1):(temp + p[i+1])
  }
  err.old <- optim.error(rbind(X.temp, Y),
                         as.matrix(U[-nrow(U),]), t(as.matrix(U[nrow(U),])), S_J, W,
                         S_i, theta2, k, obs)



  #Scale Result
  obs <- 0
  for(i in 1:k){
    U.temp <- as.matrix(U[(obs+1):(obs+p[i]),])
    for(j in 1:nrow(X[[i]])){
      temp <- var(as.vector(X[[i]][j,]))
      if(temp != 0){
        U.temp[j,] <- U.temp[j,]/sqrt(temp)
        W[[i]][j,] <- W[[i]][j,]/sqrt(temp)
        error[[i]][j,] <- error[[i]][j,]/sqrt(temp)
      }
    }
    U[(obs+1):(obs+p[i]),] <- U.temp
    obs <- obs+p[i]
  }
  #Scaling for Y
  temp <- var(as.vector(Y))
  U[nrow(U),] <- U[nrow(U),]/sqrt(temp)
  thetaS <- 0
  for(i in 1:k){
    theta2[[i]] <- theta2[[i]]/sqrt(temp)
    thetaS <- thetaS + theta2[[i]] %*% S_i[[i]]
  }
  err <- err/sqrt(temp)


  #Keep U and W as othoNORMAL
  for(j in 1:ncol(U)){
    U.norm <-norm(as.matrix(U[,j]),type="F")^2
    U[,j] <- as.matrix(U[,j])/sqrt(U.norm)
    S_J[j,] <- S_J[j,] * sqrt(U.norm)
    obs <- 0
  }
  for(i in 1:k){
    for(j in 1:ncol(W[[i]])){
      W.norm <-norm(rbind(as.matrix(W[[i]][,j]),theta2[[i]][j]),type="F")^2
      W[[i]][,j] <- as.matrix(W[[i]][,j])/sqrt(W.norm)
      theta2[[i]][j] <- theta2[[i]][j]/sqrt(W.norm)
      S_i[[i]][j,] <- S_i[[i]][j,] * sqrt(W.norm)
    }
  }

  #recreate X matrix
  bigXsig <- NULL
  bigXerr <- NULL
  for(i in 1:k){
    X[[i]] <-U[(obs+1):(obs+p[i]),] %*% S_J + W[[i]] %*% S_i[[i]] + error[[i]]
    bigXsig <- rbind(bigXsig, U[(obs+1):(obs+p[i]),] %*% S_J + W[[i]] %*% S_i[[i]] )
    bigXerr <- rbind(bigXerr, error[[i]])
    obs <- obs+p[i]
  }
  eigXsig <- svd(bigXsig)$d[1]
  eigXerr <- svd(bigXerr)$d[1]

  #recreate Y
  Y <- U[nrow(U),] %*% S_J + thetaS + err
  Ysig <- U[nrow(U),] %*% S_J + thetaS
  Ymse <- sum((Ysig-Y)^2)/length(Y)

  return(list(X=X, Y=Y,
              S_J=S_J, S_I=S_i, U_I=U[-nrow(U),], W_I=W,
              theta1=U[nrow(U),], theta2=theta2, errorX = error,
              errorY=err, error.unscl = err.old, Ymse=Ymse,
              Yexact=Ysig, eigenXsignal=eigXsig, eigenXerror=eigXerr))
}

set.seed(021722)
temp <- sim.data(k=2, p=c(30,30), n=20,
                         rankJ=1, rankI=c(1,1),
                         X.error=0.9, Y.error = 0.1)
SimData.norm <- list(c(temp$X, temp$Y))

usethis::use_data(SimData.norm, overwrite = TRUE)
