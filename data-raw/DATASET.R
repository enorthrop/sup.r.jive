

############ sJIVE dataset
##Function
sim.data <- function(k, p, n, rankJ, rankI, prop.causal=NULL,
                     eigval.J=1, eigval.I=1, X.error=.1, Y.error=0.1){
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
  U <- matrix(runif(rankJ, min=0.5, max=1), ncol=rankJ) #theta1 signal
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
    t3 <- matrix(runif(rankI[i], min=0.5, max=1), ncol=rankI[i]) #theta2_i signal
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

set.seed(031122)
dat <- sim.data(k=2, p=c(40,40),n=30,
                         rankJ = 1, rankI=c(1,1),
                         X.error = c(0.9,0.9), Y.error = 0.1)
#fit <- sJIVE(X=dat$X, Y=dat$Y)

SimData.norm <- list()
SimData.norm$X <- dat$X
SimData.norm$Y <- dat$Y

usethis::use_data(SimData.norm, overwrite = TRUE)



############ sesJIVE dataset
##Function
sim.data <- function(k, p, n, rankJ, rankI, prop.causal=NULL,
                     eigval.J=1, eigval.I=1, x.family=rep("gaussian",k),
                     y.family="gaussian", n.pred=NULL, x.mean=rep(NULL,k),
                     y.mean=NULL, scale.x = rep(.1,k), scale.y=.1, var.y=1,
                     var.x =  rep(1,k), orthogonal=T, strength=1){
  ######################################################################
  #k=integer, number of datasets
  #p=vector length k, number of predictors for each of the k datasets
  #n=integer, number of observations
  #rankJ=integer, rank for joint decomposition
  #rankI=vector length k, rank of individual decomposition
  #prop.causal=vector length k, proportion of predictors that are causal
  #eigval.J=number, weight of joint signal
  #eigval.I=number, weight of individual signal
  ######################################################################
  x.mean2 <- rep(0,k)
  if(is.null(x.mean)){x.mean[k+1] <- 0}
  for(i in 1:k){
    if(x.family[i]=="gaussian"){
      fam <- gaussian()
      if(is.na(x.mean[i])){x.mean[i] <- 0}
    }else if(x.family[i]=="binomial"){
      fam <- binomial()
      if(is.na(x.mean[i])){x.mean[i] <- 0.5}
    }else if(x.family[i]=="poisson"){
      fam <- poisson()
      if(is.na(x.mean[i])){x.mean[i] <- 1}
    }else{
      stop()
    }
    x.mean2[i] <- fam$linkfun(x.mean[i])
  }
  if(y.family=="gaussian"){
    fam <- gaussian()
    if(is.null(y.mean)){y.mean <- 0}
  }else if(y.family=="binomial"){
    fam <- binomial()
    if(is.null(y.mean)){y.mean <- 0.5}
  }else if(y.family=="poisson"){
    fam <- poisson()
    if(is.null(y.mean)){y.mean <- 1}
  }else{
    stop()
  }
  y.mean2 <- fam$linkfun(y.mean)

  if(is.null(prop.causal)){prop.causal <- rep(1,k)}
  #Simulate U and theta 1
  if(is.null(n.pred)){
    U <- matrix(runif(rankJ, min=-1, max=1), ncol=rankJ) #theta1 signal
  }else{
    t1 <-runif(n.pred, min=-1, max=1)
    U <- matrix(c(t1, rep(0,rankJ-n.pred)), ncol=rankJ) #theta1 signal
  }
  for(i in k:1){
    if(prop.causal[i]==1){
      t1 <- matrix(runif(p[i]*rankJ, min=-1, max=1), ncol=rankJ) #U_i signal
      U <- rbind(t1, U)
    }else if(prop.causal[i]==0){
      t2 <- matrix(rep(0,p[i]*rankJ), ncol = rankJ)
      U <- rbind(t2, U)
    }else{
      pc <- round(p[i]*prop.causal[i])
      t1 <- matrix(runif(pc*rankJ, min=-1, max=1), ncol=rankJ) #U_i signal
      t2 <- matrix(rep(0,(p[i]-pc)*rankJ), ncol = rankJ)
      U <- rbind(t1, t2, U)
    }
  }
  U <- U* strength
  if(orthogonal){
    U <- as.matrix(qr(U)[[1]]) #force orthogonality
  }
  if(as.vector(U)[1] < 0){U <- -1 * U} #Added 3/30/20 for identifiability

  #Simulate S_j
  S_J <-eigval.J * as.matrix(qr(matrix(rnorm(rankJ*n, 0, 1), ncol=n))[[1]])


  #Simulate W_i
  W <- theta2 <- list()
  for(i in 1:k){
    if(prop.causal[i]==1){
      t1 <- matrix(runif(p[i]*rankI[i], min=-1, max=1), ncol=rankI[i]) #W_i signal
      t2 <- NULL
    }else if(prop.causal[i]==0){
      t1 <- NULL #W_i signal
      t2 <- matrix(rep(0,p[i]*rankI[i]), ncol = rankI[i])
    }else{
      pc <- round(p[i]*prop.causal[i])
      t1 <- matrix(runif(pc*rankI[i], min=-1, max=1), ncol=rankI[i]) #W_i signal
      t2 <- matrix(rep(0,(p[i]-pc)*rankI[i]), ncol = rankI[i])
    }
    if(is.null(n.pred)){
      t3 <- matrix(runif(rankI[i], min=-1, max=1), ncol=rankI[i]) #theta2_i signal
    }else{
      tt1 <-runif(n.pred, min=-1, max=1)
      t3 <- matrix(c(tt1, rep(0,rankI[i]-n.pred)), ncol=rankI[i]) #theta2_i signal
    }
    if(orthogonal){
      temp <- qr(strength * rbind(t1, t2, t3))[[1]]
    }else{temp <- strength * rbind(t1,t2,t3)}
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
  muu <- rep(0, nrow(U)); ones <- t(as.matrix(rep(1,n)))
  for(i in 1:k){
    X[[i]] <- matrix(rep(0, p[i]*n), ncol=n)
    for(j in 1:p[i]){
      temp  <- U[(obs+j),] %*% S_J + W[[i]][j,] %*% S_i[[i]]
      scl.num <- sqrt(var(as.vector(temp))) / sqrt(var.x[i])
      if(scl.num > 0){
        U[(obs+j),] <-  U[(obs+j),]/scl.num
        W[[i]][j,]  <- W[[i]][j,] /scl.num
      }
      temp  <- U[(obs+j),] %*% S_J + W[[i]][j,] %*% S_i[[i]]
      muu[(obs+j)] <- -1*mean(temp) + x.mean2[i]
      X[[i]][j,] <- muu[(obs+j)] %*% ones + temp
    }
    obs <- obs+p[i]
  }

  #Form Y matrix
  thetaS <- 0
  for(i in 1:k){thetaS <- thetaS + theta2[[i]] %*% S_i[[i]] }
  temp <- U[nrow(U),] %*% S_J + thetaS #+ err
  scl.num <- sqrt(var(as.vector(temp))) / sqrt(var.y)
  thetaS <- 0
  if(scl.num > 0){
    U[nrow(U),]  <-  U[nrow(U),] /scl.num
    for(i in 1:k){theta2[[i]] <- theta2[[i]]/scl.num
    thetaS <- thetaS + theta2[[i]] %*% S_i[[i]]
    }}
  temp <- U[nrow(U),] %*% S_J + thetaS
  muu[length(muu)] <- -1*mean(temp) + y.mean2
  Y <- muu[length(muu)] %*% ones + temp


  #Scale Result
  obs <- 0
  for(i in 1:k){
    U.temp <- as.matrix(U[(obs+1):(obs+p[i]),])
    for(j in 1:nrow(X[[i]])){
      temp <- var(as.vector(X[[i]][j,]))
      if(temp != 0){
        U.temp[j,] <- U.temp[j,]/sqrt(temp)
        W[[i]][j,] <- W[[i]][j,]/sqrt(temp)
        muu[j] <- muu[j]/sqrt(temp)
      }
    }
    U[(obs+1):(obs+p[i]),] <- U.temp
    obs <- obs+p[i]
  }
  #Scaling for Y
  temp <- var(as.vector(Y))
  U[nrow(U),] <- U[nrow(U),]/sqrt(temp) * var.y
  thetaS <- 0
  for(i in 1:k){
    theta2[[i]] <- theta2[[i]]/sqrt(temp) * var.y
    thetaS <- thetaS + theta2[[i]] %*% S_i[[i]]
  }
  muu[length(muu)] <- muu[length(muu)]/sqrt(temp) * var.y


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
  muu.new <- NULL
  for(i in 1:k){
    temp <- U[(obs+1):(obs+p[i]),] %*% S_J + W[[i]] %*% S_i[[i]]
    muu.new <- c(muu.new, apply(temp, 1, mean) +  muu[(obs+1):(obs+p[i])])
    X[[i]] <-temp +  muu.new[(obs+1):(obs+p[i])] %*% ones
    obs <- obs+p[i]
  }


  #recreate Y
  Y <- U[nrow(U),] %*% S_J + thetaS + muu[length(muu)] %*% ones

  #### ADDED 1/7/2021: map from natural parameter space to X and Y
  natX <- X
  natY <- Y
  obs <- 0
  error <- list()
  for(i in 1:k){
    error[[i]] <- NA
    if(x.family[i]=="gaussian"){
      fam <- gaussian()
      mu <- fam$linkinv(natX[[i]])
      var_mu <- var(as.vector(mu))
      error[[i]] <- matrix(rnorm(length(mu),0,var_mu*scale.x[i]/(1-scale.x[i])), ncol=ncol(natX[[i]]))
      X[[i]] <- mu + error[[i]]
      muu.new[(obs+1):(obs+p[i])] <- apply(X[[i]], 1, mean)
    }else if(x.family[i]=="binomial"){
      fam <- binomial()
      mu <- fam$linkinv(natX[[i]])
      temp <- rbinom(length(mu), 1, as.vector(mu))
      X[[i]] <- matrix(temp, ncol=ncol(natX[[i]]))
    }else if(x.family[i]=="poisson"){
      fam <- poisson()
      mu <- fam$linkinv(natX[[i]])
      X[[i]] <- matrix(rpois(length(mu), as.vector(mu)), ncol=ncol(natX[[i]]))
    }else{
      stop()
    }
    obs <- obs+p[i]
  }
  if(y.family=="gaussian"){
    fam <- gaussian()
    mu <- fam$linkinv(natY)
    var_mu <- var(as.numeric(mu))
    err <- matrix(rnorm(length(mu), 0,var_mu*scale.y/(1-scale.y)), ncol=ncol(natY))
    Y <- mu + err
    muu.new[length(muu.new)] <- mean(Y)
  }else if(y.family=="binomial"){
    fam <- binomial()
    mu <- fam$linkinv(natY)
    temp <- rbinom(length(mu), 1, as.vector(mu))
    Y <- matrix(temp, ncol=ncol(natY))
  }else if(y.family=="poisson"){
    fam <- poisson()
    mu <- fam$linkinv(natY)
    Y <- matrix(rpois(length(mu), as.vector(mu)), ncol=ncol(natY))
  }else{
    stop()
  }

  if(x.family[1]=="gaussian"){
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
  }

  if(y.family=="gaussian"){
    #Scaling for Y
    temp <- var(as.vector(Y))
    U[nrow(U),] <- U[nrow(U),]/sqrt(temp)
    thetaS <- 0
    for(i in 1:k){
      theta2[[i]] <- theta2[[i]]/sqrt(temp)
      thetaS <- thetaS + theta2[[i]] %*% S_i[[i]]
    }
    err <- err/sqrt(temp)

    #recreate Y
    Y <- U[nrow(U),] %*% S_J + thetaS + err
    Ysig <- U[nrow(U),] %*% S_J + thetaS
    Ymse <- sum((Ysig-Y)^2)/length(Y)
  }


  return(list(X=X, Y=Y, natX=natX, natY=natY,
              S_J=S_J, S_I=S_i, U_I=U[-nrow(U),], W_I=W,
              theta1=U[nrow(U),], theta2=theta2, mu=muu))
}

set.seed(021722)
dat<- sim.data(k=2, p=c(30,30), n=20,
               y.family = "binomial",
               rankJ=1, rankI=c(1,1), prop.causal = c(.5,.5),
               scale.x=c(0.9,0.9), var.y = 8)
SimData.bin <- list()
SimData.bin$X <- dat$X
SimData.bin$Y <- dat$Y

usethis::use_data(SimData.bin, overwrite = TRUE)
