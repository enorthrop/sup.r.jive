#Graphical displays of sup.r.jive model results

show.image = function(Image,ylab=''){
  lower = mean(Image)-3*stats::sd(Image)
  upper = mean(Image)+3*stats::sd(Image)
  Image[Image<lower] = lower
  Image[Image>upper] = upper
  withr::local_par(mar=c(1,0,0,0))
  graphics::image(x=1:dim(Image)[2], y=1:dim(Image)[1],
        z=t(Image), zlim = c(lower,upper),
        axes=FALSE,col=gplots::bluered(100),
        xlab="",ylab=ylab)
}


#' Heatmap of supervised JIVE output
#'
#' Display a heatmap of JIVE.pred, sJIVE, or sesJIVE model
#'
#' @param result an object of class sJIVE
#' @param order_by specifies how to order the rows and columns
#'  of the heatmap. #If order_by=-1, orderings are determined
#'  by the outcome. If order_by=0, orderings are determined
#'  by joint structure.  Otherwise, order_by gives the number
#'  of the individual structure dataset to determine the ordering.
#'  In all cases orderings are determined by complete-linkage
#'  hierarchical clustering of Euclidean distances.
#' @param ylab a label for the outcome dataset
#' @param xlab a vector with labels for each X dataset
#' @param ycex a scalar to change the font size of the labels
#'
#'  @details This function takes a fitted sJIVE, sesJIVE, or JIVE.pred
#'  model and plots the results using heatplots. Ensure your plotting window
#'  is sufficiently large prior to running this function to get the
#'  best results visually.
#'
#'
#' @return A heatmap
#' @export
#'
#' @examples
#' data(SimData.norm)
#' fit <- sJIVE(X=SimData.norm$X,Y=SimData.norm$Y,
#'                   rankJ=1,rankA=c(1,1),eta=0.5)
#' plotHeatmap(fit,  ylab="outcome",
#'         xlab=c("Metabolomic", "Proteomic"), ycex=0.9)
#'
plotHeatmap <- function(result, order_by=-1,
                        ylab="Y", xlab=NULL, ycex=1){
  #result is an object of class sJIVE_result
  #order_by specifies how to order the rows and columns
  #  of the heatmap. #If order_by=-1, orderings are determined
  #  by the outcome. If order_by=0, orderings are determined
  #  by joint structure.  Otherwise, order_by gives the number
  #  of the individual structure dataset to determine the ordering.
  #  In all cases orderings are determined by complete-linkage
  #  hiearchichal clustering of Euclidian distances.

  dat <- result$data
  l <- length(dat$X)
  if(is.null(xlab)){xlab=paste0("X",1:l)}
  joint <- indiv <- Yindiv <- list()
  for(i in 1:l){
    joint[[i]] <- result$U_I[[i]] %*% result$S_J
    indiv[[i]] <- result$W_I[[i]] %*% result$S_I[[i]]
    Yindiv[[i]] <- result$theta2[[i]] %*% result$S_I[[i]]
  }
  Yjoint <- result$theta1 %*% result$S_J

  old.par <- graphics::par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(graphics::par(old.par))

  ####Get row/column orderings
  Mat_ColOrder <- do.call(rbind,joint)
  row.orders = list()
  if(order_by==-1){
    for(i in 1:l) row.orders[[i]]=c(dim(dat$X[[i]])[1]:1)
    col.order <- order(dat$Y)  #c(1:dim(dat$X[[i]])[2])
  }
  if(order_by>-1){
    if(order_by>0) {Mat_ColOrder = indiv[[order_by]]}
    col.dist<-stats::dist(t(Mat_ColOrder))
    rm(Mat_ColOrder)
    col.order<-stats::hclust(col.dist)$order
    for(i in 1:l){
      if(order_by==0) {row.dist <- stats::dist(joint[[i]])}
      else{row.dist <- stats::dist(indiv[[i]])}
      row.orders[[i]] <- stats::hclust(row.dist)$order
    }}

  Image_Joint = list()
  for(i in 1:l){ Image_Joint[[i]] = as.matrix(joint[[i]][row.orders[[i]],col.order]) }

  Image_Data = list()
  for(i in 1:l){ Image_Data[[i]] = as.matrix(dat$X[[i]][row.orders[[i]],col.order]) }

  Image_Indiv = list()
  for(i in 1:l){ Image_Indiv[[i]] = as.matrix(indiv[[i]][row.orders[[i]],col.order]) }

  Image_Outcome = as.matrix(dat$Y[col.order])
  Image_YJoint = as.matrix(Yjoint[col.order])
  Image_YIndiv = list()
  for(i in 1:l){ Image_YIndiv[[i]] = as.matrix(Yindiv[[i]][col.order]) }


  graphics::par(mar=c(.1,.1,.1,.1))
  ### Make heatmap of all estimates###
  Heights=c(1,1,1,rep(3,l))
  Widths=c(rep(c(1,4), l+2), .5)
  Vals=c(8,2,5,1,6,3,7,4)
  if(l>2){
    for(i in 3:l){
      n <- length(Vals)
      Vals[which((1:n)%%2==1)] <- Vals[which((1:n)%%2==1)]+1
      Vals[1] <- Vals[1]+1
      Vals <- c(Vals, Vals[n-1]+1, Vals[n]+1)
    }}
  Vals =c(Vals, length(Vals)+1)
  layoutVals =  sapply(Vals, function(x) (x*(3+l)-4):(x*(3+l)))
  graphics::layout(layoutVals, heights = Heights,widths = Widths)
  graphics::layout.show(max(Vals)*(3+l))
  #joint
  CEX <- gplots::textplot("  Joint  ", col="white") ###Plot this first, to get cex size
  graphics::text(.5,1,"Joint", cex=CEX)
  show.image(matrix(Image_YJoint,nrow=1))
  graphics::plot.new()
  graphics::segments(.5,0.1,.5,0.9,lwd=2.5)
  graphics::segments(.5,0.9,.55,0.7,lwd=2.5)
  graphics::segments(.5,0.9,.45,0.7,lwd=2.5)
  for(i in 1:l) show.image(Image_Joint[[i]])
  #data
  graphics::plot.new(); graphics::text(.5,0.3,"Data", cex=CEX)
  show.image(matrix(Image_Outcome,nrow=1),ylab="Y")
  gplots::textplot(" ")
  for(i in 1:l) show.image(Image_Data[[i]],ylab=paste0("X",i))
  #Indiv
  for(i in 1:l){
    graphics::plot.new(); graphics::text(.5,0.3,paste0("Indiv",i), cex=CEX)
    show.image(matrix(Image_YIndiv[[i]],nrow=1))
    graphics::plot.new()
    graphics::segments(.5,0.1,.5,0.9, lwd=2.5)
    graphics::segments(.5,0.9,.55,0.7,lwd=2.5)
    graphics::segments(.5,0.9,.45,0.7,lwd=2.5)
    for(j in 1:l){
      if(i==j){show.image(Image_Indiv[[i]])
      }else{gplots::textplot(" ")}
    }}
  #col3
  gplots::textplot(' ')
  plot(seq(0.1,2.9*pi,0.1),(1-cos(seq(0.1,2.9*pi,0.1)))/25+.6, ylim=c(0,1),
       type="l", lwd=2.5, yaxt="n",
       xaxt="n", bty="n", xlim=c(-pi,4*pi))
  graphics::points(seq(0.1,2.9*pi,0.1), (1-cos(seq(0.1,2.9*pi,0.1)))/25+.45, type="l",
         lwd=2.5)
  gplots::textplot(' ')
  for(i in 1:l){
    plot(seq(0.1,2.9*pi,0.1),(1-cos(seq(0.1,2.9*pi,0.1)))/25+.6, ylim=c(-2,3),
         type="l", lwd=2.5, yaxt="n",
         xaxt="n", bty="n", xlim=c(-pi,4*pi))
    graphics::points(seq(0.1,2.9*pi,0.1), (1-cos(seq(0.1,2.9*pi,0.1)))/25+.4, type="l",
           lwd=2.5)
  }
  #col 5,7,...
  for(i in 1:l){
    gplots::textplot(' ')
    graphics::plot.new()
    graphics::segments(0.2,.5,0.8,0.5,lwd=2.5)
    graphics::segments(.5,.3,.5,.7,lwd=2.5)
    gplots::textplot(' ')
    for(j in 1:l){
      if(i==j){      graphics::plot.new()
        graphics::segments(0.2,.5,0.8,0.5,lwd=2.5)
        graphics::segments(.5,.43,.5,.57,lwd=2.5)
      }else{gplots::textplot(" ")}
    }
  }
  graphics::par(mar=c(0,0,0,0))
  gplots::textplot(' ')
  graphics::plot.new(); graphics::text(1/2,1/2,ylab,srt=90,cex=CEX*0.75*ycex)
  gplots::textplot(' ')
  for(i in 1:l){ graphics::plot.new(); graphics::text(1/2,1/2,xlab[i],srt=90,cex=CEX*0.75)}
}


#' Variance Graph of supervised JIVE output
#'
#' Display a barplot given JIVE.pred, sJIVE, or sesJIVE model
#' to show the amount of variance explained by the model results
#'
#' @param x TBD
#'
#' @return A set of barplots
#' @export
plotVarExplained <- function(x){
  print("Insert Soon x2")
}



#' Plott fitted values from supervised JIVE model
#'
#' Display a fitted versus actual graph given a
#' JIVE.pred, sJIVE, or sesJIVE model
#'
#' @param x TBD
#'
#' @return A diagnostic plot
#' @export
plotFittedValues <- function(x){
  print("Insert Soon x3")
}



#' Summary of supervised JIVE output
#'
#' Display summary data of a JIVE.pred, sJIVE, or sesJIVE model
#'
#' @param object A fitted sJIVE model
#'
#' @details This function gives summary results from any
#' supervised JIVE method. Amount of variance explained
#' is expressed in terms of the standardized Frobenious
#' norm. Partial R-squared values are calculated for the
#' joint and individual components. If rank=1, a z-statistic
#' is calculated to determine the p-value. If rank>1, an F-statistic
#' is calculated.
#'
#' @return Summary measures
#' @export
summary.supJIVE <- function(object){
  #If class sJIVE:
  #cat("\n $ranks \n")
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
      if(i != ii){temp <- temp + a[[i]]}
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

  return(list(eta=object$eta, ranks=tbl_ranks,
              variance=var.table,
              pred.model=tbl))
}
