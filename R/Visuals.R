#Graphical displays of sup.r.jive model results


###### Generic Functions  ######

#' Plot Heatmap for a supervised JIVE model
#'
#' A heat map to visually display the joint and invidual components
#' in a fitted model.
#'
#' @param result An object of class "sJIVE", "JIVEpred", "sesJIVE".
#' @param ... further arguments passed to or from other methods.
#'
#' @details This function takes a fitted sJIVE, sesJIVE, or JIVE.pred
#' model and plots the results using heatplots. Ensure your plotting window
#' is sufficiently large prior to running this function to get the
#' best results visually.
#'
#' @return A heat map.
#' @export
plotHeatmap <- function(result, ...) {
  UseMethod("plotHeatmap")
}


#' Plot fitted values from supervised JIVE model
#'
#' Display adiagnostic plots given a
#' JIVE.pred, sJIVE, or sesJIVE model
#'
#' @param result A fitted JIVEpred model
#' @param ... further arguments passed to or from other methods
#'
#' @details Depending if the outcome is Gaussian or binary, different diagnostic plots
#' will be generated. If the outcome is Gaussian, a residual plot and a Q-Q plot will
#' be created. If the outcome is binary, two plots will be created for each joint and
#' individual component: the first will plot a Loess curve, and the second
#' will contain a density plot. Both plots help show the separation, or lack thereof,
#' between the component and the outcome.
#'
#' @return Diagnostic plot(s)
#' @export
#'
#' @examples
#' \dontrun{
#' #Let fit be a fitted sJIVE, JIVE.pred, or sesJIVE model
#' plotFittedValues(fit)
#' }
plotFittedValues <- function(result, ...){
  UseMethod("plotFittedValues")
}

#' Variance Graph of supervised JIVE output
#'
#' Display a barplot given JIVE.pred, sJIVE, or sesJIVE model
#' to show the amount of variance explained by the model results
#'
#' @param result An object of class "sJIVE", "JIVEpred", "sesJIVE".
#' @param col a vector containing the 3 colors that should be used for the barplot.
#' @param ... further arguments passed to or from other methods.
#'
#' @details A barplot for each \code{X} dataset and the outcome will be graphed. If \code{y}
#' is not continuous (e.g. if it is binary), a barplot for \code{y} will not be plotted.
#'
#' @return A set of barplots
#' @export
plotVarExplained <- function(result, col, ...){
  UseMethod("plotVarExplained")
}



#### Plot Heatmap for each class #####

#' plotHeatmap.sJIVE
#'
#' @describeIn plotHeatmap
#'
#' @details description of plotHeatmap.sJIVE here.
#'
#'
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
#' @param ... further arguments passed to or from other methods
#'
#' @export
#' @examples
#' data(SimData.norm)
#' fit <- sJIVE(X=SimData.norm$X,Y=SimData.norm$Y,
#'                   rankJ=1,rankA=c(1,1),eta=0.5)
#' plotHeatmap(fit,  ylab="outcome",
#'         xlab=c("Metabolomic", "Proteomic"), ycex=0.9)
#'
plotHeatmap.sJIVE <- function(result, order_by=-1,
                        ylab="Y", xlab=NULL, ycex=1, ...){
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
  #graphics::layout.show(max(Vals)*(3+l))
  #joint
  CEX <-1 #gplots::textplot("  Joint  ", col="white") ###Plot this first, to get cex size
  graphics::plot.new(); graphics::text(.5,.3,"Joint", cex=CEX)
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



#' plotHeatmap.JIVEpred
#'
#' @describeIn plotHeatmap
#'
#' @details description of plotHeatmap.JIVEpred here.
#'
#' @export
#' @examples
#' data(SimData.norm)
#' fit <- JIVE.pred(X=SimData.norm$X,Y=SimData.norm$Y,
#'                   rankJ=1,rankA=c(1,1))
#' plotHeatmap(fit,  ylab="outcome",
#'         xlab=c("Metabolomic", "Proteomic"), ycex=0.9)
#'
plotHeatmap.JIVEpred <- function(result, order_by=-1,
                              ylab="Y", xlab=NULL, ycex=1, ...){
  #result is an object of class sJIVE_result
  #order_by specifies how to order the rows and columns
  #  of the heatmap. #If order_by=-1, orderings are determined
  #  by the outcome. If order_by=0, orderings are determined
  #  by joint structure.  Otherwise, order_by gives the number
  #  of the individual structure dataset to determine the ordering.
  #  In all cases orderings are determined by complete-linkage
  #  hiearchichal clustering of Euclidian distances.

  dat <- list(X=result$jive.fit$data,
              Y=result$data.matrix$Y)
  l <- length(dat$X)
  if(is.null(xlab)){xlab=paste0("X",1:l)}
  joint <- indiv <- Yindiv <- list()
  for(i in 1:l){
    joint[[i]] <- result$jive.fit$joint[[i]]
    indiv[[i]] <- result$jive.fit$individual[[i]]

    obs <- which(grepl(paste0("I",i),names(result$mod.fit$coefficients)))
    Yindiv[[i]] <- as.matrix(result$data.matrix[,obs]) %*% result$mod.fit$coefficients[obs]
  }
  r_j <- result$jive.fit$rankJ
  Yjoint <- as.matrix(result$data.matrix[,c(2:(r_j+1))], ncol=r_j) %*% result$mod.fit$coefficients[c(2:(r_j+1))]

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
  #graphics::layout.show(max(Vals)*(3+l))
  #joint
  CEX <-1 #gplots::textplot("  Joint  ", col="white") ###Plot this first, to get cex size
  graphics::plot.new(); graphics::text(.5,.3,"Joint", cex=CEX)
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


####### Plot Variance Graph for each class #######

#' plotVarExplained.JIVEpred
#'
#' @describeIn plotVarExplained
#'
#' @details description of plotVarExplained.JIVEpred here.
#'
#' @export
#' @examples
#' \dontrun{
#' #Let fit be a fitted sJIVE, JIVE.pred, or sesJIVE model
#' plotVarExplained(fit, col=c("grey20", "grey43", "grey65"))
#' }
plotVarExplained.JIVEpred <- function(result, col=c("grey20", "grey43", "grey65"), ...){

  old.par <- graphics::par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(graphics::par(old.par))

  s.result <- summary(result)

  if(result$family == "gaussian"){
  l <- ncol(s.result$variance)-2

  graphics::par(mar=c(5,4,4,0))
  graphics::layout(matrix(c(1,2),1,2),heights=c(5,5),widths=c(5,2))
  graphics::barplot(as.matrix(s.result$variance[,-1]),col = col,main = "Variation Explained",
          names.arg=names(s.result$variance)[-1])
  graphics::par(mar=c(0,0,0,0))
  graphics::plot.new()
  graphics::legend(x=0.05,y=0.8,legend=c('Joint','Individual','Residual'),
         bty = "n",fill= col)
  }else{
    warning("Variation in Y cannot be displayed since Y is not Gaussian")
    l <- ncol(s.result$variance)-2

    graphics::par(mar=c(5,4,4,0))
    graphics::layout(matrix(c(1,2),1,2),heights=c(5,5),widths=c(5,2))
    graphics::barplot(as.matrix(s.result$variance[,-1]),col = col,main = "Variation Explained",
                      names.arg=names(s.result$variance)[-1])
    graphics::par(mar=c(0,0,0,0))
    graphics::plot.new()
    graphics::legend(x=0.05,y=0.8,legend=c('Joint','Individual','Residual'),
                     bty = "n",fill= col)
  }
}

#' plotVarExplained.sJIVE
#'
#' @describeIn plotVarExplained
#'
#' @details description of plotVarExplained.sJIVE here.
#'
#' @export
plotVarExplained.sJIVE <- function(result, col=c("grey20", "grey43", "grey65"), ...){

  old.par <- graphics::par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(graphics::par(old.par))

  s.result <- summary(result)
    l <- ncol(s.result$variance)-2

    graphics::par(mar=c(5,4,4,0))
    graphics::layout(matrix(c(1,2),1,2),heights=c(5,5),widths=c(5,2))
    graphics::barplot(as.matrix(s.result$variance[,-1]),col = col,main = "Variation Explained",
                      names.arg=names(s.result$variance)[-1])
    graphics::par(mar=c(0,0,0,0))
    graphics::plot.new()
    graphics::legend(x=0.05,y=0.8,legend=c('Joint','Individual','Residual'),
                     bty = "n",fill= col)
}


###### Plot Fitted Values for each class ######

#' plotFittedValues.JIVEpred
#'
#' @describeIn plotFittedValues
#'
#' @details description of plotFittedValues.JIVEpred here.
#'
#' @param graph A value: 0, 1, or 2
#' @export
plotFittedValues.sJIVE <- function(result, graph=0, ...){
  old.par <- graphics::par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(graphics::par(old.par))

  rsd <- result$fittedY - result$data$Y
  ylims <- max(abs(rsd))

  if(graph==0){
    graphics::par(mfrow=c(1,2), mar=c(4.5,4,3,1))
  }
  if(graph !=2){
     plot(result$fittedY, rsd, ylab="Residuals",
         xlab="Fitted Y values", main="Residuals vs. Fitted",
         ylim=c(-ylims, ylims))
    graphics::abline(h = 0, lty = 3, col = "gray")
  }
  if(graph != 1){
    stats::qqnorm(rsd, ylab="", ylim=c(-ylims, ylims))
    stats::qqline(rsd, lty = 3, col = "gray50")
  }
}



#' plotFittedValues.sJIVE
#'
#' @describeIn plotFittedValues
#'
#' @details description of plotFittedValues.sJIVE here.
#'
#' @export
plotFittedValues.JIVEpred <- function(result, graph=0, ...){
  old.par <- graphics::par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(graphics::par(old.par))

  if(result$family == "gaussian"){
  rsd <- result$mod.fit$fitted.values - result$data.matrix$Y
  ylims <- max(abs(rsd))

  if(graph==0){
    graphics::par(mfrow=c(1,2), mar=c(4.5,4,3,1))
  }
  if(graph !=2){
    plot(result$mod.fit$fitted.values, rsd, ylab="Residuals",
         xlab="Fitted Y values", main="Residuals vs. Fitted",
         ylim=c(-ylims, ylims))
    graphics::abline(h = 0, lty = 3, col = "gray")
  }
  if(graph != 1){
    stats::qqnorm(rsd, ylab="", ylim=c(-ylims, ylims), xlab="Theoretical Quantiles")
    stats::qqline(rsd, lty = 3, col = "gray50")
  }
  }else if(result$family == "binomial"){
      l <- ncol(result$data.matrix)-1
      dat <- result$data.matrix

      if(l %% 2 == 0){
      for(i in 1:(l/2)){
        x <- names(dat)[2*i-1+1]
        p1 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=Y)) +
                ggplot2::geom_point(alpha=0.33) +
                ggplot2::geom_smooth(method = "loess") + ggplot2::xlab(x)
        p3 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), group=Y, color=Y, fill=Y)) +
                ggplot2::geom_density(alpha=0.4) + ggplot2::xlab(x)

        x <- names(dat)[2*i+1]
        p2 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=Y)) +
          ggplot2::geom_point(alpha=0.33) +
          ggplot2::geom_smooth(method = "loess") + ggplot2::xlab(x)
        p4 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), group=Y, color=Y, fill=Y)) +
          ggplot2::geom_density(alpha=0.4) + ggplot2::xlab(x)
        gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)
      }
      }else{
        for(i in 1:((l-1)/2)){
          x <- names(dat)[2*i-1+1]
          p1 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=Y)) +
            ggplot2::geom_point(alpha=0.33) +
            ggplot2::geom_smooth(method = "loess") + ggplot2::xlab(x)
          p3 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), group=Y, color=Y, fill=Y)) +
            ggplot2::geom_density(alpha=0.4) + ggplot2::xlab(x)

          x <- names(dat)[2*i+1]
          p2 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=Y)) +
            ggplot2::geom_point(alpha=0.33) +
            ggplot2::geom_smooth(method = "loess") + ggplot2::xlab(x)
          p4 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), group=Y, color=Y, fill=Y)) +
            ggplot2::geom_density(alpha=0.4) + ggplot2::xlab(x)
          gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)
        }
        x <- names(dat)[l+1]
        p1 <- ggplot2::ggplot(dat, ggplot2::aes(x=get(x), y=Y)) +
          ggplot2::geom_point(alpha=0.33) +
          ggplot2::geom_smooth(method = "loess") + ggplot2::xlab(x)
        p3 <- ggplot2::ggplot(dat,ggplot2::aes(x=get(x), group=Y, color=Y, fill=Y)) +
          ggplot2::geom_density(alpha=0.4) + ggplot2::xlab(x)
        gridExtra::grid.arrange(p1, p3, ncol=1)
      }
  }
}

