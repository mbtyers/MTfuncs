#' Only function is to make an input object for caterplot.gls()
#' @description Makes an input object for caterplot.gls().  This part of the algorithm takes some time, so better to do it all once, I say!
#' @param mod A gls model object from gls()
#' @return An object to use in caterplot.gls()
#' @author Matt Tyers
#' @export
caterprep.gls<-function(mod) {
  # returns an object to put into caterplot.gls() for making a caterpillar plot of a gls model
  # input is a gls model object
  
  coef<-summary(mod)$coefficients
  se.coef<-0
  for(i in 1:length(coef)){
    se.coef[i]<-sqrt(summary(mod)$varBeta[i,i])
  }
  out<-list("coef"=coef,"se.coef"=se.coef)
}

#' Makes a caterpillar plot from a gls object
#' 
#' @description Makes a caterpillar plot from a gls object.  Plot shows estimates for each coefficient, plus 95\% CI's for each, defined as +/- 2 se's
#' @param caterprep.out An object returned from caterprep.gls()
#' @param which A vector of which elements to plot.  A value of -1 specifies all but the intercept.
#' @param colors A vector of colors to use
#' @param col.breaks A vector of breakpoints in element numbers for each color, defined as the first element of a given color.  If col=c(1,2,3) and col.breaks=c(1,4,9), then elements 1-3 get color 1, 4-8 get color 2, and 9-end get color 3.
#' @author Matt Tyers
#' @examples
#' # ------- first some data ------ #
#' Factor <- sample(c("A","B","C","D"),1000,replace=T)
#' x <- 1:1000
#' means <- 2*x+1*(Factor=="A")+3*(Factor=="B")+4*(Factor=="C")+2*(Factor=="D")
#' y <- rnorm(1000,means,4)
#' library(nlme)
#' model <- gls(y~Factor*x)
#' 
#' # --------- now the function ---------- #
#' prep.out <- caterprep.gls(model)
#' caterplot.gls(prep.out)
#' @export
caterplot.gls<-function(caterprep.out,which=-1,col.breaks=-1,colors=-1) {
  
  coef <- caterprep.out$coef
  se.coef <- caterprep.out$se.coef
  
  if(which[1]==-1) which<-2:length(coef)
  
  ci.lo <- coef-(2*se.coef)
  ci.hi <- coef+(2*se.coef)
  
  if(col.breaks[1]==-1) col.breaks <- c(0,length(coef))
  cols<-rep(NA,length(coef))
  if(colors[1]==-1) colors<-1:(length(col.breaks)-1)
  for(i in 1:(length(col.breaks)-1)) {
    cols[(col.breaks[i]+1):col.breaks[i+1]]<-colors[i]
  }
  
  plot(c(min(which),max(which)),c(min(ci.lo[which]),max(ci.hi[which]))
       ,col="white",ylab="coefficient",xlab="")
  abline(h=0,lty=3)
  for(i in which) {
    points(i,coef[i],pch='-',col=cols[i])
    lines(c(i,i),c(ci.lo[i],ci.hi[i]),col=cols[i])
  }
}

#' Proportion plots and/or tables
#' @description Creates bar charts with 95\% Jeffries confidence intervals and/or output tables from input vectors of factor data
#' @param x First vector of factor data
#' @param y Second vector of factor data.  If a second vector is used, it makes plots/tables of the proportions of variable 1 for each unique value of variable 2
#' @param numvar Number of input variables (1 or 2)
#' @param xname Name of x to use for plots/tables
#' @param yname Name of y to use for plots/tables
#' @param plot Whether a plot is produced (TRUE or FALSE)
#' @param table Whether a table is produced (TRUE or FALSE)
#' @param kable Whether to use kable format for tables, for use in R markdown (requires knitr package)
#' #' @author Matt Tyers
#' @examples
#' lake <- sample(c("lake 1","lake 2"),20,replace=T)
#' spec <- sample(c("rainbow","cutthroat"),20,replace=T)
#' data.frame(lake,spec)
#' 
#' par(mfrow=c(2,1))
#' props(spec,lake,xname="Species",yname="Lake")
#' @export
props <- function(x,y=rep(1,length(x)),numvar=2,xname=deparse(substitute(x)),yname=deparse(substitute(y)),plot=T,kable=F,table=T) {
  if(kable==T) library(knitr)
  if(numvar==1) {
    cat<-0
    count<-0
    prop<-0
    var.prop<-0
    se.prop<-0
    cv.prop<-0
    ci.jeff.lo<-0
    ci.jeff.hi<-0
    i<-1
    for(a in sort(unique(x))) {
      cat[i] <- a
      count[i] <- length(x[x==a])
      prop[i] <- length(x[x==a])/length(x)
      var.prop[i] <- prop[i]*(1-prop[i])/(length(x)-1)
      se.prop[i] <- sqrt(var.prop[i])
      cv.prop[i] <- se.prop[i]/prop[i]
      ci.jeff.lo[i] <- qbeta(.025,(length(x[x==a])+.5),(length(x)-length(x[x==a])+.5))
      ci.jeff.hi[i] <- qbeta(.975,(length(x[x==a])+.5),(length(x)-length(x[x==a])+.5))
      i<-i+1
    }
    if(table==T){
      cat(length(x),"total","\n")
      out<-data.frame(cat,count,prop,var.prop,se.prop,cv.prop)
      names(out)[1]<-xname
      if(kable==F) print(out)
      if(kable==T) kable(out)
    }
    if(plot==T) {
      se.lo <- prop-se.prop
      se.hi <- prop+se.prop
      plot(c(.5,(length(unique(x))+.5)),c(0,max(ci.jeff.hi)),
           col="white",xaxt='n',xlab=xname,ylab="proportion")
      for(q in 1:length(unique(x))) {
        #points(q,prop[q])
        #lines(c(q,q),c(se.lo[q],se.hi[q]))
        lines(c(q,q),c(ci.jeff.lo[q],ci.jeff.hi[q]))
        lines(c(q-.2,q-.2),c(0,prop[q]))
        lines(c(q-.2,q+.2),c(prop[q],prop[q]))
        lines(c(q+.2,q+.2),c(prop[q],0))
        abline(h=0)
      }
      axis(side=1,at=1:length(unique(x)),labels=sort(unique(x)))
    }
  }
  if(numvar==2) {
    j<-1
    for(b in sort(unique(y))) {
      cat<-0
      count<-0
      prop<-0
      var.prop<-0
      se.prop<-0
      cv.prop<-0
      ci.jeff.lo<-0
      ci.jeff.hi<-0
      i<-1
      for(a in sort(unique(x))) {
        cat[i] <- a
        count[i] <- length(x[x==a&y==b])
        prop[i] <- length(x[x==a&y==b])/length(x[y==b])
        var.prop[i] <- prop[i]*(1-prop[i])/(length(x[y==b])-1)
        se.prop[i] <- sqrt(var.prop[i])
        cv.prop[i] <- se.prop[i]/prop[i]
        ci.jeff.lo[i] <- qbeta(.025,(length(x[x==a&y==b])+.5),(length(x[y==b])-length(x[x==a&y==b])+.5))
        ci.jeff.hi[i] <- qbeta(.975,(length(x[x==a&y==b])+.5),(length(x[y==b])-length(x[x==a&y==b])+.5))
        i<-i+1
      }
      j<-j+1
      if(table==T){
        cat("\n","\n",yname,":",b," - ",length(x[y==b]),"total","\n")
        out<-data.frame(cat,count,prop,var.prop,se.prop,cv.prop)
        names(out)[1]<-xname
        if(kable==F) print(out)
        if(kable==T) kable(out)
      }
      if(plot==T) {
        se.lo <- prop-se.prop
        se.hi <- prop+se.prop
        plot(c(.5,(length(unique(x))+.5)),c(0,max(ci.jeff.hi)),
             col="white",xaxt='n',main=paste(yname,":",b," - ",
                                             length(x[y==b]),"total"),xlab=xname,ylab="proportion")
        for(q in 1:length(unique(x))) {
          #points(q,prop[q])
          #lines(c(q,q),c(se.lo[q],se.hi[q]))
          lines(c(q,q),c(ci.jeff.lo[q],ci.jeff.hi[q]))
          lines(c(q-.2,q-.2),c(0,prop[q]))
          lines(c(q-.2,q+.2),c(prop[q],prop[q]))
          lines(c(q+.2,q+.2),c(prop[q],0))
          abline(h=0)
        }
        axis(side=1,at=1:length(unique(x)),labels=sort(unique(x)))
      }
    }
  }
}

#' Logit and Expit
#' @description Computes the logit log(x/(1/x)) or expit exp(x)/(1+exp(x))
#' @param x Input
#' @author Matt Tyers
#' @examples
#' logit(0.42)
#' expit(0.42)
#' 
#' curve(expit(x),from=-5,to=5)
#' @export
logit <- function(x)  log(x/(1-x))


#' Logit and Expit
#' @description Computes the logit log(x/(1/x)) or expit exp(x)/(1+exp(x))
#' @param x Input
#' @author Matt Tyers
#' @examples
#' logit(0.42)
#' expit(0.42)
#' 
#' curve(expit(x),from=-5,to=5)
#' @export
expit <- function(x)  exp(x)/(1+exp(x))

#' Kernel Density Plot
#' @description Produces a kernel density plot of up to five vectors (this can be biggerized)
#' @param a Input vector
#' @param b Input vector
#' @param c Input vector
#' @param d Input vector
#' @param e Input vector
#' @param xlab Just like usual
#' @param ylab Just like usual
#' @param main Just like usual
#' @param lty Line type.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line type for each line.
#' @param lwd Line width.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line width for each line.
#' @param col Line color.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line color for each line.
#' @param labels A vector of labels for the legend.
#' @param yinfl Inflation factor for determining plotting window size.  Setting \code{yinfl=1} will use the maximum density value.  Setting it slightly higher (the default value is 1.4) allows more space for a legend.
#' @param showN Whether to include sample sizes in the legend.  Defaults to \code{TRUE}.
#' @param ... Additional plotting arguments
#' @author Matt Tyers
#' @examples
#' a <- rnorm(100)
#' b <- rnorm(100)
#' c <- rnorm(100,1,1)
#' d <- rnorm(100,0,2)
#' densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"))
#' densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1)
#' cols <- rainbow(4)
#' densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=cols,lty=1,yinfl=1)
#' 
#' a <- rnorm(10)
#' b <- rnorm(20)
#' c <- rnorm(100,1,1)
#' d <- rnorm(100,0,2)
#' densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1)
#' densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1,bwall=.4)
#' @export
densityplot <- function(a,b,c=NULL,d=NULL,e=NULL,xlab="",ylab="Density",main="",lty=0,lwd=2,col=1,labels=NULL,yinfl=1.4,showN=TRUE,bwall=NULL,...) {
  bdens <- cdens <- ddens <- edens <- NULL
  a <- a[!is.na(a)]
  if(!is.null(bwall)) adens <- density(a,bw=bwall)
  if(is.null(bwall)) adens <- density(a)
  if(!is.null(b)) {
    b <- b[!is.na(b)]
    if(!is.null(bwall)) bdens <- density(b,bw=bwall)
    if(is.null(bwall)) bdens <- density(b)
  }
  if(!is.null(c)) {
    c <- c[!is.na(c)]
    if(!is.null(bwall)) cdens <- density(c,bw=bwall)
    if(is.null(bwall)) cdens <- density(c)
  }
  if(!is.null(d)) {
    d <- d[!is.na(d)]
    if(!is.null(bwall)) ddens <- density(d,bw=bwall)
    if(is.null(bwall)) ddens <- density(d)
  }
  if(!is.null(e)) {
    e <- e[!is.na(e)]
    if(!is.null(bwall)) edens <- density(e,bw=bwall)
    if(is.null(bwall)) endens <- density(e)
  }
  
  num.var <- length(c(a[1],b[1],c[1],d[1],e[1]))
  
  if(lty[1]==0) lty <- 1:num.var
  if(length(lty)==1) lty<-rep(lty,num.var)
  
  if(lwd[1]==0) lwd <- 1:num.var
  if(length(lwd)==1) lwd<-rep(lwd,num.var)
  
  if(col[1]==0) col <- 1:num.var
  if(length(col)==1) col<-rep(col,num.var)
  
  minx <- min(c(a,b,c,d,e))
  maxx <- max(c(a,b,c,d,e))
  maxy <- yinfl*max(c(adens$y,bdens$y,cdens$y,ddens$y,edens$y))
  plot(NA,xlim=c(minx,maxx),ylim=c(0,maxy),main=main,xlab=xlab,ylab=ylab,...=...)
  
  if(is.null(labels)) labels <- 1:numvar
  
  lines(adens,lty=lty[1],lwd=lwd[1],col=col[1])
  if(showN) labels[1] <- paste0(labels[1]," (n=",getn(a),")")
  if(!is.null(b)) {
    lines(bdens,lty=lty[2],lwd=lwd[2],col=col[2])
    if(showN) labels[2] <- paste0(labels[2]," (n=",getn(b),")")
  }
  if(!is.null(c)) {
    lines(cdens,lty=lty[3],lwd=lwd[3],col=col[3])
    if(showN) labels[3] <- paste0(labels[3]," (n=",getn(c),")")
  }
  if(!is.null(d)) {
    lines(ddens,lty=lty[4],lwd=lwd[4],col=col[4])
    if(showN) labels[4] <- paste0(labels[4]," (n=",getn(d),")")
  }
  if(!is.null(e)) {
    lines(edens,lty=lty[5],lwd=lwd[5],col=col[5])
    if(showN) labels[5] <- paste0(labels[5]," (n=",getn(e),")")
  }
  legend(minx,maxy,lwd=lwd,lty=lty,col=col,legend=labels)
}

#' CDF plot
#' @description Produces an empirical CDF plot of up to five vectors (this can be biggerized)
#' @param a Input vector
#' @param b Input vector
#' @param c Input vector
#' @param d Input vector
#' @param e Input vector
#' @param lty Line type.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line type for each line.
#' @param lwd Line width.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line width for each line.
#' @param col Line color.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line color for each line.
#' @param xlab Just like usual
#' @param ylab Just like usual
#' @param main Just like usual
#' @param labels A vector of labels for the legend.
#' @param yinfl Inflation factor for determining plotting window size.  Setting \code{yinfl=1} will use the maximum CDF value (1).  Setting it slightly higher allows more space for a legend.
#' @param showN Whether to include sample sizes in the legend.  Defaults to \code{TRUE}.
#' @param ... Additional plotting arguments
#' @author Matt Tyers
#' @examples
#' a <- rnorm(100)
#' b <- rnorm(100)
#' c <- rnorm(100,1,1)
#' d <- rnorm(100,0,2)
#' cdfplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),lty=1,col=1:4)
#' @export
cdfplot <- function(a,b=NULL,c=NULL,d=NULL,e=NULL,lty=0,lwd=1,col=1,xlab="",ylab="Cumulative Proportion",main="",labels=NULL,yinfl=1,showN=TRUE,...) {
  a <- a[!is.na(a)]
  if(!is.null(b)) b <- b[!is.na(b)]
  if(!is.null(c)) c <- c[!is.na(c)]
  if(!is.null(d)) d <- d[!is.na(d)]
  if(!is.null(e)) e <- e[!is.na(e)]
  
  num.var <- length(c(a[1],b[1],c[1],d[1],e[1]))
  
  if(lty[1]==0) lty <- 1:num.var
  if(length(lty)==1) lty<-rep(lty,num.var)
  
  if(lwd[1]==0) lwd <- 1:num.var
  if(length(lwd)==1) lwd<-rep(lwd,num.var)
  
  if(col[1]==0) col <- 1:num.var
  if(length(col)==1) col<-rep(col,num.var)
  
  plot(c(min(c(a,b,c,d,e)),max(c(a,b,c,d,e))),c(0,yinfl),col="white",xlab=xlab,ylab=ylab,main=main,...=...)
  
  if(is.null(labels)) labels <- 1:num.var
  a.sort <- sort(a)
  a.cum <- seq(0,1,by=(1/(length(a)-1)))
  lines(a.sort,a.cum,lty=lty[1],lwd=lwd[1],col=col[1])
  if(showN) labels[1] <- paste0(labels[1]," (n=",getn(a),")")
  
  if(is.null(b)==F) {
    b.sort <- sort(b)
    b.cum <- seq(0,1,by=(1/(length(b)-1)))
    lines(b.sort,b.cum,lty=lty[2],lwd=lwd[2],col=col[2])
    if(showN) labels[2] <- paste0(labels[2]," (n=",getn(b),")")
  }
  
  if(is.null(c)==F) {
    c.sort <- sort(c)
    c.cum <- seq(0,1,by=(1/(length(c)-1)))
    lines(c.sort,c.cum,lty=lty[3],lwd=lwd[3],col=col[3])
    if(showN) labels[3] <- paste0(labels[3]," (n=",getn(c),")")
  }
  
  if(is.null(d)==F) {
    d.sort <- sort(d)
    d.cum <- seq(0,1,by=(1/(length(d)-1)))
    lines(d.sort,d.cum,lty=lty[4],lwd=lwd[4],col=col[4])
    if(showN) labels[4] <- paste0(labels[4]," (n=",getn(d),")")
  }
  
  if(is.null(e)==F) {
    e.sort <- sort(e)
    e.cum <- seq(0,1,by=(1/(length(e)-1)))
    lines(e.sort,e.cum,lty=lty[5],lwd=lwd[5],col=col[5])
    if(showN) labels[5] <- paste0(labels[5]," (n=",getn(e),")")
  }
  
  abline(h=c(0,1),lty=3)
  legend(min(c(a,b,c,d,e)),yinfl,legend=labels,lty=lty,lwd=lwd,col=col,bg="white")
}

#' Get n
#' @description Returns the length of a vector, excluding \code{NA} values
#' @param x input vector
#' @return Length (numeric)
#' @author Matt Tyers
#' @examples
#' x <- c(1,2,3,NA,5)
#' getn(x)
#' @export
getn <- function(x) sum(!is.na(x))

#' Cumulative frequency plot
#' @description Produces an empirical cumulative frequency (total) plot of up to five vectors (this can be biggerized)
#' @param a Input vector
#' @param b Input vector
#' @param c Input vector
#' @param d Input vector
#' @param e Input vector
#' @param lty Line type.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line type for each line.
#' @param lwd Line width.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line width for each line.
#' @param col Line color.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line color for each line.
#' @param xlab Just like usual
#' @param ylab Just like usual
#' @param main Just like usual
#' @param labels A vector of labels for the legend.
#' @param yinfl Inflation factor for determining plotting window size.  Setting \code{yinfl=1} will use the maximum frequency value.  Setting it slightly higher allows more space for a legend.
#' @param showN Whether to include sample sizes in the legend.  Defaults to \code{TRUE}.
#' @param ... Additional plotting arguments
#' @author Matt Tyers
#' @examples
#' a <- rnorm(100)
#' b <- rnorm(100)
#' c <- rnorm(100,1,1)
#' d <- rnorm(100,0,2)
#' cfplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),lty=1,col=1:4)
#' @export
cfplot <- function(a,b=NULL,c=NULL,d=NULL,e=NULL,lty=0,lwd=1,col=1,xlab="",ylab="Cumulative Frequency",main="",labels=NULL,showN=TRUE,yinfl=1,...) {
  a <- a[!is.na(a)]
  if(!is.null(b)) b <- b[!is.na(b)]
  if(!is.null(c)) c <- c[!is.na(c)]
  if(!is.null(d)) d <- d[!is.na(d)]
  if(!is.null(e)) e <- e[!is.na(e)]
  num.var <- length(c(a[1],b[1],c[1],d[1],e[1]))
  
  if(lty[1]==0) lty <- 1:num.var
  if(length(lty)==1) lty<-rep(lty,num.var)
  
  if(lwd[1]==0) lwd <- 1:num.var
  if(length(lwd)==1) lwd<-rep(lwd,num.var)
  
  if(col[1]==0) col <- 1:num.var
  if(length(col)==1) col<-rep(col,num.var)
  
  a.cum <- NULL
  b.cum <- NULL
  c.cum <- NULL
  d.cum <- NULL
  e.cum <- NULL
  
  plot(c(min(c(a,b,c,d,e)),max(c(a,b,c,d,e))),c(0,yinfl*max(c(length(a),length(b),length(c),length(d),length(e)))),col="white",xlab=xlab,ylab=ylab,main=main,...=...)
  
  a.sort <- sort(a)
  a.cum <- 1:length(a)
  lines(a.sort,a.cum,lty=lty[1],lwd=lwd[1],col=col[1])
  if(showN) labels[1] <- paste0(labels[1]," (n=",getn(a),")")
  
  if(is.null(b)==F) {
    b.sort <- sort(b)
    b.cum <- 1:length(b)
    lines(b.sort,b.cum,lty=lty[2],lwd=lwd[2],col=col[2])
    if(showN) labels[2] <- paste0(labels[2]," (n=",getn(b),")")
  }
  
  if(is.null(c)==F) {
    c.sort <- sort(c)
    c.cum <- 1:length(c)
    lines(c.sort,c.cum,lty=lty[3],lwd=lwd[3],col=col[3])
    if(showN) labels[3] <- paste0(labels[3]," (n=",getn(c),")")
  }
  
  if(is.null(d)==F) {
    d.sort <- sort(d)
    d.cum <- 1:length(d)
    lines(d.sort,d.cum,lty=lty[4],lwd=lwd[4],col=col[4])
    if(showN) labels[4] <- paste0(labels[4]," (n=",getn(d),")")
  }
  
  if(is.null(e)==F) {
    e.sort <- sort(e)
    e.cum <- 1:length(e)
    lines(e.sort,e.cum,lty=lty[5],lwd=lwd[5],col=col[5])
    if(showN) labels[5] <- paste0(labels[5]," (n=",getn(e),")")
  }
  
  if(is.null(labels)) labels <- 1:num.var
  legend(min(c(a,b,c,d,e)),yinfl*max(c(a.cum,b.cum,c.cum,d.cum,e.cum)),legend=labels,lty=lty,lwd=lwd,col=col)
}

#' Factor combiner thingy
#' @description Combines vectors of factors so you don't have to deal with specifying the levels
#' @param a Input vector
#' @param b Input vector
#' @param c Input vector
#' @param d Input vector
#' @param e Input vector
#' @param f Input vector
#' @author Matt Tyers
#' @examples
#' a <- as.factor(c("Apples","Bananas","Apples","Bananas"))
#' b <- as.factor(c("Cucumbers","Bananas","Cucumbers","Durians"))
#' c <- com.fac(a,b)
#' summary(c)
#' @export
com.fac <- function(a,b,c=NULL,d=NULL,e=NULL,f=NULL) {
  out <- as.factor(c(as.character(a),as.character(b),as.character(c),as.character(d),as.character(e),as.character(f)))
  return(out)
}

#' Time format conversion
#' @description pass in a character vector of times formatted hh:mm and it will output a vector of numeric times
#' @param string An input vector of times formatted as text (hh:mm)
#' @param time Value 24 or 12.  Specifies whether input is 24-hour time, or 12-hour + AM/PM
#' @param out Value "dd" or "dh".  Specifies whether output is in decimal days or decimal hours.
#' @author Matt Tyers
#' @examples
#' times <- c("11:00","12:15","22:30","0:15")
#' hhmm2dec(times,time=24,out="dd")
#' hhmm2dec(times,time=24,out="dh")
#' 
#' times12 <- c("11:00AM","12:15PM","2:30PM")
#' hhmm2dec(times12,time=12)
#' @export
hhmm2dec <- function(string,time=24,out="dd") {
  if(is.factor(string)) string <- as.character(string)
  if(time==24) {
    h <- 0
    m <- 0
    split <- strsplit(string,":")
    for(i in 1:length(string)) {
      h[i] <- split[[i]][1]
      m[i] <- split[[i]][2]
    }
    h <- as.numeric(h)
    m <- as.numeric(m)
    if(out=="dd") out.vec <- (h/24) + (m/(24*60))
    if(out=="dh") out.vec <- h + (m/60)
  }
  if(time==12) {
    h <- 0
    pm <- 0
    m1 <- 0
    m10 <- 0
    m01 <- 0
    split <- strsplit(string,":")
    for(i in 1:length(string)) {
      h[i] <- split[[i]][1]
      m1[i] <- split[[i]][2]
      m1split <- strsplit(m1[i],"")
      m10[i] <- m1split[[1]][1]
      m01[i] <- m1split[[1]][2]
      pm[i] <- ifelse((m1split[[1]][3]=="p"|m1split[[1]][3]=="P"),TRUE,FALSE)
    }
    h <- as.numeric(h)
    pm[h>=12] <- F
    m10 <- as.numeric(m10)
    m01 <- as.numeric(m01)
    m <- 10*m10 + m01
    if(out=="dd") out.vec <- (h/24) + (m/(24*60)) + .5*pm
    if(out=="dh") out.vec <- h + (m/60) + 12*pm
  }
  return(out.vec)
}

#' Clears all plots
#' @description Clears all plots to free up memory. 
#' @author David Evans
#' @examples
#' clear.plots()
#' @export
clear.plots <- function() {
  g=dev.list()
  if(is.null(g)==FALSE){dev.off(dev.list()["RStudioGD"])}  
}

#' Color-ramp scatterplot
#' @description Creates a 2-d scatterplot with a third variable plotted on a color ramp
#' @param x X-coordinate vector
#' @param y Y-coordinate vector
#' @param z vector of values to use for the color ramp
#' @param ramp Type of color ramp to use, essentially possible ramps from red, green, and blue.  Possible values are "gr" (default), "rg", "gb", "bg", "br", "rb", "r", "g", "b", "rr", "gg", "bb"
#' @param neut An RGB value (between 0 and 1) for the neutral color(s), those not used on the rgb ramp.  It might take a little twiddling.
#' @param invert A value of TRUE inverts the current color ramp
#' @param legend A value of FALSE suppresses the legend
#' @param cex Universal character expansion factor for plotting (between 0 and 1)
#' @param legend.title Title for the legend
#' @param pch Point character to use (defaults to 16, or solid circles)
#' @param add A value of TRUE adds the points to an active plot, rather than creating a new one
#' @param ... xlab, ylab, and main work like usual
#' @author Matt Tyers
#' @examples
#' x <- rnorm(1000)
#' y <- rnorm(1000)
#' z <- x+y
#' dotsplot.rgb(x,y,z,ramp="rg",legend.title="Legend Title")
#' dotsplot.rgb(x,y,z,ramp="rg",cex=.6,legend.title="Legend Title",invert=T)
#' @export
dotsplot.rgb <- function(x,y,z,ramp="gr",pch=16,add=F,legend=T,cex=1,neut=.5,legend.title=NULL,invert=F,xlab="",ylab="",main="",...) {
  z.s <- (max(z,na.rm=T)-z)/(max(z,na.rm=T)-min(z,na.rm=T))
  legend.z <- seq(max(z,na.rm=T),min(z,na.rm=T),by=((min(z,na.rm=T)-max(z,na.rm=T))/5))
  legend.s <- seq(0,1,by=.2)
  if(invert==T) {
    legend.s <- 1-legend.s
    z.s <- 1-z.s
  }
  if(ramp=="gr") {
    cols <-rgb((1-z.s),z.s,neut)
    cols.legend <- rgb((1-legend.s),legend.s,neut)
  }
  if(ramp=="rg") {
    cols <-rgb(z.s,(1-z.s),neut)
    cols.legend <- rgb(legend.s,(1-legend.s),neut)
  }
  if(ramp=="gb") {
    cols <-rgb(neut,z.s,(1-z.s))
    cols.legend <- rgb(neut,legend.s,(1-legend.s))
  }
  if(ramp=="bg") {
    cols <-rgb(neut,(1-z.s),z.s)
    cols.legend <- rgb(neut,(1-legend.s),legend.s)
  }
  if(ramp=="br") {
    cols <-rgb((1-z.s),neut,z.s)
    cols.legend <- rgb((1-legend.s),neut,legend.s)
  }
  if(ramp=="rb") {
    cols <-rgb(z.s,neut,(1-z.s))
    cols.legend <- rgb(legend.s,neut,(1-legend.s))
  }
  if(ramp=="r") {
    cols <- rgb(z.s,neut,neut)
    cols.legend <- rgb(legend.s,neut,neut)
  }
  if(ramp=="g") {
    cols <- rgb(neut,z.s,neut)
    cols.legend <- rgb(neut,legend.s,neut)
  }
  if(ramp=="b") {
    cols <- rgb(neut,neut,z.s)
    cols.legend <- rgb(neut,neut,legend.s)
  }
  if(ramp=="rr") {
    cols <- rgb(neut,z.s,z.s)
    cols.legend <- rgb(neut,legend.s,legend.s)
  }
  if(ramp=="gg") {
    cols <- rgb(z.s,neut,z.s)
    cols.legend <- rgb(legend.s,neut,legend.s)
  }
  if(ramp=="bb") {
    cols <- rgb(z.s,z.s,neut)
    cols.legend <- rgb(legend.s,legend.s,neut)
  }
  if(add==F) plot(x,y,pch=pch,col=cols,xlab=xlab,ylab=ylab,main=main,cex.axis=cex,cex.main=cex,cex.lab=cex,...=...)
  if(add==T) points(x,y,pch=pch,col=cols)
  if(legend==T) legend(par("usr")[1],par("usr")[4],legend=round(legend.z,2),pch=pch,col=cols.legend,cex=cex,title=legend.title)
}

#' Bootstrap CI
#' @description Creates a bootstrap CI for the mean
#' @param x Vector of data
#' @param B Number of bootstrap replicates to run (defaults to 10000)
#' @param conf Confidence level (defaults to 0.95)
#' @author Matt Tyers
#' @examples
#' x <- rpois(10,10)
#' boot.ci(x)
#' @export
boot.ci <- function(x,B=10000,conf=0.95) {
  boot.mn <- 0
  for(i in 1:B) {
    boot.mn[i] <- mean(sample(x,length(x),replace=T))
  }
  q <- c((1-conf)/2,1-(1-conf)/2)
  ci <- quantile(boot.mn,q)
  return(ci)
}
#cov <- 0
#for(j in 1:10000) {
#  y <- rexp(100,1)
#  ci <- boot.ci(y,B=1000)
#  cov[j] <- 1*(ci[1]<=1 & ci[2]>=1)
#  if(j/100 == floor(j/100)) cat(100*j/10000,"% ... ")
#}
#mean(cov)

#' Creates an assignment statement from a vector
#' @description Creates a c( , , ...) assignment statement from an existing vector.  Useful for copying data directly into BUGS, or hard-coding values into a routine.
#' @param x A vector of data
#' @author Matt Tyers
#' @examples
#' x <- sample(1:10,5)
#' make.c(x)
#' @export
make.c <- function(x) {
  cat("c(")
  for(i in 1:(length(x)-1)) {
    cat(x[i])
    cat(",")
  }
  cat(x[length(x)])
  cat(")")
}

chisq.perm <- function(a,b,reps=10000,plot=F) {
  # permutation test for when sample sizes are too small for a chi-squared contingency table independence/homogeneity test.
  # pass in vectors of both factors.
  
  library(mosaic)
  tally <- tally(~a+b)
  expected <- matrix(NA,nrow=dim(tally)[1],ncol=dim(tally)[2])
  for(i in 1:dim(tally)[1]) {
    for(j in 1:dim(tally)[2]) {
      expected[i,j] <- sum(tally[i,])*sum(tally[,j])/sum(tally)
    }
  }
  xsq <- sum(((tally-expected)^2)/expected)
  xsq.perm <- rep(NA,reps)
  for(k in 1:reps) {
    tally1 <- tally(~sample(a,length(a))+b)
    expected1 <- matrix(NA,nrow=dim(tally1)[1],ncol=dim(tally1)[2])
    for(i in 1:dim(tally1)[1]) {
      for(j in 1:dim(tally1)[2]) {
        expected1[i,j] <- sum(tally1[i,])*sum(tally1[,j])/sum(tally1)
      }
    }
    xsq.perm[k] <- sum(((tally1-expected1)^2)/expected1)
  }
  if(plot==T) {
    hist(xsq.perm)
    abline(v=xsq,lwd=2,col=2)
  }
  pval <- sum(xsq.perm>=xsq)/reps
  return(pval)
}

loop.prog <-function(i,n.rep,n.updates) {
  if(n.updates*i/n.rep == floor(n.updates*i/n.rep)) cat(100*i/reps,"% .. ")
}

#' Numeric integration using Simpson's Rule
#' @description Numeric integration of a univariate function using Simpson's Rule.  NOTE: do NOT use variables with names f, a, b, m, n, h, sum, xx, or simp.
#' @param f An expression to be integrated as a function of t, e.g. "sin(2*pi*t)"
#' @param a Lower limit of integration
#' @param b Upper limit of integration
#' @param m Simpson's parameter (larger number is more accurate but slower)
#' @author Matt Tyers
#' @examples
#' Simp("t^2",0,1)
#' Simp("1.5+3*sin(2*pi*t)",0,0.0833333333,m=2000)
#' @export
Simp <- function(f,a,b,m=100) {n <- 2*m
  h <- (b-a)/n
  t <- a-h
  sum <- 0
  for(i in 0:n) {
    t <- t+h
    xx <- eval(parse(text=f))
    sum <- ifelse(i==0|i==n, sum + xx, ifelse(floor(i/2)==i/2,sum+(2*xx),sum+(4*xx)))
  }
  simp <- sum*h/3
  return(simp)
}

#' Permutation Alternative to 2-sample t-test
#' @description Permutation Alternative to 2-sample 2-test, for inference to a difference in means.
#' @param x1 Data vector
#' @param x2 Data vector
#' @param alternative The alternative hypothesis to test.  Options are \code{"2sided"} (the default),\code{"greater"}, or \code{"less"}.
#' @param reps The number of permutations to run.  Defaults to 10000.
#' @return The permutation p-value (numeric)
#' @author Matt Tyers
#' @examples
#' x1 <- rexp(100,3)
#' x2 <- rexp(100,2)
#' tperm(x1,x2)
#' @export
tperm <- function(x1,x2,alternative="2sided",reps=10000) {
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  datadiff <- mean(x1)-mean(x2)
  pvec <- NA
  x <- c(x1,x2)
  isx1 <- c(rep(T,length(x1)),rep(F,length(x2)))
  for(i in 1:reps) {
    isx1perm <- sample(isx1,length(isx1),replace=F)
    x1perm <- mean(x[isx1perm])
    x2perm <- mean(x[!isx1perm])
    if(alternative=="2sided") pvec[i] <- (abs(x1perm-x2perm)>=abs(datadiff))
    if(alternative=="greater") pvec[i] <- ((x1perm-x2perm)>=(datadiff))
    if(alternative=="less") pvec[i] <- ((x1perm-x2perm)<=(datadiff))
  }
  pval <- mean(pvec)
  return(pval)
}