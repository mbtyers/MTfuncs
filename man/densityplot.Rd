% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/MTfuncs_all.R
\name{densityplot}
\alias{densityplot}
\title{Kernel Density Plot}
\usage{
densityplot(a, b, c = NULL, d = NULL, e = NULL, xlab = "",
  ylab = "Density", main = "", lty = 0, lwd = 2, col = 1,
  labels = NULL, yinfl = 1.4, showN = TRUE, bwall = NULL, ...)
}
\arguments{
\item{a}{Input vector}

\item{b}{Input vector}

\item{c}{Input vector}

\item{d}{Input vector}

\item{e}{Input vector}

\item{xlab}{Just like usual}

\item{ylab}{Just like usual}

\item{main}{Just like usual}

\item{lty}{Line type.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line type for each line.}

\item{lwd}{Line width.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line width for each line.}

\item{col}{Line color.  Accepts one value or a vector of values for each line.  Specifying zero will automatically generate one line color for each line.}

\item{labels}{A vector of labels for the legend.}

\item{yinfl}{Inflation factor for determining plotting window size.  Setting \code{yinfl=1} will use the maximum density value.  Setting it slightly higher (the default value is 1.4) allows more space for a legend.}

\item{showN}{Whether to include sample sizes in the legend.  Defaults to \code{TRUE}.}

\item{...}{Additional plotting arguments}
}
\description{
Produces a kernel density plot of up to five vectors (this can be biggerized)
}
\examples{
a <- rnorm(100)
b <- rnorm(100)
c <- rnorm(100,1,1)
d <- rnorm(100,0,2)
densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"))
densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1)
cols <- rainbow(4)
densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=cols,lty=1,yinfl=1)

a <- rnorm(10)
b <- rnorm(20)
c <- rnorm(100,1,1)
d <- rnorm(100,0,2)
densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1)
densityplot(a,b,c,d,labels=c("apples","bananas","cucumbers","durians"),col=0,lty=1,bwall=.4)
}
\author{
Matt Tyers
}

