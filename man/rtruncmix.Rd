% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/MTfuncs_all.R
\name{rtruncmix}
\alias{rtruncmix}
\title{Random Draws from a Truncated Normal Mixture}
\usage{
rtruncmix(length, mu = 0, sig = 1, p = 1, lwr = -1, upr = 1)
}
\arguments{
\item{length}{The number of observations to generate}

\item{mu}{Vector of mean parameters, one for each normal component}

\item{sig}{Vector of standard deviation parameters, one for each normal component}

\item{p}{Vector of proportion parameters, one for each normal component.  These do not have to add to one, and can even represent counts.}

\item{lwr}{Lower truncation bound}

\item{upr}{Upper truncation bound}
}
\value{
A vector of random draws (numeric)
}
\description{
Random Draws from a Truncated Normal Mixture
}
\examples{
draws <- rtruncmix(length=100000, mu=c(0,3,10), sig=c(.5,1,2),p=c(1,1,1), lwr=0, upr=8)
plot(density(draws,bw=.1))
}
\author{
Matt Tyers
}

