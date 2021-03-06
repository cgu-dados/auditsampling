% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nPilot.R
\name{nPilot}
\alias{nPilot}
\title{Calculates the remaing sample size when using pilot sample}
\usage{
nPilot(s2 = NULL, pq = NULL, n1, alpha, moe, N)
}
\arguments{
\item{s2}{The variance, as \eqn{\sigma^2}. Not necessary if \code{pq} is not \code{NULL}.}

\item{pq}{The variance, as \eqn{P(1-P)}. Not necessary if \code{s2} is not \code{NULL}.}

\item{n1}{The size of the pilot sample.}

\item{alpha}{1 - confidence level.}

\item{moe}{Margin of error.}

\item{N}{Population size.}
}
\value{
The final sample is
\deqn{n = n_1 + n_2}
Where \eqn{n_2} is the value returned by this function.
}
\description{
Based on a pilot sample, this function calculates the remain elements to the final sample
given the variance considered.
}
\examples{
nPilot(s2=13.5,n1=30,alpha=0.05,moe=5, N = 306)
nPilot(s2=13.5,n1=30,alpha=0.05,moe=5, N = 5000)
nPilot(pq=0.15,n1=30,alpha=0.1,moe=0.05, N = 306)
nPilot(pq=0.15,n1=30,alpha=0.1,moe=0.05, N = 5000)
}
\references{
COCHRAN, William Gemmell. Sampling techniques-3. 1977.
}
