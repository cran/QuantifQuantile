\name{plot.select.N.QuantifQuantile}
\alias{plot.select.N}
\alias{plot.select.N.QuantifQuantile}
\title{Plot of estimated MSE as a funcion of N}
\usage{
\method{plot.select.N}{QuantifQuantile}(x, type = "l", ...)
}
\arguments{
  \item{x}{An object of class \code{QuantifQuantile}, which
  is the result of \code{\link{QuantifQuantile}}.}

  \item{type}{What type of plot should be drawn.}

  \item{\dots}{Arguments to be passed to
  \code{\link{par}}.}
}
\description{
This function illustrates our data driven selection
criterion for \code{N}. It provides the plot of the
bootstrap estimated values of MSE(N) versus N.
}
\details{
This graph allows to adapt the choice of the grid for
\code{N}, called \code{testN}. For example, if the curve is
decreasing with \code{N}, it indicates that the values in
\code{testN} are too small and the optimal \code{N} is
larger.
}
\examples{
#for a univariate X
set.seed(644936)
n <- 300
X <- runif(300,-2,2)
Y <- X^2+rnorm(n)
res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
plot.select.N(res)

\dontrun{
#for a bivariate X
#(a few minuts to execute)
set.seed(253664)
n <- 1000
X <- matrix(runif(n*2,-2,2),nc=n)
Y <- apply(X^2,2,sum)+rnorm(n)
res <- QuantifQuantile.d2(X,Y,testN=seq(80,130,by=10),B=20,tildeB=15)
plot.select.N(res)
}
}
\references{
Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
\emph{Conditional quantiles estimation through optimal
quantization}, Manuscript in preparation

Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
\emph{Numerical study of a conditional quantile estimator
based on optimal quantization}, Manuscript in preparation
}
\seealso{
\code{\link{QuantifQuantile}},
\code{\link{QuantifQuantile.d2}} and
\code{\link{QuantifQuantile.d}}
}

