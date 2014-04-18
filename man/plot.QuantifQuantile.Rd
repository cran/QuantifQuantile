\name{plot.QuantifQuantile}
\alias{plot.QuantifQuantile}
\title{Plot of estimated conditional quantiles using optimal quantization.}
\usage{
\method{plot}{QuantifQuantile}(x, col.plot = c(1:(length(x$alpha) + 1)), ...)
}
\arguments{
  \item{x}{An object of class \code{QuantifQuantile}, which
  is the result of \code{\link{QuantifQuantile}} or
  \code{\link{QuantifQuantile.d2}}.}

  \item{col.plot}{Vector of size \code{length(x$alpha)+1}.
  The first entry corresponds to the color of the data
  points while the other colors are for the conditional
  quantiles curves or points.}

  \item{\dots}{Arguments to be passed to
  \code{\link{par}}.}
}
\description{
This function plots the estimated conditional quantiles.
}
\details{
If \code{X} is univariate, the graph is two-dimensional and
if \code{X} is bivariate, it provides a 3D-graph using the
\code{\link{rgl}} package. When only one value for \code{x}
is considered, estimated conditional quantiles are plotted
as points. When \code{x} is a grid of values, they are
plotted as curves if \code{d}=1 and surfaces if \code{d}=2.
}
\examples{
#for a univariate X
set.seed(644936)
n <- 300
X <- runif(300,-2,2)
Y <- X^2+rnorm(n)
res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
plot(res)

set.seed(92536)
n <- 300
X <- runif(300,-2,2)
Y <- X^2+rnorm(n)
res <- QuantifQuantile(X,Y,testN=seq(10,25,by=5),x=1)
plot(res)

\dontrun{
#for a bivariate X
#(a few minuts to execute)
set.seed(253664)
d <- 2
n <- 1000
X<-matrix(runif(d*n,-2,2),nr=d)
Y<-apply(X^2,2,sum)+rnorm(n)
res <- QuantifQuantile.d2(X,Y,testN=seq(80,130,by=10),B=20,tildeB=15)
plot(res)

set.seed(193854)
d <- 2
n <- 1000
X<-matrix(runif(d*n,-2,2),nr=d)
Y<-apply(X^2,2,sum)+rnorm(n)
res <- QuantifQuantile.d2(X,Y,testN=seq(80,150,by=10),x=as.matrix(c(1,0)),
B=30,tildeB=20)
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

