\name{QuantifQuantile.d}
\alias{QuantifQuantile.d}
\title{QuantifQuantile for general X}
\usage{
QuantifQuantile.d(X, Y, x, alpha = c(0.05, 0.25, 0.5, 0.75, 0.95),
  testN = c(35, 40, 45, 50, 55), p = 2, B = 50, tildeB = 20)
}
\arguments{
  \item{X}{matrix of covariables.}

  \item{Y}{vector of response variables.}

  \item{x}{matrix of values for \code{x} in q_alpha(x).}

  \item{alpha}{vector of order of the quantiles.}

  \item{testN}{grid of values of \code{N} that will be
  tested.}

  \item{p}{L_p norm optimal quantization.}

  \item{B}{number of bootstrap replications for the
  bootstrap estimator.}

  \item{tildeB}{number of bootstrap replications for the
  choice of \code{N}.}
}
\value{
An object of class \code{QuantifQuantile} which is a list
with the following components:

\item{hatq_opt}{A matrix containing the estimated
conditional quantiles. The number of rows is the number of
considered values for \code{x} and the number of columns
the size of the order vector \code{alpha}.}

\item{N_opt}{The optimal selected value for \code{N}.}

\item{hatMSEmean_N}{The vector of estimated MSE provided by
our selection criterion for \code{N} (after taking the mean
according to \code{alpha}). Its size is
\code{length(testN)}. The optimal value \code{N_opt}
minimizes this vector.}

\item{hatMSE_N}{The matrix of estimated MSE provided by our
selection criterion for \code{N} before taking the mean
according to \code{alpha}. The number of columns is then
\code{length(testN)} and the number of rows
\code{length(alpha)}.}

\item{hatq_N}{A 3-dimensional array containing the
estimated conditional quantiles for each considered value
for \code{N}.}

\item{X}{The matrix of covariables.}

\item{Y}{The vector of response variables.}

\item{x}{The considered vector of values for \code{x} in
q_alpha(x).}

\item{alpha}{The considered vector of order for the
quantiles.}

\item{testN}{The considered grid of values for \code{N}
that were tested.}
}
\description{
Estimation of conditional quantiles using optimal
quantization when \code{X} is d-dimensional.
}
\details{
\itemize{\item This function calculates estimated
conditional quantiles with a method based on optimal
quantization for any dimension for the covariable. The
matrix of covariable \code{X} must have \code{d} rows
(dimension). For particular cases of \code{d} =1 or 2, it
is strongly recommended to use
\code{\link{QuantifQuantile}} and
\code{\link{QuantifQuantile.d2}} respectively
(computationally faster). The argument \code{x} must also
have \code{d} rows. \item The criterion for selecting the
number of quantizers is implemented in this function. The
user has to choose a grid of possible values \code{testN}
in which \code{N} will be selected. It actually minimizes
some bootstrap estimate version of the MSE (Mean Squared
Error). The function
\code{\link{plot.select.N.QuantifQuantile}} illustrates it.
If the graph is not globally convex, the arguments
\code{testN} should be adapted.}
}
\examples{
\dontrun{
set.seed(491625)
n <- 500
X <- runif(n,-2,2)
Y <- X^2+rnorm(n)
x <- seq(min(X),max(X),length=100)
res <- QuantifQuantile.d(X,Y,x,testN=seq(10,30,by=5))
}
\dontrun{
set.seed(272422)
n <- 1000
X <- matrix(runif(n*2,-2,2),ncol=n)
Y <- apply(X^2,2,sum)+rnorm(n)
x1 <- seq(min(X[1,]),max(X[1,]),length=20)
x2 <- seq(min(X[2,]),max(X[2,]),length=20)
x <- matrix(c(rep(x1,20),sort(rep(x2,20))),nrow=nrow(X),byrow=TRUE)
res <- QuantifQuantile.d(X,Y,x,testN=seq(90,140,by=10),B=20,tildeB=15)
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
\code{\link{QuantifQuantile}} and
\code{\link{QuantifQuantile.d2}} for particular dimensions
one and two.

\code{\link{plot.select.N.QuantifQuantile}} for the
\code{N} selection criterion.

\code{\link{plot.QuantifQuantile}},
\code{\link{print.QuantifQuantile}},
\code{\link{summary.QuantifQuantile}}
}

