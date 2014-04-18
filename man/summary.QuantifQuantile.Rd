\name{summary.QuantifQuantile}
\alias{summary.QuantifQuantile}
\title{Summary of QuantifQuantile results}
\usage{
\method{summary}{QuantifQuantile}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{QuantifQuantile},
  which is the result of the \code{\link{QuantifQuantile}},
  \code{\link{QuantifQuantile.d2}} or
  \code{\link{QuantifQuantile.d}} functions.}

  \item{\dots}{Not used.}
}
\description{
This function display a summary of QuantifQuantile results.
}
\details{
This function prints the estimated conditional quantiles
q_alpha(x) for each \code{x} and \code{alpha} considered,
as an array.
}
\examples{
set.seed(644936)
n <- 300
X <- runif(300,-2,2)
Y <- X^2+rnorm(n)
res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
summary(res)
}
\author{
Isabelle Charlier, Davy Paindaveine, Jerome Saracco
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

\code{\link{plot.QuantifQuantile}},
\code{\link{print.QuantifQuantile}}
}

