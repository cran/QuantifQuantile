#' @export
#' @title Plot of estimated conditional quantiles using optimal quantization.
#' @name plot.QuantifQuantile
#' @description This function plots the estimated conditional quantiles. 
#' @param x An object of class \code{QuantifQuantile}, which is the result of 
#' \code{\link{QuantifQuantile}} or \code{\link{QuantifQuantile.d2}}.
#' @param col.plot Vector of size \code{length(x$alpha)+1}. The first entry 
#' corresponds to the color of the data points while the other colors are for 
#' the conditional quantiles curves or points.
#' @param \dots Arguments to be passed to \code{\link{par}}.
#' 
#' @details If \code{X} is univariate, the graph is two-dimensional and if 
#' \code{X} is bivariate, it provides a 3D-graph using the \code{\link{rgl}} 
#' package. When only one value for \code{x} is considered, estimated 
#' conditional quantiles are plotted as points. When \code{x} is a grid of 
#' values, they are plotted as curves if \code{d}=1 and surfaces if \code{d}=2.

#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Conditional quantiles estimation through optimal quantization}, 
#' Submitted.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Numerical study of a conditional quantile estimator based on optimal 
#' quantization}, Submitted.

#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and
#'\code{\link{QuantifQuantile.d}}
#' 
#' @examples
#' #for a univariate X
#' set.seed(644936)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
#' plot(res)
#' 
#' set.seed(92536)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,25,by=5),x=1)
#' plot(res)
#' 
#' \dontrun{
#' #for a bivariate X
#' #(a few minuts to execute)
#' set.seed(253664)
#' d <- 2
#' n <- 1000
#' X<-matrix(runif(d*n,-2,2),nr=d)
#' Y<-apply(X^2,2,sum)+rnorm(n)
#' res <- QuantifQuantile.d2(X,Y,testN=seq(80,130,by=10),B=20,tildeB=15)
#' plot(res)
#' 
#' set.seed(193854)
#' d <- 2
#' n <- 1000
#' X<-matrix(runif(d*n,-2,2),nr=d)
#' Y<-apply(X^2,2,sum)+rnorm(n)
#' res <- QuantifQuantile.d2(X,Y,testN=seq(80,150,by=10),x=as.matrix(c(1,0)),
#' B=30,tildeB=20)
#' }
#' 
#' @import rgl
#' @method plot QuantifQuantile
#' @S3method plot QuantifQuantile

plot.QuantifQuantile <- function(x, col.plot = c(1:(length(x$alpha) + 
    1)), ...) {
    stopifnot(class(x)=="QuantifQuantile")
    if (is.vector(x$X)) {
        col.plot[1] <- "grey"
        plot(x$X, x$Y, col = col.plot[1], cex = 0.7, ...)
        if (length(x$x) > 1) {
            for (j in 1:length(x$alpha)) {
                lines(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                  1])
            }
        }
        if (length(x$x) == 1) {
            for (j in 1:length(x$alpha)) {
                points(x$x, x$hatq_opt[j, ], col = col.plot[j + 
                  1], pch = 20, cex = 1)
            }
        }
    } else {
        d <- nrow(x$X)
        if (d == 1) {
            col.plot[1] <- "grey"
            plot(x$X, x$Y, col = col.plot[1], cex = 0.7, ...)
            if (length(x$x) > 1) {
                for (j in 1:length(x$alpha)) {
                  lines(x$x, x$hatq_opt[ j,], col = col.plot[j + 
                    1])
                }
            }
            if (length(x$x) == 1) {
                for (j in 1:length(x$alpha)) {
                  points(x$x, x$hatq_opt[ j, ], col = col.plot[j + 
                    1], pch = 20, cex = 1)
                }
            }
        }
        if(d ==2) {
          plot3d(x$X[1, ], x$X[2, ], x$Y, col = col.plot[1], ...)
          hatq_matrix <- array(0, dim = c(sqrt(dim(x$x)[2]), sqrt(dim(x$x)[2]), 
                                          length(x$alpha)))
          for (i in 1:length(x$alpha)) {
            hatq_matrix[, , i] <- matrix(x$hatq_opt[ i,], ncol = dim(x$x)[2])
          }
          if (length(x$x)/d > 1) {
            for (j in 1:length(x$alpha)) {
              surface3d(unique(x$x[1, ]), unique(x$x[2, ]), 
                        hatq_matrix[, , j], col = col.plot[j + 1])
            }
          }
          if (length(x$x)/d == 1) {
            for (j in 1:length(x$alpha)) {
              points3d(x$x[1, ], x$x[2, ], hatq_matrix[, , 
                                                       j], col = col.plot[j + 1], size = 5)
            }
          }
        }
        if (d > 2) {
            stop("No graphical illustration available when d>2")
        }
    }   
} 
