#' @export plot.select.N.QuantifQuantile
#' @title Plot of estimated ISE as a function of N
#' @name plot.select.N.QuantifQuantile
#' @description This function illustrates our data driven selection criterion 
#' for \code{N}. It provides the plot of the bootstrap estimated values of
#' integrated squared error ISE(N) versus N.
#' @param x An object of class \code{QuantifQuantile}, which is the result of 
#' \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} or \code{\link{QuantifQuantile.d}}.
#' @param type What type of plot should be drawn.
#' @param \dots Arguments to be passed to \code{\link{par}}.
#' @details This graph allows to adapt the choice of the grid for \code{N}, 
#' called \code{testN}. For example, if the curve is decreasing with \code{N}, it 
#' indicates that the values in \code{testN} are too small and the optimal 
#' \code{N} is larger. 

#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimation through optimal quantization}, 
#' Journal of Statistical Planning and Inference, to appear.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J.,
#' \emph{Conditional quantile estimator based on optimal 
#' quantization: from theory to practice}, Submitted.

#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and 
#' \code{\link{QuantifQuantile.d}}

#' @examples
#' #for a univariate X
#' set.seed(644936)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
#' plot.select.N(res)
#' 
#' \dontrun{
#' #for a bivariate X
#' #(a few minuts to execute)
#' set.seed(253664)
#' n <- 1000
#' X <- matrix(runif(n*2,-2,2),nc=n)
#' Y <- apply(X^2,2,sum)+rnorm(n)
#' res <- QuantifQuantile.d2(X,Y,testN=seq(80,130,by=10),B=20,tildeB=15)
#' plot.select.N(res)
#' }

#' @export plot.select.N
#' @method plot.select.N QuantifQuantile 
#' @aliases plot.select.N plot.select.N.QuantifQuantile
plot.select.N <- function(x, type = "l", ...) UseMethod("plot.select.N")

plot.select.N.QuantifQuantile <- function(x, type = "l", ...) {
    if(length(x$N_opt)==1){
      hatISEmean_N <- apply(x$hatISE_N, 2, mean)
      plot(x$testN, hatISEmean_N, type = type, xlab="N")
      abline(v = x$N_opt, col = 2)
    }else{
      plot(x$testN, x$hatISE_N[1,], type = type, ylim = c(min(x$hatISE_N),
            max(x$hatISE_N)),col=1,lty=1,xlab="N",ylab="hatISE_N")
      j=which(x$N_opt[1]==x$testN)
      lines(c(x$N_opt[1],x$N_opt[1]),c(-1,x$hatISE_N[1,j]),
            col=2)
      if(length(x$alpha)>= 1){
        if(type=="l"){
          for(i in 2:length(x$alpha)){
            lines(x$testN, x$hatISE_N[i,], type = type,col=1,lty=i,lwd=1)
            j=which(x$N_opt[i]==x$testN)
            lines(c(x$N_opt[i],x$N_opt[i]),c(-1,x$hatISE_N[i,j]),
                  col=2)
          }
        }else{
          for(i in 2:length(x$alpha)){
            points(x$testN, x$hatISE_N[i,], type = type, col=1)
            j=which(x$N_opt[i]==x$testN)
            lines(c(x$N_opt[i],x$N_opt[i]),c(-1,x$hatISE_N[i,j]),
                  col=2)
          }
        }
      }
    }
}
 
