#' Summary of QuantifQuantile results
#' 
#' This function display a summary of QuantifQuantile results.
#' 
#' This function prints the estimated conditional quantiles q_alpha(x) for each 
#' \code{x} and \code{alpha} considered, as an array.
#'
#' @param object An object of class \code{QuantifQuantile}, which is the result 
#' of the
#' \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} or 
#' \code{\link{QuantifQuantile.d}} functions.
#' @param \dots Not used.
#'
#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Conditional quantiles estimation through optimal quantization}, 
#' Submitted
#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Numerical study of a conditional quantile estimator based on optimal 
#' quantization}, Manuscript in preparation.
#' @seealso \code{\link{QuantifQuantile}}, \code{\link{QuantifQuantile.d2}} and 
#' \code{\link{QuantifQuantile.d}}
#' @seealso \code{\link{plot.QuantifQuantile}}, \code{\link{print.QuantifQuantile}}
#' 
#' @author Isabelle Charlier, Davy Paindaveine, Jerome Saracco
#'
#' @examples
#' set.seed(644936)
#' n <- 300
#' X <- runif(300,-2,2)
#' Y <- X^2+rnorm(n)
#' res <- QuantifQuantile(X,Y,testN=seq(10,30,by=5))
#' summary(res)
#' 
#'
#' @method summary QuantifQuantile
#' @S3method summary QuantifQuantile
#' @export summary.QuantifQuantile
summary.QuantifQuantile <- function(object, ...) {
  stopifnot(class(object)=="QuantifQuantile")
  if(length(object$N_opt)==1){
    cat(paste("** Resulting estimated conditional quantiles with N_opt=",object$N_opt," **"),fill=TRUE)
  }else{
    cat(paste("** Resulting estimated conditional quantiles with N_opt depending on alpha"," **"),fill=TRUE)
    cat(paste("** N_opt="),paste(object$N_opt),paste("**"),fill=TRUE)
  }
  
  if(is.vector(object$X)){
    cat(paste("For each x, the corresponding estimated q_alpha(x) for each alpha\n"))
    result <- array(c(object$x,t(object$hatq_opt)),dim=c(length(object$x),length(object$alpha)+1),dimnames=list(1:length(object$x),c("x",object$alpha)))
    result <- t(result)
  }else{
    cat(paste("For each x (one component by column), the corresponding estimated q_alpha(x) for each alpha\n"))
    result <- array(c(t(object$x),t(object$hatq_opt)),dim=c(ncol(object$x),length(object$alpha)+nrow(object$x)),dimnames=list(1:ncol(object$x),c(rep("x",nrow(object$x)),object$alpha)))  
    result <- t(result)
  }
  print(result)
} 