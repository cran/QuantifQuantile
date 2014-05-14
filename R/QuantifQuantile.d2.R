#' @export
#' @title QuantifQuantile for X bivariate
#' @name QuantifQuantile.d2
#' @description Estimation of conditional quantiles using optimal quantization 
#' when \code{X} is bivariate.
#' 
#' @details \itemize{\item This function calculates estimated conditional 
#' quantiles with a method based on optimal quantization when the covariate is 
#' bivariate. The matrix of covariate \code{X} must have two rows (dimension). 
#' For other dimensions, see \code{\link{QuantifQuantile}} or 
#' \code{\link{QuantifQuantile.d}}. The argument \code{x} must also have two rows.
#' \item The criterion for selecting the number of quantizers is implemented in 
#' this function. The user has to choose a grid \code{testN} of possible values 
#' in which \code{N} will be selected. It actually minimizes some bootstrap 
#' estimated version of the MSE (Mean Squared Error). More precisely, for 
#' \code{N} fixed, it calculates the sum according to \code{alpha} of 
#' \code{hatMSE_N} and then minimizes the resulting vector to get \code{N_opt}.
#'  However, the user can choose to select a different value of \code{N_opt} for
#'  each \code{alpha} by setting \code{same_N=FALSE}. In this case, the vector 
#'  \code{N_opt} is obtained by minimizing each column of \code{hatME_N} 
#'  separately. The reason why \code{same_N=TRUE} by default is that taking 
#'  \code{N_opt} according to \code{alpha} could provide crossing condtional 
#'  quantile curves (rarely observed for not too close values of \code{alpha}. 
#'  The function \code{\link{plot.select.N.QuantifQuantile}} 
#'  illustrates the selection of \code{N_opt}. If the graph is not globally convex, the arguments 
#'  \code{testN} should be adapted.}
#'  
#' @param X matrix of covariates.
#' @param Y vector of response variables.
#' @param alpha vector of order of the quantiles.
#' @param x matrix of values for \code{x} in q_alpha(x).
#' @param testN grid of values of \code{N} that will be tested.
#' @param p L_p norm optimal quantization.
#' @param B number of bootstrap replications for the bootstrap estimator.
#' @param tildeB number of bootstrap replications for the choice of \code{N}.
#' @param same_N whether to use the same value of \code{N} for each \code{alpha}
#' (\code{TRUE} by default).
#'  
#' @return An object of class \code{QuantifQuantile} which is a list with the 
#' following components:
#' @return \item{hatq_opt}{A matrix containing the estimated conditional 
#' quantiles. The number of columns is the number of considered values for \code{x} 
#' and the number of rows the size of the order vector \code{alpha}.}
#' @return \item{N_opt}{Optimal selected value for \code{N}. An integer if 
#' \code{same_N}=TRUE and a vector of integers of length \code{length(alpha)} 
#' otherwise.}
#' @return \item{hatMSE_N}{The matrix of estimated MSE provided by our selection
#'  criterion for \code{N} before taking the mean according to \code{alpha}. The
#'   number of columns is then \code{length(testN)} and the number of rows 
#'   \code{length(alpha)}.}
#' @return \item{hatq_N}{A 3-dimensional array containing the estimated 
#' conditional quantiles for each considered value for \code{N}.}
#' @return \item{X}{The matrix of covariates.}
#' @return \item{Y}{The vector of response variables.}
#' @return \item{x}{The considered vector of values for \code{x} in q_alpha(x).}
#' @return \item{alpha}{The considered vector of order for the quantiles.}
#' @return \item{testN}{The considered grid of values for \code{N} that were tested.}

#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Conditional quantiles estimation through optimal quantization}, 
#' Submitted.
#' @references Charlier, I. and Paindaveine, D. and Saracco, J. (2014),
#' \emph{Numerical study of a conditional quantile estimator based on optimal 
#' quantization}, Manuscript in preparation.
#' 
#' @seealso \code{\link{QuantifQuantile}} and \code{\link{QuantifQuantile.d}} 
#' for other dimensions.
#' @seealso \code{\link{plot.select.N.QuantifQuantile}} for the \code{N} 
#' selection criterion.
#' @seealso \code{\link{plot.QuantifQuantile}}, 
#' \code{\link{print.QuantifQuantile}}, \code{\link{summary.QuantifQuantile}}
#' 
#' @examples
#' \dontrun{
#' #(a few minuts to execute)
#' set.seed(253664)
#' n <- 1000
#' X <- matrix(runif(n*2,-2,2),ncol=n)
#' Y <- apply(X^2,2,sum)+rnorm(n)
#' res3 <- QuantifQuantile.d2(X,Y,testN=seq(90,140,by=10),B=20,tildeB=15)
#' res4 <- QuantifQuantile.d2(X,Y,testN=seq(90,150,by=10),B=20,tildeB=15,same_N=FALSE)
#' }

QuantifQuantile.d2 <- function(X, Y, alpha = c(0.05, 0.25, 0.5, 
    0.75, 0.95), x = matrix(c(rep(seq(min(X[1, ]), max(X[1, ]), 
    length = 20), 20), sort(rep(seq(min(X[2, ]), max(X[2, ]), 
    length = 20), 20))), nrow = 2, byrow = TRUE), testN = c(110, 
    120, 130, 140, 150), p = 2, B = 50, tildeB = 20, same_N=TRUE) {
    if (!is.numeric(X)) 
        stop("X must be numeric")
    if (!is.numeric(Y)) 
        stop("Y must be numeric")
    if (!is.numeric(x)) 
        stop("x must be numeric")
    if (!is.matrix(X)) 
        stop("X must be a matrix with d rows")
    if (!is.vector(Y)) 
        stop("Y must be a vector")
    if (!is.matrix(x)) 
        stop("X must be a matrix with d rows")
    if (!all(floor(testN) == testN & testN > 0)) 
        stop("testN must have entire positive entries")
    if (all(alpha > 0 & alpha < 1) == FALSE) 
        stop("alpha must be strictly between 0 and 1")
    if ((!(floor(B) == B)) | (B <= 0)) 
        stop("B must be a positive entire")
    if ((!(floor(tildeB) == tildeB)) | (tildeB <= 0)) 
        stop("tildeB must be a positive entire")
    if (p < 1) 
        stop("p must be at least 1")
    if (!is.logical(same_N))
      stop("same_N must be logical")
    n <- ncol(X)
    d <- nrow(X)
    
    if (nrow(x) != d) 
        stop("X must be a matrix with d rows")
    
    m <- length(x)/d  #number of vectors x for which we estimate q_alpha(x)
    
    hatMSE_N <- array(0, dim = c(length(alpha), length(testN)))
    hatq_N <- array(0, dim = c(m, length(alpha), length(testN)))

    if (nrow(X) == n) 
        {
            X = t(X)
        }  #X doit avoir d lignes
    if (nrow(x) == m) 
        {
            x = t(x)
        }  #X doit avoir d lignes
    
    primeX <- array(X[, sample(c(1:n), n * (B + tildeB), replace = T)], 
        dim = c(d, n, (B + tildeB)))
    
    for (jj in 1:length(testN)) {
        N <- testN[jj]
        
        hatX <- choice.grid(X, N, B = B, tildeB = tildeB)$opti_grid
        
        # projection of the sample X on the B+tildeB optimal grids
        
        projXboot <- array(0, dim = c(d, n, B + tildeB))
        # index of the grid on which X is projected
        iminx <- array(0, dim = c(n, B + tildeB))
        for (i in 1:n) {
            RepX <- array(rep(X[, i], N * (B + tildeB)), dim = c(d, 
                N, B + tildeB))
            Ax <- sqrt(apply((RepX - hatX)^2, c(2, 3), sum))
            iminx[i, ] <- apply(Ax, 2, which.min)
            mx <- array(0, dim = c(d, B + tildeB, 3))
            for (k in 1:d) {
                mx[k, , ] <- matrix(c(rep(k, B + tildeB), iminx[i, 
                  ], c(1:(B + tildeB))), ncol = 3)
                projXboot[k, i, ] <- hatX[mx[k, , ]]
            }
        }
        # estimation of q_alpha(x) for N fixed 
        # save the B+tildeB estimation of q_alpha(x)
        Hatq <- array(0, dim = c(m, length(alpha), B + tildeB))
        # save by Voronoi cell
        Hatq_cell <- array(0, dim = c(N, length(alpha), (B + 
            tildeB)))
        
        proj_gridx_boot = function(z) {
            proj <- array(0, dim = c(d, B + tildeB))
            Repz <- array(rep(z, N * (B + tildeB)), dim = c(d, 
                N, B + tildeB))
            A <- sqrt(apply((Repz - hatX)^2, c(2, 3), sum))
            i <- apply(A, 2, which.min)
            mx <- array(0, dim = c(d, B + tildeB, 3))
            for (k in 1:d) {
                mx[k, , ] <- matrix(c(rep(k, B + tildeB), i, 
                  c(1:(B + tildeB))), ncol = 3)
                proj[k, ] <- hatX[mx[k, , ]]
            }
            proj
        }
        
        projection_x <- apply(x, FUN = proj_gridx_boot, MARGIN = 2)
        projection_x <- array(projection_x, dim = c(d, B + tildeB, 
            m))
        
        repY <- matrix(rep(Y, B + tildeB), nrow = n)
        
        # calculation of the conditional quantile for each cell. 
        # Since any point of a cell is projected on the center of this cell, 
        # the corresponding conditional quantiles are equal
        
        for (i in 1:N) {
            for (j in 1:(B + tildeB)) {
                init = proc.time()
                a <- which(projXboot[1, , j] == hatX[1, i, j] & 
                  projXboot[2, , j] == hatX[2, i, j])
                proc.time() - init
                if (length(a) > 0) {
                  Hatq_cell[i, , j] <- quantile(repY[a, j], probs = alpha)
                }
            }
        }
        
        # we now identify the cell in which belongs each x to associate the
        # corresponding value of conditional quantiles
        
        identification <- function(z) {
            identification <- array(0, dim = c(B + tildeB, 1))
            i <- which(z[1] == x[1, ] & z[2] == x[2, ])
            for (j in 1:(B + tildeB)) {
                identification[j, ] <- which(projection_x[1, 
                  j, i] == hatX[1, , j] & projection_x[2, j, 
                  i] == hatX[2, , j])
            }
            identification
        }
        
        identification_projection_x <- apply(x, FUN = identification, 
            2)
        
        
        for (i in 1:length(alpha)) {
            for (j in 1:m) {
                r <- matrix(c(identification_projection_x[, j], 
                  rep(i, B + tildeB), c(1:(B + tildeB))), ncol = 3)
                Hatq[j, i, ] <- Hatq_cell[r]
            }
        }
        
        # the final estimation is the mean of the B estimations
        hatq <- array(0, dim = c(m, length(alpha)))
        hatq <- apply(Hatq[, , c(1:B), drop = FALSE], c(1, 2), 
            mean)
        
        i <- which(N == testN)
        hatq_N[, , i] <- hatq
        
        # the last tilde B are used to estimate the MSE
        HATq <- array(rep(hatq, tildeB), dim = c(m, length(alpha), 
            tildeB))
        hatMSE <- (HATq - Hatq[, , c((1 + B):(B + tildeB)), drop = FALSE])^2
        hatMSE <- apply(hatMSE, 2, sum)/(m * tildeB)
        hatMSE_N[, i] <- hatMSE
        
        print(N)
    }
    
    if(same_N){
      #choice of optimal N
      hatMSEmean_N <- apply(hatMSE_N, 2, mean)
      i_opt <- which.min(hatMSEmean_N)
      #optimal value for N chosen as minimizing the sum of hatMSE for the 
      # different alpha's
      N_opt <- testN[i_opt]
      
      # table of the associated estimated conditional quantiles
      hatq_opt <- hatq_N[, , i_opt, drop = F]
      hatq_opt <- matrix(hatq_opt, ncol = length(alpha))
      hatq_opt <- t(hatq_opt)
    }else{
      #choice of optimal N
      i_opt <- apply(hatMSE_N, 1, which.min)
      #optimal value for N chosen as minimizing the sum of hatMSE for the 
      # different alpha's
      N_opt <- testN[i_opt]
      # table of the associated estimated conditional quantiles
      hatq_opt <- array(0, dim = c(length(alpha), dim(x)[2]))
      for(i in 1:length(alpha)){
        hatq_opt[i, ] <- hatq_N[,i , i_opt[i]]
      }
    }
    
    output <- list(hatq_opt = hatq_opt, N_opt = N_opt, 
        hatMSE_N = hatMSE_N, hatq_N = hatq_N, X = X, Y = Y, x = x, 
        alpha = alpha, testN = testN)
    class(output) <- "QuantifQuantile"
    output
    ############################ 
    
}  # fin de la fonction 
