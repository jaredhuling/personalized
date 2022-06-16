
#' Fit weighted kernel svm model.
#'
#' @description Fits weighted kernel SVM.  To be used for OWL with hinge loss (but can be used more generally)
#'
#' @param y The response vector (either a character vector, factor vector, or numeric vector with values in {-1, 1})
#' @param x The design matrix (not including intercept term)
#' @param weights vector of sample weights for weighted SVM
#' @param C cost of constraints violation, see \code{\link[kernlab]{ksvm}}
#' @param kernel kernel function used for training and prediction. See \code{\link[kernlab]{ksvm}} and \code{\link[kernlab]{kernels}}
#' @param kpar list of hyperparameters for the kernel function. See \code{\link[kernlab]{ksvm}}
#' @param nfolds number of cross validation folds for selecting value of C
#' @param foldid optional vector of values between 1 and nfolds specifying which fold each observation is in. If specified, it will
#' override the \code{nfolds} argument.
#' @param eps penalty nugget parameter. Defaults to \code{1e-8}
#' @param ... extra arguments to be passed to \code{\link[kernlab]{ipop}} from the kernlab package
#' @seealso \code{\link[personalized]{predict.wksvm}} for predicting from fitted \code{weighted.ksvm} objects
#' @importFrom kernlab ipop primal dual kernelMatrix sigest kpar
#' @importFrom kernlab rbfdot polydot tanhdot vanilladot laplacedot besseldot anovadot splinedot
#' @importFrom methods new
#'
#' @export
#' @examples
#'
#' library(kernlab)
#'
#' x <- matrix(rnorm(200 * 2), ncol = 2)
#'
#' y <- 2 * (sin(x[,2]) ^ 2 * exp(-x[,2]) - 0.2 > rnorm(200, sd = 0.1)) - 1
#'
#' weights <- runif(100, max = 1.5, min = 0.5)
#'
#' wk <- weighted.ksvm(x = x[1:100,], y = y[1:100],
#'                     C = c(0.1, 0.5, 1, 2),
#'                     nfolds = 5,
#'                     weights = weights[1:100])
#'
#' pr <- predict(wk, newx = x[101:200,])
#'
#' mean(pr == y[101:200])
#'
weighted.ksvm <- function(y,
                          x,
                          weights,
                          C = c(0.1, 0.5, 1, 2, 10),
                          kernel = "rbfdot",
                          kpar = "automatic",
                          nfolds = 10,
                          foldid = NULL,
                          eps = 1e-8,
                          ...)
{

    nfolds <- as.integer(nfolds[1])
    if (nfolds < 1)
    {
        warning("nfolds must be a positive integer; defaulting to 10")
        nfolds <- 10
    }

    if (is.null(foldid) & nfolds > 1)
    {
        foldid <- sample(rep(seq(nfolds), length = NROW(y)))
    } else
    {
        if (!is.null(foldid))
        {
            foldid <- as.integer(foldid)
            nfolds <- max(foldid)
        }
    }

    if (nfolds >= NROW(y))
    {
        stop("nfolds must be less than the number of observations")
    }

    if (is.factor(y))
    {
        y.vals <- levels(y)
    } else
    {
        y.vals <- sort(unique(y))
    }

    dots <- list(...)

    epsil <- eps

    if (length(y.vals) != 2)
    {
        stop("y must contain only two distinct values.")
    } else
    {
        if (is.character(y.vals))
        {
            y.vec <- ifelse(y == y.vals[2], 1, -1)
        } else if (is.factor(y))
        {
            y.vec <- ifelse(y.vals[y] == y.vals[2], 1, -1)
        } else
        {
            if (any(y.vals != c(-1, 1)))
            {
                stop("y must either be a factor, string vector, or take values {-1, 1}.")
            }
            y.vec <- y
        }
    }

    x    <- as.matrix(x)
    n    <- NROW(y)
    dims <- dim(x)
    n2   <- NROW(weights)

    if (dims[1] != n)
    {
        stop("x must have same number of rows as length of y")
    }
    if (n != n2)
    {
        stop("weights must have same number of rows as length of y")
    }

    possible.kernels <- c("rbfdot",
                          "polydot",
                          "tanhdot",
                          "vanilladot",
                          "laplacedot",
                          "besseldot",
                          "anovadot",
                          "splinedot")
    if (is.character(kernel))
    {
        kernel <- match.arg(kernel, possible.kernels)
        if(is.character(kpar))
        {
            if((kernel == "tanhdot"    ||
                kernel == "vanilladot" ||
                kernel == "polydot"    ||
                kernel == "besseldot"  ||
                kernel == "anovadot"   ||
                kernel == "splinedot")    &&
               kpar == "automatic" )
            {
                kpar <- list()
            }
        }
    }

    ## taken from kernlab
    if (!is.function(kernel))
    {
        if (!is.list(kpar) &&
            is.character(kpar) && (inherits(kernel, "rbfkernel")  ||
                                   inherits(kernel, "laplacedot") ||
                                   kernel        == "laplacedot" ||
                                   kernel        == "rbfdot"))
        {
            kp <- match.arg(kpar,"automatic")
            if(kp == "automatic")
            {
                kpar <- list(sigma = mean(sigest(x,scaled=FALSE)[c(1,3)]))
            }
        }
    }
    if(!is(kernel, "kernel"))
    {
        if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
        kernel <- do.call(kernel, kpar)
    }

    if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")


    # primal: min ||beta||^2 + C\sum_i w_i * max(0, 1 - T_i * f(x_i, beta))

    # T_i in {-1, 1}

    # dual:       max (alpha_i >= 0) sum_i alpha_i - 0.5 * sum_{jk}alpha_j * alpha_k * T_j * T_k K(x_j,x_k)
    # subject to: 0 <= alpha_i <= C * w_i and sum_i alpha_i * T_i = 0

    # which is equivalent to:
    # min (alpha_i >= 0) -sum_i alpha_i + 0.5 * sum_{jk}alpha_j * alpha_k * T_j * T_k K(x_j,x_k)
    # subject to the same constraints

    # which can be solved with ipop() function from kernlab

    K   <- kernelMatrix(kernel, x)
    A   <- matrix(y.vec, ncol = n, nrow = 1)

    best.idx <- 1L
    best.C <- C[1]

    cv.mat <- cv.res <- NULL

    if (nfolds > 1 & length(C) > 1)
    {
        cv.mat <- matrix(0, nrow = length(C), ncol = nfolds)

        for (k in 1:nfolds)
        {
            which.test <- which(foldid == k)

            Hsub <- outer(y.vec[-which.test], y.vec[-which.test], FUN = "*") *
                K[-which.test, -which.test]
            nsub <- nrow(Hsub)

            K.test <- K[-which.test, which.test]

            for (l in 1:length(C))
            {

                cat <- try(fitsub <- ipop_conditioned(c = rep(-1, nsub),
                                                      H = Hsub + epsil * diag(NCOL(Hsub)),
                                                      A = A[,-which.test],
                                                      b = 0,
                                                      r = 0,
                                                      l = rep(0, nsub),
                                                      u = drop(C[l] * weights[-which.test]),
                                                      ...), silent = TRUE)

                if (exists("fitsub"))
                {
                    prim <- kernlab::primal(fitsub)
                    du   <- kernlab::dual(fitsub)

                    pred.test <- sign(unname(colSums(prim * y.vec[-which.test] * K.test) - du))


                    cv.mat[l, k] <- mean(weights[which.test] * (y.vec[which.test] == pred.test))
                } else
                {
                    cv.mat[l, k] <- NA
                }
            }
        }
        cv.res   <- rowMeans(cv.mat, na.rm = TRUE)
        best.idx <- which.max(cv.res)
        best.C   <- C[best.idx]
    }

    fit <- ipop_conditioned(c = rep(-1, n),
                            H = outer(y.vec, y.vec, FUN = "*") * K + epsil * diag(NCOL(K)),
                            A = A,
                            b = 0,
                            r = 0,
                            l = rep(0, n),
                            u = drop(best.C * weights),
                            ...)

    if (fit@how != "converged")
    {
        warning(paste("svm optimization not converged. \nOptimization Status (from kernlab::ipop):\n", fit@how) )
    }

    ret <- list(x          = x,
                y          = y.vec,
                kernel     = kernel,
                primal     = kernlab::primal(fit),
                dual       = kernlab::dual(fit),
                C          = C,
                best.C     = best.C,
                which.best = best.idx,
                cv.mat     = cv.mat,
                cv.res     = cv.res)
    class(ret) <- "wksvm"
    ret
}

#' Prediction function for weighted ksvm objects
#' @description Function to obtain predictions for weighted ksvm objects
#' @rdname predict
#' @seealso \code{\link[personalized]{weighted.ksvm}} for fitting \code{weighted.ksvm} objects
#'
#' @export
predict.wksvm <- function(object, newx, type = c("class", "linear.predictor"), ...)
{
    type <- match.arg(type)
    K    <- kernelMatrix(object$kernel, object$x, newx)

    pred <- unname(colSums(object$primal * object$y * K) - object$dual)

    if (type == "class")
    {
        pred <- sign(pred)
    }
    pred
}

#' Summary of results for weighted ksvm
#'
#' @description Prints summary of results for estimated weighted ksvm
#' @rdname summary
#' @export
summary.wksvm <- function(object, digits = max(getOption('digits')-3, 3), ...)
{
    cat("\n")
    # taken from kernlab
    switch(class(object$kernel),
           "rbfkernel"     = cat(paste("Gaussian Radial Basis kernel function.", "\n","Hyperparameter :" ,"sigma = ",
                                       kpar(object$kernel)$sigma,"\n")),
           "laplacekernel" = cat(paste("Laplace kernel function.", "\n","Hyperparameter :" ,"sigma = ",
                                       kpar(object$kernel)$sigma,"\n")),
           "besselkernel"  = cat(paste("Bessel kernel function.", "\n","Hyperparameter :" ,"sigma = ",
                                       kpar(object$kernel)$sigma,"order = ",
                                       kpar(object$kernel)$order, "degree = ",
                                       kpar(object$kernel)$degree,"\n")),
           "anovakernel"   = cat(paste("Anova RBF kernel function.", "\n","Hyperparameter :" ,"sigma = ",
                                       kpar(object$kernel)$sigma, "degree = ",
                                       kpar(object$kernel)$degree,"\n")),
           "tanhkernel"    = cat(paste("Hyperbolic Tangent kernel function.", "\n","Hyperparameters :","scale = ",
                                       kpar(object$kernel)$scale," offset = ",
                                       kpar(object$kernel)$offset,"\n")),
           "polykernel"    = cat(paste("Polynomial kernel function.", "\n","Hyperparameters :","degree = ",
                                       kpar(object$kernel)$degree," scale = ",
                                       kpar(object$kernel)$scale," offset = ",
                                       kpar(object$kernel)$offset,"\n")),
           "vanillakernel" = cat(paste("Linear (vanilla) kernel function.", "\n")),
           "splinekernel"  = cat(paste("Spline kernel function.", "\n"))
    )
    cat("\n")
    pmat <- rbind(object$C,
                  object$cv.res)
    rownames(pmat) <- c("C", "CV weighted accuracy")
    colnames(pmat) <- c(rep("", length(object$C)))

    print.default(round(pmat, digits), quote = FALSE, right = TRUE, na.print = "NA", ...)
}

## this is taken from kernlab, but fixes
## some of ipop's numerical instabilities
## with a tiny nugget term

##ipop solves the quadratic programming problem
##minimize   c' * primal + 1/2 primal' * H * primal
##subject to b <= A*primal <= b + r
##           l <= x <= u
##           d is the optimizer itself
##returns primal and dual variables (i.e. x and the Lagrange
##multipliers for b <= A * primal <= b + r)
##for additional documentation see
##     R. Vanderbei
##     LOQO: an Interior Point Code for Quadratic Programming, 1992
## Author:      R version Alexandros Karatzoglou, orig. matlab Alex J. Smola
## Created:     12/12/97
## R Version:   12/08/03
## Updated:     13/10/05
## This code is released under the GNU Public License



setGeneric("ipop_conditioned",function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, eps = 1e-6, verb=0) standardGeneric("ipop_conditioned"))
setMethod("ipop_conditioned",signature(H="matrix"),
  function(c, H, A, b, l, u, r, sigf=7, maxiter=40, margin=0.05, bound=10, eps = 1e-6, verb=0)
  {

      if(!is.matrix(H)) stop("H must be a matrix")
      if(!is.matrix(A)&&!is.vector(A)) stop("A must be a matrix or a vector")
      if(!is.matrix(c)&&!is.vector(c)) stop("c must be a matrix or a vector")
      if(!is.matrix(l)&&!is.vector(l)) stop("l must be a matrix or a vector")
      if(!is.matrix(u)&&!is.vector(u)) stop("u must be a matrix or a vector")

      n <- dim(H)[1]

      ## check for a decomposed H matrix
      if(n == dim(H)[2])
          smw <- 0
      if(n > dim(H)[2])
          smw <- 1
      if(n < dim(H)[2])
      {
          smw <- 1
          n <- dim(H)[2]
          H <- t(H)
      }

      if (is.vector(A)) A <- matrix(A,1)
      m <- dim(A)[1]
      primal <- rep(0,n)
      if (missing(b))
          bvec <- rep(0, m)
      ## if(n !=nrow(H))
      ##   stop("H matrix is not symmetric")
      if (n != length(c))
          stop("H and c are incompatible!")
      if (n != ncol(A))
          stop("A and c are incompatible!")
      if (m != length(b))
          stop("A and b are incompatible!")
      if(n !=length(u))
          stop("u is incompatible with H")
      if(n !=length(l))
          stop("l is incompatible with H")

      c <- matrix(c)
      l <- matrix(l)
      u <- matrix(u)

      m <- nrow(A)
      n <- ncol(A)
      H.diag <- diag(H)
      if(smw == 0)
          H.x <- H
      else if (smw == 1)
          H.x <- t(H)
      b.plus.1 <- max(svd(b)$d) + 1
      c.plus.1 <- max(svd(c)$d) + 1
      one.x <- -matrix(1,n,1)
      one.y <- -matrix(1,m,1)
      ## starting point
      if(smw == 0)
          diag(H.x) <- H.diag + 1
      else
          smwn <- dim(H)[2]
      H.y <- diag(1,m)
      c.x <- c
      c.y <- b
      ## solve the system [-H.x A' A H.y] [x, y] = [c.x c.y]
      if(smw == 0)
      {
          AP <- matrix(0,m+n,m+n)
          xp <- 1:(m+n) <= n
          AP[xp,xp] <- -H.x
          AP[xp == FALSE,xp] <- A
          AP[xp,xp == FALSE] <- t(A)
          AP[xp == FALSE, xp== FALSE] <- H.y
          s.tmp <- solve(AP + eps * diag(NCOL(AP)),c(c.x,c.y))
          x <- s.tmp[1:n]
          y <- s.tmp[-(1:n)]
      }
      else
      {
          V <- diag(smwn)
          smwinner <- chol(V + crossprod(H))
          smwa1 <- t(A)
          smwc1 <- c.x
          smwa2 <- smwa1 - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwa1))))
          smwc2 <- smwc1 - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1))))
          y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
          x <- smwa2 %*% y - smwc2
      }

      g <- pmax(abs(x - l), bound)
      z <- pmax(abs(x), bound)
      t <- pmax(abs(u - x), bound)
      s <- pmax(abs(x), bound)
      v <- pmax(abs(y), bound)
      w <- pmax(abs(y), bound)
      p <- pmax(abs(r - w), bound)
      q <- pmax(abs(y), bound)
      mu <- as.vector(crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
      sigfig <- 0
      counter <- 0
      alfa <- 1
      if (verb > 0)	                       # print at least one status report
          cat("Iter    PrimalInf  DualInf  SigFigs  Rescale  PrimalObj  DualObj","\n")

      while (counter < maxiter)
      {
          ## update the iteration counter
          counter <- counter + 1
          ## central path (predictor)
          if(smw == 0)
              H.dot.x <- H %*% x
          else if (smw == 1)
              H.dot.x <- H %*% crossprod(H,x)
          rho <- b - A %*% x + w
          nu <- l - x + g
          tau <- u - x - t
          alpha <- r - w - p
          sigma <- c  - crossprod(A, y) - z + s + H.dot.x
          beta <- y + q - v
          gamma.z <- - z
          gamma.w <- - w
          gamma.s <- - s
          gamma.q <- - q
          ## instrumentation
          x.dot.H.dot.x <-  crossprod(x, H.dot.x)
          primal.infeasibility <- max(svd(rbind(rho, tau, matrix(alpha), nu))$d)/ b.plus.1
          dual.infeasibility <- max(svd(rbind(sigma,t(t(beta))))$d) / c.plus.1
          primal.obj <- crossprod(c,x) + 0.5 * x.dot.H.dot.x
          dual.obj <- crossprod(b,y) - 0.5 * x.dot.H.dot.x + crossprod(l, z) - crossprod(u,s) - crossprod(r,q)
          old.sigfig <- sigfig
          sigfig <- max(-log10(abs(primal.obj - dual.obj)/(abs(primal.obj) + 1)), 0)
          if (sigfig >= sigf) break
          if (verb > 0)		      	# final report
              cat( counter, "\t", signif(primal.infeasibility,6), signif(dual.infeasibility,6), sigfig, alfa, primal.obj, dual.obj,"\n")
          ## some more intermediate variables (the hat section)
          hat.beta <- beta - v * gamma.w / w
          hat.alpha <- alpha - p * gamma.q / q
          hat.nu <- nu + g * gamma.z / z
          hat.tau <- tau - t * gamma.s / s
          ## the diagonal terms
          d <- z / g + s / t
          e <- 1 / (v / w + q / p)
          ## initialization before the big cholesky
          if (smw == 0)
              diag(H.x) <- H.diag + d
          diag(H.y) <- e
          c.x <- sigma - z * hat.nu / g - s * hat.tau / t
          c.y <- rho - e * (hat.beta - q * hat.alpha / p)
          ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
          if(smw == 0){
              AP[xp,xp] <- -H.x
              AP[xp == FALSE, xp== FALSE] <- H.y
              s1.tmp <- solve(AP + eps * diag(NCOL(AP)), c(c.x,c.y))
              delta.x<-s1.tmp[1:n] ; delta.y <- s1.tmp[-(1:n)]
          }
          else
          {
              V <- diag(smwn)
              smwinner <- chol(V + chunkmult(t(H),2000,d))
              smwa1 <- t(A)
              smwa1 <- smwa1 / d
              smwc1 <- c.x / d
              smwa2 <- t(A) - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwa1))))
              smwa2 <- smwa2 / d
              smwc2 <- (c.x - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1)))))/d
              delta.y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
              delta.x <- smwa2 %*% delta.y - smwc2
          }

          ## backsubstitution
          delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
          delta.s <- s * (delta.x - hat.tau) / t
          delta.z <- z * (hat.nu - delta.x) / g
          delta.q <- q * (delta.w - hat.alpha) / p
          delta.v <- v * (gamma.w - delta.w) / w
          delta.p <- p * (gamma.q - delta.q) / q
          delta.g <- g * (gamma.z - delta.z) / z
          delta.t <- t * (gamma.s - delta.s) / s
          ## compute update step now (sebastian's trick)
          alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
          newmu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
          newmu <- mu * ((alfa - 1) / (alfa + 10))^2
          gamma.z <- mu / g - z - delta.z * delta.g / g
          gamma.w <- mu / v - w - delta.w * delta.v / v
          gamma.s <- mu / t - s - delta.s * delta.t / t
          gamma.q <- mu / p - q - delta.q * delta.p / p
          ## some more intermediate variables (the hat section)
          hat.beta <- beta - v * gamma.w / w
          hat.alpha <- alpha - p * gamma.q / q
          hat.nu <- nu + g * gamma.z / z
          hat.tau <- tau - t * gamma.s / s
          ## initialization before the big cholesky
          ##for (  i  in  1 : n H.x(i,i) <- H.diag(i) + d(i) ) {
          ##H.y <- diag(e)
          c.x <- sigma - z * hat.nu / g - s * hat.tau / t
          c.y <- rho - e * (hat.beta - q * hat.alpha / p)

          ## and solve the system [-H.x A' A H.y] [delta.x, delta.y] <- [c.x c.y]
          if (smw == 0)
          {
              AP[xp,xp] <- -H.x
              AP[xp == FALSE, xp== FALSE] <- H.y
              s1.tmp <- solve(AP + eps * diag(NCOL(AP)), c(c.x,c.y))
              delta.x<-s1.tmp[1:n] ; delta.y<-s1.tmp[-(1:n)]
          }
          else if (smw == 1)
          {
              smwc1 <- c.x / d
              smwc2 <- (c.x - (H %*% solve(smwinner,solve(t(smwinner),crossprod(H,smwc1))))) / d
              delta.y <- solve(A %*% smwa2 + H.y , c.y + A %*% smwc2)
              delta.x <- smwa2 %*% delta.y - smwc2
          }
          ## backsubstitution
          delta.w <- - e * (hat.beta - q * hat.alpha / p + delta.y)
          delta.s <- s * (delta.x - hat.tau) / t
          delta.z <- z * (hat.nu - delta.x) / g
          delta.q <- q * (delta.w - hat.alpha) / p
          delta.v <- v * (gamma.w - delta.w) / w
          delta.p <- p * (gamma.q - delta.q) / q
          delta.g <- g * (gamma.z - delta.z) / z
          delta.t <- t * (gamma.s - delta.s) / s
          ## compute the updates
          alfa <- - (1 - margin) / min(c(delta.g / g, delta.w / w, delta.t / t, delta.p / p, delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1))
          x <- x + delta.x * alfa
          g <- g + delta.g * alfa
          w <- w + delta.w * alfa
          t <- t + delta.t * alfa
          p <- p + delta.p * alfa
          y <- y + delta.y * alfa
          z <- z + delta.z * alfa
          v <- v + delta.v * alfa
          s <- s + delta.s * alfa
          q <- q + delta.q * alfa
          ## these two lines put back in ?
          ## mu <- (crossprod(z,g) + crossprod(v,w) + crossprod(s,t) + crossprod(p,q))/(2 * (m + n))
          ## mu <- mu * ((alfa - 1) / (alfa + 10))^2
          mu <- newmu
      }
      if (verb > 0)		      	## final report
          cat( counter, primal.infeasibility, dual.infeasibility, sigfig, alfa, primal.obj, dual.obj)

      ret <- new("ipop_conditioned")               ## repackage the results
      primal(ret) <- x
      dual(ret)   <- drop(y)
      if ((sigfig > sigf) & (counter < maxiter))
          ret@how    <- 'converged'
      else
      {					## must have run out of counts
          if ((primal.infeasibility > 10e5) & (dual.infeasibility > 10e5))
              ret@how    <- 'primal and dual infeasible'
          if (primal.infeasibility > 10e5)
              ret@how    <- 'primal infeasible'
          if (dual.infeasibility > 10e5)
              ret@how    <- 'dual infeasible'
          else					## don't really know
              ret@how    <- 'slow convergence, change bound?'
      }
      ret
  })


setGeneric("chunkmult",function(Z, csize, colscale) standardGeneric("chunkmult"))
setMethod("chunkmult",signature(Z="matrix"),
  function(Z, csize, colscale)
  {
      n <- dim(Z)[1]
      m <- dim(Z)[2]
      d <- sqrt(colscale)
      nchunks <- ceiling(m/csize)
      res <- matrix(0,n,n)

      for( i in 1:nchunks)
      {
          lowerb <- (i - 1) * csize + 1
          upperb <- min(i * csize, m)
          buffer <- t(Z[,lowerb:upperb,drop = FALSE])
          bufferd <- d[lowerb:upperb]
          buffer <- buffer / bufferd
          res <- res + crossprod(buffer)
      }
      return(res)
  })



setClass("ipop_conditioned", representation(primal = "vector",
                                            dual = "numeric",
                                            how = "character"
))

if(!isGeneric("primal")){
    if (is.function("primal"))
        fun <- primal
    else fun <- function(object) standardGeneric("primal")
    setGeneric("primal", fun)
}
setMethod("primal", "ipop_conditioned", function(object) object@primal)
setGeneric("primal<-", function(x, value) standardGeneric("primal<-"))
setReplaceMethod("primal", "ipop_conditioned", function(x, value) {
    x@primal <- value
    x
})

if(!isGeneric("dual")){
    if (is.function("dual"))
        fun <- dual
    else fun <- function(object) standardGeneric("dual")
    setGeneric("dual", fun)
}
setMethod("dual", "ipop_conditioned", function(object) object@dual)
setGeneric("dual<-", function(x, value) standardGeneric("dual<-"))
setReplaceMethod("dual", "ipop_conditioned", function(x, value) {
    x@dual <- value
    x
})


