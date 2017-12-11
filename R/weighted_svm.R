
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
#' @param ... extra arguments to be passed to \code{\link[kernlab]{ipop}} from the kernlab package
#' @seealso \code{\link[personalized]{predict.wksvm}} for predicting from fitted \code{weighted.ksvm} objects
#' @importFrom kernlab ipop primal dual kernelMatrix sigest kpar
#' @importFrom kernlab rbfdot polydot tanhdot vanilladot laplacedot besseldot anovadot splinedot
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
#' wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(0.1, 0.5, 1, 2, 10),
#'                     weights = weights[1:100])
#'
#' pr <- predict(wk, newx = x[101:200,])
#'
#' mean(pr == y[101:200])
#'
weighted.ksvm <- function(y,
                          x,
                          weights,
                          C = c(0.1, 0.5, 1, 5, 10),
                          kernel = "rbfdot",
                          kpar = "automatic",
                          nfolds = 10,
                          foldid = NULL,
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
            is.character(kpar) && (class(kernel) == "rbfkernel"  ||
                                   class(kernel) == "laplacedot" ||
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

                fitsub <- kernlab::ipop(c = rep(-1, nsub),
                                        H = Hsub,
                                        A = A[,-which.test],
                                        b = 0,
                                        r = 0,
                                        l = rep(0, nsub),
                                        u = drop(C[l] * weights[-which.test]),
                                        ...)
                prim <- kernlab::primal(fitsub)
                du   <- kernlab::dual(fitsub)

                pred.test <- sign(unname(colSums(prim * y.vec[-which.test] * K.test) - du))


                cv.mat[l, k] <- mean(weights[which.test] * (y.vec[which.test] == pred.test))
            }
        }
        cv.res   <- rowMeans(cv.mat)
        best.idx <- which.max(cv.res)
        best.C   <- C[best.idx]
    }

    fit <- kernlab::ipop(c = rep(-1, n),
                         H = outer(y.vec, y.vec, FUN = "*") * K,
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
