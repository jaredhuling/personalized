


#' Creation of augmentation functions
#'
#' @description Creates an augmentation function that optionally utilizes cross-fitting
#'
#' @param family The response type (see options in \code{\link[glmnet]{glmnet}} help file)
#' @param crossfit A logical value indicating whether to use cross-fitting (\code{TRUE}) or not (\code{FALSE}).
#' Cross-fitting is more computationally intensive, but helps to prevent overfitting, see Chernozhukov, et al. (2018)
#' @param nfolds.crossfit An integer specifying the number of folds to use for cross-fitting. Must be greater than 1
#' @param cv.glmnet.args A list of NAMED arguments to pass to the \code{\link[glmnet]{cv.glmnet}} function. For
#' example, \code{cv.glmnet.args = list(type.measure = "mse", nfolds = 10)}. See \code{\link[glmnet]{cv.glmnet}} and \code{\link[glmnet]{glmnet}}
#' for all possible options.
#'
#' @seealso \code{\link[personalized]{fit.subgroup}} for estimating ITRs and \code{\link[personalized]{create.propensity.function}} for creation of propensity functions
#' @return A function which can be passed to the \code{augment.func} argument of the \code{\link[personalized]{fit.subgroup}} function.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018).
#' Double/debiased machine learning for treatment and structural parameters \url{https://arxiv.org/abs/1608.00060}
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 15
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,7] - 0.5 * x[,9]
#' trt.prob <- exp(xbetat) / (1 + exp(xbetat))
#' trt01    <- rbinom(n.obs, 1, prob = trt.prob)
#'
#' trt      <- 2 * trt01 - 1
#'
#' # simulate response
#' # delta below drives treatment effect heterogeneity
#' delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
#' xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
#' xbeta <- xbeta + delta * trt
#'
#' # continuous outcomes
#' y <- drop(xbeta) + rnorm(n.obs, sd = 2)
#'
#' aug.func <- create.augmentation.function(family = "gaussian",
#'                                          crossfit = TRUE,
#'                                          nfolds.crossfit = 10,
#'                                          cv.glmnet.args = list(type.measure = "mae",
#'                                                                nfolds = 5))
#'
#' prop.func <- create.propensity.function(crossfit = TRUE,
#'                                         nfolds.crossfit = 10,
#'                                         cv.glmnet.args = list(type.measure = "auc",
#'                                                               nfolds = 5))
#'
#' subgrp.model <- fit.subgroup(x = x, y = y,
#'                              trt = trt01,
#'                              propensity.func = prop.func,
#'                              augment.func = aug.func,
#'                              loss   = "sq_loss_lasso",
#'                              nfolds = 10)    # option for cv.glmnet (for ITR estimation)
#'
#' summary(subgrp.model)
#'
#' @importFrom stats model.matrix
#' @export
create.augmentation.function <- function(family, crossfit = TRUE, nfolds.crossfit = 10, cv.glmnet.args = NULL)
{
    if (family == "binomial")
    {
        tm <- "auc"
    } else
    {
        tm <- "mse"
    }

    nfolds.crossfit <- as.integer(nfolds.crossfit[1])
    stopifnot(nfolds.crossfit > 1)

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    cv.glmnet.args

    augment.func <- function(x, y, trt)
    {
        glmnet_aug_kfold_crossfit(x = x, y = y, trt = trt, use.crossfitting = crossfit,
                                  K = nfolds.crossfit, cv.glmnet.args = cv.glmnet.args, family = family,
                                  predtype = "link", interactions = TRUE)
    }

    augment.func
}


#' Creation of propensity fitting function
#'
#' @description Creates an propensity function that optionally utilizes cross-fitting
#'
#' @param crossfit A logical value indicating whether to use cross-fitting (\code{TRUE}) or not (\code{FALSE}).
#' Cross-fitting is more computationally intensive, but helps to prevent overfitting, see Chernozhukov, et al. (2018)
#' @param nfolds.crossfit An integer specifying the number of folds to use for cross-fitting. Must be greater than 1
#' @param cv.glmnet.args A list of NAMED arguments to pass to the \code{\link[glmnet]{cv.glmnet}} function. For
#' example, \code{cv.glmnet.args = list(type.measure = "mse", nfolds = 10)}. See \code{\link[glmnet]{cv.glmnet}} and \code{\link[glmnet]{glmnet}}
#' for all possible options.
#'
#' @seealso \code{\link[personalized]{fit.subgroup}} for estimating ITRs and \code{\link[personalized]{create.propensity.function}} for creation of propensity functions
#' @return A function which can be passed to the \code{augment.func} argument of the \code{\link[personalized]{fit.subgroup}} function.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018).
#' Double/debiased machine learning for treatment and structural parameters \url{https://arxiv.org/abs/1608.00060}
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 15
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,7] - 0.5 * x[,9]
#' trt.prob <- exp(xbetat) / (1 + exp(xbetat))
#' trt01    <- rbinom(n.obs, 1, prob = trt.prob)
#'
#' trt      <- 2 * trt01 - 1
#'
#' # simulate response
#' # delta below drives treatment effect heterogeneity
#' delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
#' xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
#' xbeta <- xbeta + delta * trt
#'
#' # continuous outcomes
#' y <- drop(xbeta) + rnorm(n.obs, sd = 2)
#'
#' aug.func <- create.augmentation.function(family = "gaussian",
#'                                          crossfit = TRUE,
#'                                          nfolds.crossfit = 10,
#'                                          cv.glmnet.args = list(type.measure = "mae",
#'                                                                nfolds = 5))
#'
#' prop.func <- create.propensity.function(crossfit = TRUE,
#'                                         nfolds.crossfit = 10,
#'                                         cv.glmnet.args = list(type.measure = "mae",
#'                                                               nfolds = 5))
#'
#' subgrp.model <- fit.subgroup(x = x, y = y,
#'                              trt = trt01,
#'                              propensity.func = prop.func,
#'                              augment.func = aug.func,
#'                              loss   = "sq_loss_lasso",
#'                              nfolds = 10)    # option for cv.glmnet (for ITR estimation)
#'
#' summary(subgrp.model)
#'
#' @export
create.propensity.function <- function(crossfit = TRUE, nfolds.crossfit = 10, cv.glmnet.args = NULL)
{

    tm <- "auc"

    nfolds.crossfit <- as.integer(nfolds.crossfit[1])
    stopifnot(nfolds.crossfit > 1)

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    cv.glmnet.args

    propensity.func <- function(x, trt)
    {
        glmnet_propensity_kfold_crossfit(x = x, trt = trt, use.crossfitting = crossfit,
                                         K = nfolds.crossfit, cv.glmnet.args = cv.glmnet.args)
    }

    propensity.func
}





glmnet_aug_kfold_crossfit <- function(x, y, trt, wts = NULL,
                                      use.crossfitting = TRUE,
                                      K = 10,
                                      predtype = c("link", "response"),
                                      family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                                      interactions = TRUE, cv.glmnet.args = NULL)
{

    predtype <- match.arg(predtype)
    family   <- match.arg(family)


    if (family == "binomial")
    {
        tm <- "auc"
    } else
    {
        tm <- "mse"
    }

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    if (is.null(wts))
    {
        wts <- rep(1, NROW(x))
    }

    if (interactions)
    {
        ## full model for nonzeroness
        df_all <- data.frame(x, trt = trt)
        df_1   <- data.frame(x, trt = 1)
        df_0   <- data.frame(x, trt = -1)

        mm_all <- model.matrix(~x*trt-1, data = df_all)
        mm_1   <- model.matrix(~x*trt-1, data = df_1)
        mm_0   <- model.matrix(~x*trt-1, data = df_0)
    } else
    {
        mm_all <- x
    }

    n <- NROW(mm_all)

    predvec <- numeric(n)

    if (use.crossfitting)
    {
        foldid = sample(rep(seq(K), length = n))

        for (i in seq(K))
        {
            which <- foldid == i

            ## full model for nonzeroness
            # glmfit_zero_main  <- cv.glmnet(y            = y[!which],
            #                                x            = mm_all[!which,,drop=FALSE],
            #                                weights      = wts[!which],
            #                                family       = family,
            #                                parallel     = FALSE,
            #                                type.measure = type.measure)

            glmfit_zero_main <- do.call(cv.glmnet, c(list(y = y[!which], x = mm_all[!which,,drop=FALSE],
                                                          weights = wts[!which], family = family), cv.glmnet.args))

            if (interactions)
            {
                ## get predictions for trt = 1 & -1
                pred1_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_1[which,,drop=FALSE], s = "lambda.min", type = predtype)))
                pred0_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_0[which,,drop=FALSE], s = "lambda.min", type = predtype)))

                predvec[which] <- 0.5 * (pred1_zerr + pred0_zerr)
            } else
            {
                ## get predictions for trt = 1 & -1
                pred_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_all[which,,drop=FALSE], s = "lambda.min", type = predtype)))

                predvec[which] <- pred_zerr
            }

        }
    } else
    {
        glmfit_zero_main <- do.call(cv.glmnet, c(list(y = y, x = mm_all,
                                                      weights = wts, family = family), cv.glmnet.args))

        if (interactions)
        {
            ## get predictions for trt = 1 & -1
            pred1_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_1, s = "lambda.min", type = predtype)))
            pred0_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_0, s = "lambda.min", type = predtype)))

            predvec <- 0.5 * (pred1_zerr + pred0_zerr)
        } else
        {
            ## get predictions for trt = 1 & -1
            pred_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_all, s = "lambda.min", type = predtype)))

            predvec <- pred_zerr
        }
    }

    predvec
}



glmnet_propensity_kfold_crossfit <- function(x, trt, use.crossfitting = TRUE, K = 10, cv.glmnet.args = NULL)
{

    n <- NROW(x)

    tm <- "auc"

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    propensvec <- numeric(n)

    if (use.crossfitting)
    {
        foldid = sample(rep(seq(K), length = n))

        for (i in seq(K))
        {
            which <- foldid == i

            ## propensity score model fit on K-1 folds
            # glmfit_propens  <- cv.glmnet(y = trt[!which], x = x[!which,,drop=FALSE],
            #                              family = "binomial", #parallel = TRUE,
            #                              type.measure = type.measure)

            glmfit_propens <- do.call(cv.glmnet, c(list(y = trt[!which], x = x[!which,,drop=FALSE],
                                                        family = "binomial"), cv.glmnet.args))

            ## get propensity scores for the held out fold
            propensvec[which] <- unname(drop(predict(glmfit_propens, newx = x[which,,drop=FALSE],
                                                     s = "lambda.1se", type = "response")))
        }

    } else
    {
        glmfit_propens <- do.call(cv.glmnet, c(list(y = trt, x = x, family = "binomial"), cv.glmnet.args))

        ## get propensity scores for the held out fold
        propensvec <- unname(drop(predict(glmfit_propens, newx = x, s = "lambda.1se", type = "response")))
    }


    ## propensity scores will never be outside of 0 or 1 and
    ## shouldn't have missing values, but this code is a safety
    ## check just in case
    propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
    propensvec[propensvec <= 0] <- 1e-5
    propensvec[propensvec >= 1] <- 1 - 1e-5

    propensvec
}
