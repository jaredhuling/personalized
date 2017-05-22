
#' Fitting subgroup identification models
#'
#' @description Fits subgroup identification model class of Chen, et al (2017)
#'
#' @param x The design matrix (not including intercept term)
#' @param y The response vector
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param propensity.func function that inputs the design matrix x and the treatment vector trt and outputs
#' the propensity score, ie Pr(trt = 1 | X = x). Function should take two arguments 1) x and 2) trt. See example below.
#' For a randomized controlled trial this can simply be a function that returns a constant equal to the proportion
#' of patients assigned to the treatment group, i.e.:
#' \code{propensity.func = function(x, trt) 0.5}.
#' @param family family for the response. \code{gaussian} for continuous outcomes, \code{binomial} for binomial outcomes,
#' and \code{cox} for time-to-event outcomes
#' @param loss choice of both the M function from Chen, et al (2017) and potentially the penalty used for variable selection.
#' All \code{loss} options starting with \code{sq_loss} use M(y, v) = (v - y) ^ 2, all options starting with \code{logistic_loss} use
#' the logistic loss: M(y, v) = y * log(1 + exp{-v}), and all options starting with \code{cox_loss} use the negative partial likelihood loss for the Cox PH model.
#' All options ending with \code{lasso} have a lasso penalty added to the loss for variable selection
#' @param method subgroup ID model type. Either the weighting or A-learning method
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group
#' @param larger.outcome.better boolean value of whether a larger outcome is better. Set to \code{TRUE}
#' if a larger outcome is better and set to \code{FALSE} if a smaller outcome is better. Defaults to \code{TRUE}.
#' @param retcall boolean value. if \code{TRUE} then the passed arguments will be saved. Do not set to \code{FALSE}
#' if the \code{validate.subgrp()} function will later be used for your fitted subgroup model. Only set to \code{FALSE}
#' if memory is limited as setting to code{TRUE} saves the design matrix to the fitted object
#' @param ... options to be passed to underlying fitting function. For all \code{loss} options with \code{lasso},
#' this will be passed to \code{cv.glmnet} and for all \code{loss} options with \code{mcp} this will be passed
#' to \code{cv.ncvreg}.
#' @seealso \code{\link[personalized]{validate.subgrp}} for function which creates validation results for subgroup
#' identification models, \code{\link[personalized]{predict.subgroup_fitted}} for a prediction function for fitted models
#' from \code{fit.subgrp}, and \code{\link[personalized]{plot.subgroup_fitted}} for a function which plots
#' results from fitted models
#' from \code{fit.subgrp}.
#' @references Chen, S., Tian, L., Cai, T. and Yu, M. (2017), A general statistical framework for subgroup identification
#' and comparative treatment scoring. Biometrics. doi:10.1111/biom.12676 \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12676/abstract}
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 50
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
#' trt.prob <- exp(xbetat) / (1 + exp(xbetat))
#' trt01    <- rbinom(n.obs, 1, prob = trt.prob)
#'
#' trt      <- 2 * trt01 - 1
#'
#' # simulate response
#' delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
#' xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13]
#' xbeta <- xbeta + delta * trt
#'
#' # continuous outcomes
#' y <- drop(xbeta) + rnorm(n.obs, sd = 2)
#'
#' # binary outcomes
#' y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )
#'
#' # time-to-event outcomes
#' surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
#' cens.time <- exp(rnorm(n.obs, sd = 3))
#' y.time.to.event  <- pmin(surv.time, cens.time)
#' status           <- 1 * (surv.time <= cens.time)
#'
#' # create function for fitting propensity score model
#' prop.func <- function(x, trt)
#' {
#'     # fit propensity score model
#'     propens.model <- cv.glmnet(y = trt,
#'                                x = x, family = "binomial")
#'     pi.x <- predict(propens.model, s = "lambda.min",
#'                     newx = x, type = "response")[,1]
#'     pi.x
#' }
#'
#' subgrp.model <- fit.subgrp(x = x, y = y,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            family = "gaussian",
#'                            loss   = "sq_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model$subgroup.trt.effects
#'
#' subgrp.model.bin <- fit.subgrp(x = x, y = y.binary,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            family = "binomial",
#'                            loss   = "logistic_loss_lasso",
#'                            type.measure = "auc",    # option for cv.glmnet
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model.bin$subgroup.trt.effects
#'
#' library(survival)
#' subgrp.model.cox <- fit.subgrp(x = x, y = Surv(y.time.to.event, status),
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            family = "cox",
#'                            loss   = "cox_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model.cox$subgroup.trt.effects
#'
#'
#' @export
fit.subgrp <- function(x,
                       y,
                       trt,
                       propensity.func,
                       family     = c("gaussian", "binomial", "cox"),
                       loss       = c("sq_loss_lasso",
                                      "logistic_loss_lasso",
                                      "cox_loss_lasso"),
                       method     = c("weighting", "a_learning"),
                       cutpoint   = 0,
                       larger.outcome.better = TRUE,
                       retcall    = TRUE,
                       ...)
{

    family <- match.arg(family)
    loss   <- match.arg(loss)
    method <- match.arg(method)

    y <- drop(y)

    # make sure outcome is consistent with
    # other options selected
    if (class(y) == "Surv")
    {
        if (length(grep("cox", loss))   == 0 |
            length(grep("cox", family)) == 0 )
        {
            stop("loss and family must correspond to a cox model for time-to-event outcomes")
        }
    }

    larger.outcome.better <- as.logical(larger.outcome.better[1])
    retcall <- as.logical(retcall[1])

    # save the passed arguments for later use in validate.subgrp()
    # and plot.subgroup_fitted() functions
    if (retcall)
    {
        this.call <- mget(names(formals()), sys.frame(sys.nframe()))

        this.call$... <- NULL
        this.call <- c(this.call, list(...))
    } else
    {
        this.call <- NULL
    }

    trt         <- as.integer(trt)
    unique.trts <- sort(unique(trt))

    if (length(unique.trts) != 2)    stop("trt must have 2 distinct levels")
    if (any(unique.trts != c(0, 1))) stop("trt should be coded as 0 and 1")

    # check to make sure arguments of propensity.func are correct
    propfunc.names <- sort(names(formals(propensity.func)))
    if (length(propfunc.names) == 2)
    {
        if (any(propfunc.names != c("trt", "x")))
        {
            stop("arguments of propensity.func() should be 'x' and 'trt'")
        }
    } else
    {
        stop("propensity.func() should only have two arguments: 'trt' and 'x'")
    }

    # compute propensity scores
    pi.x <- propensity.func(x = x, trt = trt)

    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity.func should return values between 0 and 1")

    trt2 <- 2 * trt - 1

    # construct modified design matrices
    # and weights
    # depending on what method is used
    if (method == "weighting")
    {
        x.tilde <- trt2 * cbind(1, x)
        wts     <- 1 / (pi.x * (trt == 1) + (1 - pi.x) * (trt == 0))
    } else
    {   # A-learning method
        x.tilde <- (trt - pi.x) * cbind(1, x)
        wts     <- rep(1, nrow(x))
    }

    # identify correct fitting function and call it
    fit_fun      <- paste0("fit_", loss)
    fitted.model <- do.call(fit_fun, list(x = x.tilde, y = y, wts = wts, family = family, ...))

    # save extra results
    fitted.model$call   <- this.call
    fitted.model$family <- family
    fitted.model$loss   <- loss
    fitted.model$method <- method
    fitted.model$benefit.scores <- fitted.model$predict(x)

    # calculate sizes of subgroups and the
    # subgroup treatment effects based on the
    # benefit scores and specified benefit score cutpoint
    if (family == "cox")
    {
        fitted.model$subgroup.trt.effects <- subgrp.benefit(fitted.model$benefit.scores,
                                                            y[,1], trt, cutpoint,
                                                            larger.outcome.better)
    } else
    {
        fitted.model$subgroup.trt.effects <- subgrp.benefit(fitted.model$benefit.scores,
                                                            y, trt, cutpoint,
                                                            larger.outcome.better)
    }

    class(fitted.model) <- "subgroup_fitted"

    fitted.model
}
