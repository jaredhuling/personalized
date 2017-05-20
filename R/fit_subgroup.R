
#' Fitting subgroup identification models
#'
#' @description Fits subgroup identification model class of Chen, et al (2017)
#'
#' @param x The design matrix (not including intercept term)
#' @param y The response vector
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param pi.x vector of treatment status probabilities (ie Pr(trt = 1 | X = x)). If data is from a
#' randomized controlled trial, these probabilities are constants. In an observational study \code{pi.x}
#' is the propensity score.
#' @param family family for the response. \code{gaussian} for continuous outcomes, \code{binomial} for binomial outcomes,
#' and \code{cox} for time-to-event outcomes
#' @param loss choice of both the M function from Chen, et al (2017) and potentially the penalty used for variable selection.
#' All \code{loss} options starting with \code{sq_loss} use M(y, v) = (v - y) ^ 2, all options starting with \code{logistic_loss} use
#' the logistic loss: M(y, v) = y * log(1 + exp{-v}), and all options starting with \code{cox_loss} use the negative partial likelihood loss for the Cox PH model.
#' All options ending with \code{lasso} have a lasso penalty added to the loss for variable selection
#' @param method subgroup ID model type. Either the weighting or A-learning method
#' @param ... options to be passed to underlying fitting function. For all \code{loss} options with \code{lasso},
#' this will be passed to \code{cv.glmnet} and for all \code{loss} options with \code{mcp} this will be passed
#' to \code{cv.ncvreg}.
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 100
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
#' delta <- 2 * (0.25 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12])
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
#' surv.time <- exp(-xbeta + rnorm(n.obs, sd = 1))
#' cens.time <- exp(rnorm(n.obs, sd = 3))
#' y.time.to.event  <- pmin(surv.time, cens.time)
#' status           <- 1 * (surv.time <= cens.time)
#'
#' # fit propensity score model
#' propens.model <- cv.glmnet(y = (1+trt)/2,
#'                            x = x, family = "binomial")
#' pi.x <- predict(propens.model, s = "lambda.min",
#'                 newx = x, type = "response")[,1]
#'
#' subgrp.model <- fit.subgrp(x = x, y = y, trt = trt01, pi.x = pi.x,
#'                            family = "gaussian",
#'                            loss   = "sq_loss_lasso",
#'                            gamma  = 2,              # option for cv.ncvreg
#'                            nfolds = 5)              # option for cv.ncvreg
#'
#' @export
fit.subgrp <- function(x,
                       y,
                       trt,
                       pi.x,
                       family     = c("gaussian", "binomial", "cox"),
                       loss       = c("sq_loss_lasso",
                                      "logistic_loss_lasso",
                                      "cox_loss_lasso"),
                       method     = c("weighting", "a_learning"),
                       ...)
{

    family <- match.arg(family)
    loss   <- match.arg(loss)
    method <- match.arg(method)

    this.call <- match.call()

    trt         <- as.integer(trt)
    unique.trts <- sort(unique(trt))

    if (length(unique.trts) != 2) stop("trt must have 2 distinct levels")
    if (any(unique.trts != c(0, 1))) stop("trt should be coded as 0 and 1")

    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("pi.x should be between 0 and 1")

    trt2 <- 2 * trt - 1

    if (method == "weighting")
    {
        x.tilde <- trt2 * cbind(1, x)
        wts     <- pi.x * (trt == 1) + (1 - pi.x) * (trt == 0)
    } else
    {
        x.tilde <- (trt - pi.x) * cbind(1, x)
        wts     <- rep(1, nrow(x))
    }

    fit_fun      <- paste0("fit_", loss)
    fitted.model <- do.call(fit_fun, list(x.tilde, y, wts, family, ...))

    fitted.model$loss   <- loss
    fitted.model$method <- method
    fitted.model$benefit.scores <- fitted.model$predict(x)

    fitted.model
}
