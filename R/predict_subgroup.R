
#' Predict function for fitted subgroup identification models
#'
#' @description Predicts benefit score based on a fitted subgroup identification model
#'
#' @param object fitted object returned by \code{validate.subgrp()} function
#' @param newx new design matrix for which predictions will be made
#' @param type type of prediction. \code{"benefit.score"} results in predicted benefit scores and
#' \code{"trt.group"} results in prediction of recommended treatment group
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group
#' @param ... not used
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname predict
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 1000
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
#' subgrp.model <- fit.subgroup(x = x, y = y,
#'                             trt = trt01,
#'                             propensity.func = prop.func,
#'                             loss   = "sq_loss_lasso",
#'                             nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model$subgroup.trt.effects
#' benefit.scores <- predict(subgrp.model, newx = x, type = "benefit.score")
#'
#' rec.trt.grp <- predict(subgrp.model, newx = x, type = "trt.group")
#' @export
predict.subgroup_fitted <- function(object,
                                    newx,
                                    type     = c("benefit.score", "trt.group"),
                                    cutpoint = 0,
                                    ...)
{
    type <- match.arg(type)

    retval <- drop(object$predict(newx))

    if (type == "trt.group")
    {
        if (object$larger.outcome.better)
        {
            retval <- 1 * (retval > cutpoint)
        } else
        {
            retval <- 1 * (retval < cutpoint)
        }
    }
    retval
}
