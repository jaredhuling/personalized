
#' Predict function for fitted subgroup identification models
#'
#' @description Predicts benefit score based on a fitted subgroup identification model
#'
#' @param object fitted object returned by \code{validate.subgrp()} function.
#'
#' For \code{predict.wksvm()}, this should be a fitted \code{wksvm} object from the \code{weighted.ksvm()} function
#' @param newx new design matrix for which predictions will be made
#' @param type type of prediction. \code{type = "benefit.score"} results in predicted benefit scores and
#' \code{type = "trt.group"} results in prediction of recommended treatment group.
#'
#' For \code{predict.wksvm()}, \code{type = 'class'} yields predicted
#' class and \code{type = 'linear.predictor'} yields estimated function (the sign of which is the estimated class)
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group. Can also set \code{cutpoint = "median"}, which will
#' use the median value of the benefit scores as the cutpoint or can set specific quantile values via \code{"quantx"}
#' where \code{"x"} is a number between 0 and 100 representing the quantile value; e.g. \code{cutpoint = "quant75"}
#' will use the 75th perent upper quantile of the benefit scores as the quantile.
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

    # simply call prediction function
    # defined by the loss function used
    if (grepl("owl_", object$loss) & object$n.trts > 2 & type == "trt.group")
    {
        retval <- drop(object$predict(newx, type = "class"))
    } else
    {
        retval <- drop(object$predict(newx))
    }

    cutpoint <- convert.cutpoint(cutpoint, retval)

    # need to make predicted (ie recommended)
    # treatment behavior different if larger
    # outcomes are better
    if (type == "trt.group")
    {
        if (object$n.trts > 2)
        {
            if (grepl("owl_", object$loss))
            {
                # nothing to be done
            } else
            {
                # meaning of larger vs smaller benefit score
                # is different depending on whether larger means
                # better or not for the outcome
                if (object$larger.outcome.better)
                {
                    best.comp.idx   <- apply(retval, 1, which.max)
                    recommended.trt <- 1 * (retval > cutpoint)
                    rec.ref         <- rowSums(recommended.trt) == 0

                    retval <- ifelse(rec.ref, object$reference.trt, object$comparison.trts[best.comp.idx])
                } else
                {
                    best.comp.idx   <- apply(retval, 1, which.min)
                    recommended.trt <- 1 * (retval < cutpoint)
                    rec.ref         <- rowSums(recommended.trt) == 0

                    retval <- ifelse(rec.ref, object$reference.trt, object$comparison.trts[best.comp.idx])
                }
            }

        } else
        {
            # meaning of larger vs smaller benefit score
            # is different depending on whether larger means
            # better or not for the outcome
            if (object$larger.outcome.better)
            {
                retval <- ifelse(retval > cutpoint, object$comparison.trts, object$reference.trt)
            } else
            {
                retval <- ifelse(retval < cutpoint, object$comparison.trts, object$reference.trt)
            }

        }

    }
    retval
}
