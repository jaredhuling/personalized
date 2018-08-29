#' Calculation of covariate-conditional treatment effects
#'
#' @description Calculates covariate conditional treatment effects using estimated benefit scores
#' @param ... not used
#' @return A List with elements \code{delta} (if the treatment effects are a difference/contrast,
#' i.e. \eqn{E[Y|T=1, X] - E[Y|T=-1, X]}) and \code{gamma} (if the treatment effects are a ratio,
#' i.e. \eqn{E[Y|T=1, X] / E[Y|T=-1, X]})
#' @export
treatment.effects <- function(x, ...) UseMethod("treatment.effects")


all_losses <- c("sq_loss_lasso",
                "logistic_loss_lasso",
                "poisson_loss_lasso",
                "cox_loss_lasso",
                "owl_logistic_loss_lasso",
                "owl_logistic_flip_loss_lasso",
                "owl_hinge_loss",
                "owl_hinge_flip_loss",
                "sq_loss_lasso_gam",
                "poisson_loss_lasso_gam",
                "logistic_loss_lasso_gam",
                "sq_loss_gam",
                "poisson_loss_gam",
                "logistic_loss_gam",
                "owl_logistic_loss_gam",
                "owl_logistic_flip_loss_gam",
                "owl_logistic_loss_lasso_gam",
                "owl_logistic_flip_loss_lasso_gam",
                "sq_loss_gbm",
                "poisson_loss_gbm",
                "logistic_loss_gbm",
                "cox_loss_gbm",
                "custom")

calc_treatment_effects <- function(benefit.scores,
                                   loss = all_losses,
                                   method = c("weighting", "a_learning"),
                                   pi.x = NULL)
{
    loss   <- match.arg(loss)
    method <- match.arg(method)

    benefit.scores <- drop(benefit.scores)
    if (!is.null(pi.x))
    {
        pi.x <- drop(pi.x)
    }

    if (method == "a_learning" & is.null(pi.x) &
        !(grepl("sq_loss", loss) | grepl("cox_loss", loss)) )
    {
        stop("propensity score vector 'pi.x' must be supplied for method = 'a_learning'")
    }

    trt_eff_delta <- trt_eff_gamma <- NA

    if (grepl("sq_loss", loss))
    {
        if (method == "weighting")
        {
            trt_eff_delta <- 2 * benefit.scores
        } else
        {
            trt_eff_delta <- benefit.scores
        }
    } else if (grepl("owl_logistic_loss", loss) )
    {
        if (method == "weighting")
        {
            trt_eff_gamma <- exp(benefit.scores)
        } else
        {
            trt_eff_gamma <- (1 + exp((1 - pi.x) * benefit.scores)) / (1 + exp(-pi.x * benefit.scores))
        }
    } else if (grepl("cox_loss", loss) )
    {
        trt_eff_gamma <- exp(-benefit.scores)

    }  else if (grepl("poisson_loss", loss) )
    {
        if (method == "weighting")
        {
            trt_eff_delta <- exp(benefit.scores) - exp(-benefit.scores)
        } else
        {
            trt_eff_delta <- exp((1 - pi.x) * benefit.scores) - exp(-pi.x * benefit.scores)
        }
    } else if (grepl("logistic_loss", loss) & !grepl("owl_", loss))
    {
        if (method == "weighting")
        {
            exp_bs <- exp(benefit.scores)
            trt_eff_delta <- (exp_bs - 1) / (exp_bs + 1)
        } else
        {
            exp_bs    <- exp(benefit.scores)
            exp_pi_bs <- exp(pi.x * benefit.scores)
            trt_eff_delta <- ((exp_bs - 1) / (exp_pi_bs + 1)) * (exp_pi_bs / (exp_pi_bs + exp_bs))
        }
    } else
    {
        warning(paste("treatment effects not available for loss:", loss) )
    }
    list(delta = trt_eff_delta,
         gamma = trt_eff_gamma)
}

#' @rdname treatment.effects
#' @export
treatment.effects.default <- function(x, ...)
{
    treatment.effects.subgroup_fitted(x, ...)
}


#' @param benefit.scores vector of estimated benefit scores
#' @param loss choice USED TO CALCULATE \code{benefit.scores} of both the M function from Chen, et al (2017) and
#' potentially the penalty used for variable selection. See \code{\link[personalized]{fit.subgroup}} for more details.
#' @param method method USED TO CALCULATE \code{benefit.scores}. Either the \code{"weighting"} method or
#' \code{"a_learning"} method. See \code{\link[personalized]{fit.subgroup}} for more details
#' @param pi.x The propensity score for each observation
#' @rdname treatment.effects
#' @export
treat.effects <- function(benefit.scores,
                          loss = all_losses,
                          method = c("weighting", "a_learning"),
                          pi.x = NULL, ...)
{
    loss   <- match.arg(loss)
    method <- match.arg(method)
    calc_treatment_effects(benefit.scores, loss, method, pi.x)
}



#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @param x a fitted object from \code{fit.subgroup()} or a matrix of covariate values
#' @rdname treatment.effects
#' @export
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
#' subgrp.model <- fit.subgroup(x = x, y = y,
#'                              trt = trt01,
#'                              propensity.func = prop.func,
#'                              loss   = "sq_loss_lasso",
#'                              nfolds = 5)    # option for cv.glmnet
#'
#' trt_eff <- treatment.effects(subgrp.model)
#' str(trt_eff)
#'
#' print(summary(trt_eff$delta))
#'
#'
#' library(survival)
#' subgrp.model.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "cox_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' trt_eff_c <- treatment.effects(subgrp.model.cox)
#' str(trt_eff_c)
#'
#' print(summary(trt_eff_c$gamma))
#'
treatment.effects.subgroup_fitted <- function(x, ...)
{
    treat.effects(x$benefit.scores,
                  x$loss,
                  x$method,
                  x$pi.x)
}
