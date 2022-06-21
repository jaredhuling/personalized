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
                "sq_loss_xgboost",
                "custom")

calc_treatment_effects <- function(benefit.scores,
                                   loss = c("sq_loss_lasso",
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
                                            "sq_loss_xgboost",
                                            "custom"),
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

    effects <- list(delta = trt_eff_delta,
                    gamma = trt_eff_gamma)
    class(effects) <- c("individual_treatment_effects", class(effects) )

    effects
}

#' @rdname treatment.effects
#' @export
treatment.effects.default <- function(x, ...)
{
    treatment.effects.subgroup_fitted(x, ...)
}


#' @param benefit.scores vector of estimated benefit scores
#' @param loss loss choice USED TO CALCULATE \code{benefit.scores} of both the M function from Chen, et al (2017) and
#' potentially the penalty used for variable selection. See \code{\link[personalized]{fit.subgroup}} for more details.
#' @param method method choice USED TO CALCULATE \code{benefit.scores}. Either the \code{"weighting"} method or
#' \code{"a_learning"} method. See \code{\link[personalized]{fit.subgroup}} for more details
#' @param pi.x The propensity score for each observation
#' @rdname treatment.effects
#' @export
treat.effects <- function(benefit.scores,
                          loss = c("sq_loss_lasso",
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
                                   "sq_loss_xgboost",
                                   "custom"),
                          method = c("weighting", "a_learning"),
                          pi.x = NULL, ...)
{
    loss   <- match.arg(loss)
    method <- match.arg(method)
    calc_treatment_effects(benefit.scores, loss, method, pi.x)
}



#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @param x a fitted object from \code{fit.subgroup()}
#' @seealso \code{\link[personalized]{print.individual_treatment_effects}} for printing of objects returned by
#' \code{treat.effects} or \code{treatment.effects}
#' @rdname treatment.effects
#' @export
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 25
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,11]
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
#'                              nfolds = 3)    # option for cv.glmnet
#'
#' trt_eff <- treatment.effects(subgrp.model)
#' str(trt_eff)
#'
#' trt_eff
#'
#'
#' library(survival)
#' subgrp.model.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "cox_loss_lasso",
#'                            nfolds = 3)              # option for cv.glmnet
#'
#' trt_eff_c <- treatment.effects(subgrp.model.cox)
#' str(trt_eff_c)
#'
#' trt_eff_c
#'
treatment.effects.subgroup_fitted <- function(x, ...)
{
    treat.effects(x$benefit.scores,
                  x$loss,
                  x$method,
                  x$pi.x)
}







#' Printing individualized treatment effects
#'
#' @description Prints results for estimated subgroup treatment effects
#'
#' @param x a fitted object from either \code{\link[personalized]{treat.effects}} or \code{\link[personalized]{treatment.effects}}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @export
print.individual_treatment_effects <- function(x, digits = max(getOption('digits')-3, 3), ...)
{

    if (!is.na(x$delta[1]))
    {
        comp <- attr(x$delta, "comparison.trts")
        ref  <- attr(x$delta, "reference.trt")
        trts <- attr(x$delta, "trts")
        if (NCOL(x$delta) == 1)
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", comp, ", X] - E[Y|T=", ref, ", X]\n\n"))
            print(summary(x$delta), digits = digits, ...)
        } else if (NCOL(x$delta) == length(trts))  ##multinomial model case
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", "trt", ", X]\nwhere 'trt' is ", paste(trts, collapse = " and "), "\n\n"))

            print(summary(x$delta), digits = digits, ...)
        } else
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", "trt", ", X] - E[Y|T=", ref, ", X]\nwhere 'trt' is ", paste(comp, collapse = " and "), "\n\n"))

            print(summary(x$delta), digits = digits, ...)
        }
    }


    if (!is.na(x$gamma[1]))
    {
        comp <- attr(x$gamma, "comparison.trts")
        ref  <- attr(x$gamma, "reference.trt")
        trts <- attr(x$gamma, "trts")
        if (NCOL(x$gamma) == 1)
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", comp, ", X] / E[Y|T=", ref, ", X]\n\n"))

            cat(paste0("Note: for survival outcomes, the above ratio is \nE[g(Y)|T=", comp, ", X] / E[g(Y)|T=", ref, ", X], \nwhere g() is a monotone increasing function of Y, \nthe survival time\n\n"))

            print(summary(x$gamma), digits = digits, ...)
        } else if (NCOL(x$gamma) == length(trts))  ##multinomial model case
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", "trt", ", X]\nwhere 'trt' is ", paste(trts, collapse = " and "), "\n\n"))

            print(summary(x$gamma), digits = digits, ...)
        } else
        {
            cat(paste0("Summary of individual treatment effects: \nE[Y|T=", "trt", ", X] / E[Y|T=", ref, ", X]\nwhere 'trt' is ", paste(comp, collapse = " and "), "\n\n"))

            cat("Note: for survival outcomes, the above ratio is \nE[g(Y)|T=1, X] / E[g(Y)|T=-1, X], \nwhere g() is a monotone increasing function of Y, \nthe survival time\n\n")

            print(summary(x$gamma), digits = digits, ...)
        }
    }

}

