
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
#' @param loss choice of both the M function from Chen, et al (2017) and potentially the penalty used for variable selection.
#' All \code{loss} options starting with \code{sq_loss} use M(y, v) = (v - y) ^ 2, all options starting with \code{logistic_loss} use
#' the logistic loss: M(y, v) = y * log(1 + exp\{-v\}), and all options starting with \code{cox_loss} use the negative partial likelihood loss for the Cox PH model.
#' All options ending with \code{lasso} have a lasso penalty added to the loss for variable selection. \code{sq_loss_lasso_gam}
#' and \code{logistic_loss_lasso_gam} first use the lasso to select variables and then fit a generalized additive model
#' with nonparametric additive terms for each selected variable. \code{sq_loss_gam} involves a squared error loss with a generalized additive model and no variable selection.
#' \code{sq_loss_gbm} involves a squared error loss with a gradient-boosted decision trees model for the benefit score; this
#' allows for flexible estimation using machine learning and can be useful when the underlying treatment-covariate interaction
#' is complex.
#' \itemize{
#'     \item{Continuous Outcomes}
#'     \itemize{
#'         \item{\code{"sq_loss_lasso"}}{ - M(y, v) = (v - y) ^ 2 with linear model and lasso penalty}
#'         \item{\code{"owl_logistic_loss_lasso"}} { - M(y, v) = ylog(1 + exp\{-v\}) (method of Regularized Outcome Weighted Subgroup Identification)}
#'         \item{\code{"owl_logistic_flip_loss_lasso"}} { - M(y, v) = |y|log(1 + exp\{-sign(y)v\})}
#'         \item{\code{"owl_hinge_loss"}} { - M(y, v) = ymax(0, 1 - v) (method of Estimating individualized treatment rules using outcome weighted learning)}
#'         \item{\code{"owl_hinge_flip_loss"}} { - M(y, v) = |y|max(0, 1 - sign(y)v) }
#'         \item{\code{"sq_loss_lasso_gam"}}{ - M(y, v) = (v - y) ^ 2 with variables selected by lasso penalty and generalized additive model fit on the selected variables}
#'         \item{\code{"sq_loss_gam"}}{ - M(y, v) = (v - y) ^ 2 with generalized additive model fit on all variables}
#'         \item{\code{"owl_logistic_loss_gam"}} { - M(y, v) = ylog(1 + exp\{-v\}) with generalized additive model fit on all variables}
#'         \item{\code{"owl_logistic_flip_loss_gam"}} { - M(y, v) = |y|log(1 + exp\{-sign(y)v\}) with generalized additive model fit on all variables}
#'         \item{\code{"owl_logistic_loss_lasso_gam"}} { - M(y, v) = ylog(1 + exp\{-v\}) with variables selected by lasso penalty and generalized additive model fit on the selected variables}
#'         \item{\code{"owl_logistic_flip_loss_lasso_gam"}} { - M(y, v) = |y|log(1 + exp\{-sign(y)v\}) with variables selected by lasso penalty and generalized additive model fit on the selected variables}
#'         \item{\code{"sq_loss_gbm"}}{ - M(y, v) = (v - y) ^ 2 with gradient-boosted decision trees model}
#'         \item{\code{"abs_loss_gbm"}}{ - M(y, v) = |v - y| with gradient-boosted decision trees model}
#'     }
#'     \item{Binary Outcomes}
#'     \itemize{
#'         \item{All losses for continuous outcomes can be used plus the following:}
#'         \item{\code{"logistic_loss_lasso"}}{ - M(y, v) = -[yv - log(1 + exp\{-v\})] with with linear model and lasso penalty}
#'         \item{\code{"logistic_loss_lasso_gam"}}{ - M(y, v) = y * log(1 + exp\{-v\}) with variables selected by lasso penalty and generalized additive model fit on the selected variables}
#'         \item{\code{"logistic_loss_gam"}}{ - M(y, v) = y * log(1 + exp\{-v\}) with generalized additive model fit on all variables}
#'         \item{\code{"logistic_loss_gbm"}}{ - M(y, v) = -[yv - log(1 + exp\{-v\})] with gradient-boosted decision trees model}
#'     }
#'     \item{Time-to-Event Outcomes}
#'     \itemize{
#'         \item{\code{"cox_loss_lasso"}}{ - M corresponds to the negative partial likelihood of the cox model with linear model and additionally a lasso penalty}
#'         \item{\code{"cox_loss_gbm"}}{ - M corresponds to the negative partial likelihood of the cox model with gradient-boosted decision trees model}
#'     }
#' }
#' @param method subgroup ID model type. Either the weighting or A-learning method of Chen et al, (2017)
#' @param match.id a (character, factor, or integer) vector with length equal to the number of observations in \code{x} indicating using integers or
#' levels of a factor vector which patients are
#' in which matched groups. Defaults to \code{NULL} and assumes the samples are not from a matched cohort. Matched
#' case-control groups can be created using any method (propensity score matching, optimal matching, etc). If each case
#' is matched with a control or multiple controls, this would indicate which case-control pairs or groups go together.
#' If \code{match.id} is supplied, then it is unecessary to specify a function via the \code{propensity.func} argument.
#' A quick usage example: if the first patient is a case and the second and third are controls matched to it, and the
#' fouth patient is a case and the fifth through seventh patients are matched with it, then the user should specify
#' \code{match.id = c(1,1,1,2,2,2,2)} or \code{match.id = c(rep("Grp1", 3),rep("Grp2", 4)) }
#' @param augment.func function which inputs the response \code{y}, the covariates \code{x}, and \code{trt} and outputs
#' predicted values (on the link scale) for the response using a model constructed with \code{x}. \code{augment.func()} can also be simply
#' a function of \code{x} and \code{y}. This function is used for efficiency augmentation.
#' When the form of the augmentation function is correct, it can provide efficient estimation of the subgroups
#' Example 1: \code{augment.func <- function(x, y) {lmod <- lm(y ~ x); return(fitted(lmod))}}
#'
#' Example 2: \code{augment.func <- function(x, y, trt) {lmod <- lm(y ~ x + trt); return(fitted(lmod))}}
#'
#' For binary and time-to-event outcomes, make sure that predictions are returned on the scale of the predictors
#'
#' Example 3:
#'
#' \code{augment.func <- function(x, y) {
#'     bmod <- glm(y ~ x, family = binomial());
#'     return(predict(bmod, type = "link"))}
#'  }
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group. Can also set \code{cutpoint = "median"}, which will
#' use the median value of the benefit scores as the cutpoint or can set specific quantile values via \code{"quantx"}
#' where \code{"x"} is a number between 0 and 100 representing the quantile value; e.g. \code{cutpoint = "quant75"}
#' will use the 75th perent upper quantile of the benefit scores as the quantile.
#' @param larger.outcome.better boolean value of whether a larger outcome is better/preferable. Set to \code{TRUE}
#' if a larger outcome is better/preferable and set to \code{FALSE} if a smaller outcome is better/preferable. Defaults to \code{TRUE}.
#' @param reference.trt which treatment should be treated as the reference treatment. Defaults to the first level of \code{trt}
#' if \code{trt} is a factor or the first alphabetical or numerically first treatment level. Not used for multiple treatment fitting with OWL-type losses.
#' @param retcall boolean value. if \code{TRUE} then the passed arguments will be saved. Do not set to \code{FALSE}
#' if the \code{validate.subgroup()} function will later be used for your fitted subgroup model. Only set to \code{FALSE}
#' if memory is limited as setting to \code{TRUE} saves the design matrix to the fitted object
#' @param ... options to be passed to underlying fitting function. For all \code{loss} options with 'lasso',
#' this will be passed to \code{\link[glmnet]{cv.glmnet}}. For all \code{loss} options with 'gam', this will
#' be passed to \code{\link[mgcv]{gam}} from the \pkg{mgcv} package
#' Note that for all \code{loss} options that use \code{gam()}
#' from the \pkg{mgcv} package,
#' the user cannot supply the \code{gam} argument \code{method} because it is also an argument of \code{fit.subgroup}, so
#' instead, to change the \code{gam method} argument, supply \code{method.gam}, ie \code{method.gam = "REML"}.
#'
#' For all \code{loss} options with 'hinge', this will be passed to both \code{\link[personalized]{weighted.ksvm}} and
#' \code{\link[kernlab]{ipop}} from the \pkg{kernlab} package
#' @seealso \code{\link[personalized]{validate.subgroup}} for function which creates validation results for subgroup
#' identification models, \code{\link[personalized]{predict.subgroup_fitted}} for a prediction function for fitted models
#' from \code{fit.subgroup}, \code{\link[personalized]{plot.subgroup_fitted}} for a function which plots
#' results from fitted models, and \code{\link[personalized]{print.subgroup_fitted}}
#' for arguments for printing options for \code{fit.subgroup()}.
#' from \code{fit.subgroup}.
#' @references Chen, S., Tian, L., Cai, T. and Yu, M. (2017), A general statistical framework for subgroup identification
#' and comparative treatment scoring. Biometrics. doi:10.1111/biom.12676 \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12676/abstract}
#'
#' Xu, Y., Yu, M., Zhao, Y. Q., Li, Q., Wang, S., & Shao, J. (2015),
#'  Regularized outcome weighted subgroup identification for differential treatment effects. Biometrics, 71(3), 645-653.
#'  doi: 10.1111/biom.12322 \url{http://onlinelibrary.wiley.com/doi/10.1111/biom.12322/full}
#'
#'  Zhao, Y., Zeng, D., Rush, A. J., & Kosorok, M. R. (2012),
#'   Estimating individualized treatment rules using outcome weighted learning.
#'   Journal of the American Statistical Association, 107(499), 1106-1118. doi: 10.1080/01621459.2012.695674
#'   \url{http://dx.doi.org/10.1080/01621459.2012.695674}
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
#' delta <- 2 * (0.5 + x[,2] - x[,3] - x[,11] + x[,1] * x[,12] )
#' xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
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
#' subgrp.model <- fit.subgroup(x = x, y = y,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "sq_loss_lasso",
#'                            nfolds = 10)              # option for cv.glmnet
#'
#' summary(subgrp.model)
#'
#' # fit lasso + gam model with REML option for gam
#'
#' subgrp.modelg <- fit.subgroup(x = x, y = y,
#'                             trt = trt01,
#'                             propensity.func = prop.func,
#'                             loss   = "sq_loss_lasso_gam",
#'                             method.gam = "REML",     # option for gam
#'                             nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.modelg
#'
#' subgrp.model.bin <- fit.subgroup(x = x, y = y.binary,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "logistic_loss_lasso",
#'                            type.measure = "auc",    # option for cv.glmnet
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model.bin
#'
#' library(survival)
#' subgrp.model.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "cox_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model.cox
#'
#'
#' @export
fit.subgroup <- function(x,
                         y,
                         trt,
                         propensity.func = NULL,
                         loss       = c("sq_loss_lasso",
                                        "logistic_loss_lasso",
                                        "cox_loss_lasso",
                                        "owl_logistic_loss_lasso",
                                        "owl_logistic_flip_loss_lasso",
                                        "owl_hinge_loss",
                                        "owl_hinge_flip_loss",
                                        "sq_loss_lasso_gam",
                                        "logistic_loss_lasso_gam",
                                        "sq_loss_gam",
                                        "logistic_loss_gam",
                                        "owl_logistic_loss_gam",
                                        "owl_logistic_flip_loss_gam",
                                        "owl_logistic_loss_lasso_gam",
                                        "owl_logistic_flip_loss_lasso_gam",
                                        "sq_loss_gbm",
                                        "abs_loss_gbm",
                                        "logistic_loss_gbm",
                                        "cox_loss_gbm"),
                         method     = c("weighting", "a_learning"),
                         match.id = NULL,
                         augment.func = NULL,
                         cutpoint   = 0,
                         larger.outcome.better = TRUE,
                         reference.trt = NULL,
                         retcall    = TRUE,
                         ...)
{

    loss   <- match.arg(loss)
    method <- match.arg(method)


    # make sure outcome is consistent with
    # other options selected if there is any
    # indication of a survival outcome or model
    if ( xor(class(y) == "Surv", grepl("cox_loss", loss)) )
    {
        ifelse(
            grepl("cox_loss", loss),
            stop("Must provide 'Surv' object if loss/family corresponds to a Cox model. See\n
                 '?Surv' for more information about 'Surv' objects."),
            stop("Loss and family must correspond to a Cox model for time-to-event outcomes.")
            )
    }

    if (grepl("cox_loss", loss))
    {
        family <- "cox"
    } else if (grepl("logistic_loss", loss) | grepl("huberized_loss", loss))
    {
        family <- "binomial"
    } else
    {
        family <- "gaussian"
    }


    dims   <- dim(x)
    if (is.null(dims)) stop("x must be a matrix object.")

    y      <- drop(y)
    vnames <- colnames(x)

    if (!is.null(match.id))
    {
        if (length(match.id) != dims[1L]) stop("match.id must be same length as number of observations.")
    }

    # set variable names if they are not set
    if (is.null(vnames)) vnames <- paste0("V", 1:dims[2])

    ## will be a flag for later use if
    ## I decide to make outcome-weighted learning
    ## (ie flipping loss) an option
    outcome.weighted <- FALSE





    # check to make sure arguments of augment.func are correct
    if (!is.null(augment.func))
    {
        augmentfunc.names <- sort(names(formals(augment.func)))
        if (length(augmentfunc.names) == 3)
        {
            if (any(augmentfunc.names != c("trt", "x", "y")))
            {
                stop("arguments of augment.func() should be 'trt', 'x', and 'y'")
            }
        } else if (length(augmentfunc.names) == 2)
        {
            if (any(augmentfunc.names != c("x", "y")))
            {
                stop("arguments of augment.func() should be 'x' and 'y'")
            }
            augment.func2 <- augment.func
            augment.func  <- function(trt, x, y) augment.func2(x = x, y = y)
        } else
        {
            stop("augment.func() should only have either two arguments: 'x' and 'y', or three arguments:
                 'trt', 'x', and 'y'")
        }
    }

    ## this can be done with 'grepl("owl_", loss)' easily
    ## but is left as below in case there are more loss options
    ## in the future
    augment.method <- switch(loss,
                             "sq_loss_lasso"                    = "offset",
                             "logistic_loss_lasso"              = "offset",
                             "cox_loss_lasso"                   = "offset",
                             "owl_logistic_loss_lasso"          = "adj",
                             "owl_logistic_flip_loss_lasso"     = "adj",
                             "owl_hinge_loss"                   = "adj",
                             "owl_hinge_flip_loss"              = "adj",
                             "sq_loss_lasso_gam"                = "offset",
                             "logistic_loss_lasso_gam"          = "offset",
                             "sq_loss_gam"                      = "offset",
                             "logistic_loss_gam"                = "offset",
                             "owl_logistic_loss_gam"            = "adj",
                             "owl_logistic_flip_loss_gam"       = "adj",
                             "owl_logistic_loss_lasso_gam"      = "adj",
                             "owl_logistic_flip_loss_lasso_gam" = "adj",
                             "sq_loss_gbm"                      = "offset",
                             "abs_loss_gbm"                     = "offset",
                             "logistic_loss_gbm"                = "offset",
                             "cox_loss_gbm"                     = "offset")

    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    if (n.trts > 2 & grepl("owl_", loss) & grepl("hinge_", loss))
    {
        stop("OWL hinge loss not available for multiple treatments.")
    }

    if (grepl("owl_", loss))
    {
        if (method != "weighting")
        {
            warning("Only method = 'weighting' available for OWL-type losses; defaulting
                    to 'weighting' method.")
            method <- "weighting"
        }
    }

    if (n.trts < 2)           stop("trt must have at least 2 distinct levels")
    if (n.trts > dims[1] / 3) stop("trt must have no more than n.obs / 3 distinct levels")

    if (!is.null(reference.trt))
    {
        if (grepl("owl_", loss))
        {
            warning("reference.trt specification not applicable for owl losses.")
            reference.idx   <- 1L
            comparison.idx  <- (1:n.trts)[-reference.idx]
            comparison.trts <- unique.trts[-reference.idx]
            reference.trt   <- unique.trts[reference.idx]
        } else
        {
            if (!(reference.trt %in% unique.trts))
            {
                stop("reference.trt must be one of the treatment levels")
            }

            reference.idx   <- which(unique.trts == reference.trt)
            comparison.idx  <- (1:n.trts)[-reference.idx]
            comparison.trts <- unique.trts[-reference.idx]
        }
        refnull <- FALSE
    } else
    {
        reference.idx   <- 1L
        comparison.idx  <- (1:n.trts)[-reference.idx]
        comparison.trts <- unique.trts[-reference.idx]
        reference.trt   <- unique.trts[reference.idx]
        refnull <- TRUE
    }

    if (n.trts > 2 & (grepl("_gbm", loss) | grepl("_gam", loss)) )
    {
        stop("gbm and gam based losses not supported for multiple treatments (number of total treatments > 2)")
    }

    # Check match.id validity and convert it to a factor, if supplied
    if (!is.null(match.id)) {
        match.id <- tryCatch(expr=as.factor(match.id), error = function(e) {stop("match.id must be a factor or capable of being coerced to a factor.")})
        if (length(levels(match.id)) < 2) {stop("match.id must have at least 2 levels")}
    }

    # defaults to constant propensity score within trt levels
    # the user will almost certainly want to change this
    if (is.null(propensity.func))
    {
        if (is.null(match.id))
        { # No propensity score supplied and no match.id supplied
            if (n.trts == 2)
            {
                mean.trt <- mean(trt == unique.trts[2L])
                propensity.func <- function(trt, x) rep(mean.trt, length(trt))
            } else
            {
                mean.trt <- numeric(n.trts)
                for (t in 1:n.trts)
                {
                    mean.trt[t] <- mean(trt == unique.trts[t])
                }
                propensity.func <- function(trt, x)
                {
                    pi.x <- numeric(length(trt))
                    for (t in 1:n.trts)
                    {
                        which.t       <- trt == unique.trts[t]
                        pi.x[which.t] <- mean(which.t)
                    }
                    pi.x
                }
            }
        } else
        { # No propensity score supplied but match.id supplied
            if (n.trts == 2)
            {
                # default to pct in treatment group
                mean.trt <- mean(trt == unique.trts[2L])
                pf <- function(trt, x, match.id) rep(mean.trt, NROW(trt))
                propensity.func <- pf
            } else
            {
                mean.trt <- numeric(n.trts)
                for (t in 1:n.trts)
                {
                    mean.trt[t] <- mean(trt == unique.trts[t])
                }
                pf <- function(trt, x, match.id)
                {
                    pi.x <- numeric(length(trt))
                    for (t in 1:n.trts)
                    {
                        which.t       <- trt == unique.trts[t]
                        pi.x[which.t] <- mean(which.t)
                    }
                    pi.x
                }
                propensity.func <- pf
            }
        }
    }


    extra.args <- NULL
    # check to make sure arguments of augment.func are correct
    if (!is.null(augment.func))
    {
        B.x   <- unname(drop(augment.func(trt = trt, x = x, y = y)))

        if (NROW(B.x) != NROW(y))
        {
            stop("augment.func() should return the same number of predictions as observations in y")
        }

        if (augment.method == "adj")
        {
            y.adj <- y - B.x
        } else
        {
            y.adj      <- y
            extra.args <- list(offset = B.x)
        }

    } else
    {
        y.adj <- y
    }

    # stop if augmentation function provided
    # for non-gaussian outcomes.
    # has not been developed yet
    #if (!is.null(augment.func) & family != "gaussian")
    #{
    #    warning("Efficiency augmentation not available for non-continuous outcomes yet. No augmentation applied.")
    #}

    larger.outcome.better <- as.logical(larger.outcome.better[1])
    retcall               <- as.logical(retcall[1])

    # save the passed arguments for later use in validate.subgroupu()
    # and plot.subgroup_fitted() functions
    if (retcall)
    {
        this.call     <- mget(names(formals()), sys.frame(sys.nframe()))

        this.call$... <- NULL
        this.call     <- c(this.call, list(...))
    } else
    {
        this.call     <- NULL
    }

    # check to make sure arguments of propensity.func are correct
    propfunc.names <- sort(names(formals(propensity.func)))
    if (length(propfunc.names) == 3)
    {
        if (any(propfunc.names != c("match.id", "trt", "x")))
        {
            stop("arguments of propensity.func() should be 'trt','x', and (optionally) 'match.id'")
        }
    } else if (length(propfunc.names) == 2)
    {
        if (any(propfunc.names != c("trt", "x")))
        {
            stop("arguments of propensity.func() should be 'trt','x', and (optionally) 'match.id'")
        }
    } else
    {
        stop("propensity.func() should only have two or three arguments: 'trt' and 'x', or: 'trt', 'x', and 'match.id'")
    }

    # compute propensity scores
    if (is.null(match.id) | length(propfunc.names) == 2)
    {
        pi.x <- drop(propensity.func(x = x, trt = trt))
    } else
    {
        pi.x <- drop(propensity.func(x = x, trt = trt, match.id = match.id))
    }


    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity.func() should return values between 0 and 1")

    ## if returned propensity score
    ## is a matrix, then pick out the
    ## right column for each row so we
    ## always get Pr(T = T_i | X = x)

    dim.pi.x <- dim(pi.x)
    if (!is.null(dim.pi.x))
    {
        if (length(dim.pi.x) == 1)
        {
            pi.x <- as.vector(pi.x)
        } else if (length(dim.pi.x) == 2)
        {
            if (ncol(pi.x) != n.trts)
            {
                stop("Number of columns in the matrix returned by propensity.func() is not the same
                     as the number of levels of 'trt'.")
            }
            if (is.factor(trt))
            {
                values <- levels(trt)[trt]
            } else
            {
                values <- trt
            }

            levels.pi.mat <- colnames(pi.x)
            if (is.null(levels.pi.mat))
            {
                levels.pi.mat <- unique.trts
            }

            # return the probability corresponding to the
            # treatment that was observed
            pi.x <- pi.x[cbind(1:nrow(pi.x), match(values, levels.pi.mat))]
            } else
            {
                stop("propensity.func() returns a multidimensional array; it can only return a vector or matrix.")
            }
    }

    if (grepl("owl_", loss))
    {
        x.tilde <- cbind(1, x)
    } else
    {
        # construct design matrix to be passed to fitting function
        x.tilde <- create.design.matrix(x             = x,
                                        pi.x          = pi.x,
                                        trt           = trt,
                                        method        = method,
                                        reference.trt = reference.trt)
    }

    # construct observation weight vector
    wts     <- create.weights(pi.x   = pi.x,
                              trt    = trt,
                              method = method)


    if (n.trts > 2 & !grepl("owl_", loss))
    {
        all.cnames <- numeric(ncol(x.tilde))
        len.names  <- length(vnames) + 1
        for (tr in 1:(n.trts - 1))
        {
            idx.cur <- ((len.names * (tr - 1)) + 1):(len.names * tr)
            if (comparison.trts[tr] != 1 & comparison.trts[tr] != "1")
            {
                trt.name.cur <- comparison.trts[tr]
            } else
            {
                trt.name.cur <- "Trt1"
            }
            all.cnames[idx.cur] <- c( trt.name.cur,
                                      paste(vnames, 1:(n.trts - 1), sep = ".") )
        }
    } else
    {
        if (comparison.trts[1] != 1 & comparison.trts[1] != "1")
        {
            trt.name.cur <- comparison.trts[1]
        } else
        {
            trt.name.cur <- "Trt1"
        }
        all.cnames <- c( trt.name.cur[1],
                         vnames )

        if (grepl("owl_", loss) & n.trts > 2)
        {
            all.cnames[1] <- "Trt"
        }
    }


    colnames(x.tilde) <- all.cnames


    # extra preparation needed for OWL-type losses
    if (grepl("owl_", loss))
    {
        intercept <- FALSE
        if (n.trts == 2)
        {
            family <- "binomial"
        } else
        {
            family <- "multinomial"
        }


        if (grepl("logistic_", loss) & grepl("lasso", loss) & !grepl("gam$", loss))
        {
            fit_fun <- "fit_logistic_loss_lasso"

            intercept <- TRUE
        } else if (grepl("logistic_", loss) & grepl("gam$", loss))
        {
            if (grepl("lasso", loss))
            {
                fit_fun <- "fit_logistic_loss_lasso_gam"
            } else
            {
                fit_fun <- "fit_logistic_loss_gam"
            }
        } else if (grepl("hinge_", loss))
        {
            fit_fun <- "fit_owl_hinge_loss"
        }


        if (grepl("flip_", loss))
        {
            if (n.trts == 2)
            {
                trt.y <- trt
                y.wt  <- abs(y.adj)

                # flip the trtmnt for negative response values
                idx.flip.1 <- y.adj < 0 & trt == unique.trts[1]
                idx.flip.2 <- y.adj < 0 & trt == unique.trts[2]
                trt.y[idx.flip.1] <- unique.trts[2]
                trt.y[idx.flip.2] <- unique.trts[1]
            } else
            {
                stop("Flipping loss not available for multiple treatments scenarios.")
            }
        } else
        {
            trt.y <- trt
            y.wt  <- y.adj - min(y.adj)
        }

        fitted.model <- do.call(fit_fun, c(list(x = x.tilde, trt = trt, n.trts = n.trts,
                                                y = trt.y, wts = drop(wts) * drop(y.wt),
                                                family = family, intercept = intercept,
                                                match.id = match.id, ...), extra.args) )
    } else
    {
        # identify correct fitting function and call it
        fit_fun      <- paste0("fit_", loss)

        fitted.model <- do.call(fit_fun, c(list(x = x.tilde, trt = trt, n.trts = n.trts,
                                                y = y.adj, wts = wts, family = family,
                                                match.id = match.id, ...), extra.args)  )
    }


    if (refnull) this.call$reference.trt <- NULL


    # save extra results
    fitted.model$call                  <- this.call
    fitted.model$family                <- family
    fitted.model$loss                  <- loss
    fitted.model$method                <- method
    fitted.model$propensity.func       <- propensity.func
    fitted.model$augment.func          <- augment.func
    fitted.model$larger.outcome.better <- larger.outcome.better
    fitted.model$var.names             <- vnames
    fitted.model$n.trts                <- n.trts
    fitted.model$comparison.trts       <- comparison.trts
    fitted.model$reference.trt         <- reference.trt
    fitted.model$trts                  <- unique.trts
    fitted.model$trt.received          <- trt
    fitted.model$pi.x                  <- pi.x
    fitted.model$y                     <- y

    fitted.model$benefit.scores        <- fitted.model$predict(x)

    fitted.model$recommended.trts      <- predict.subgroup_fitted(fitted.model, newx = x,
                                                                  type = "trt.group",
                                                                  cutpoint = cutpoint)

    # calculate sizes of subgroups and the
    # subgroup treatment effects based on the
    # benefit scores and specified benefit score cutpoint
    fitted.model$subgroup.trt.effects <- subgroup.effects(fitted.model$benefit.scores,
                                                          y, trt,
                                                          pi.x,
                                                          cutpoint,
                                                          larger.outcome.better,
                                                          reference.trt = reference.trt)

    class(fitted.model) <- "subgroup_fitted"

    fitted.model
}
