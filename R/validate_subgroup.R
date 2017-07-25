#' Validating fitted subgroup identification models
#'
#' @description Validates subgroup treatment effects for fitted
#'  subgroup identification model class of Chen, et al (2017)
#'
#' @param model fitted model object returned by \code{fit.subgroup()} function
#' @param method validation method. \code{"boot_bias_correction"} for the bootstrap
#' bias correction method of Harrell, et al (1996) or \code{"training_test_replication"}
#' for repeated training and test splitting of the data (\code{train.fraction} should be specified
#' for this option)
#' @param B integer. number of bootstrap replications or refitting replications.
#' @param train.fraction fraction (between 0 and 1) of samples to be used for training in
#' training/test replication. Only used for \code{method = "training_test_replication"}
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models,
#' \code{\link[personalized]{plot.subgroup_validated}} for plotting of validation results, and
#' \code{\link[personalized]{print.subgroup_validated}} for arguments for printing options for \code{validate.subgroup()}.
#' @references Chen, S., Tian, L., Cai, T. and Yu, M. (2017), A general statistical framework for subgroup identification
#' and comparative treatment scoring. Biometrics. doi:10.1111/biom.12676
#'
#' Harrell, F. E., Lee, K. L., and Mark, D. B. (1996). Tutorial in biostatistics multivariable prognostic models: issues in developing models,
#' evaluating assumptions and adequacy, and measuring and reducing errors. Statistics in medicine, 15, 361-387.
#' doi:10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4
#' @importFrom stats predict sd
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 500
#' n.vars <- 20
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,11] - 0.5 * x[,13]
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
#'                              trt = trt01,
#'                              propensity.func = prop.func,
#'                              loss   = "sq_loss_lasso",
#'                              nfolds = 5)    # option for cv.glmnet
#'
#' subgrp.model$subgroup.trt.effects
#'
#' x.test <- matrix(rnorm(10 * n.obs * n.vars, sd = 3), 10 * n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat.test   <- 0.5 + 0.5 * x.test[,11] - 0.5 * x.test[,13]
#' trt.prob.test <- exp(xbetat.test) / (1 + exp(xbetat.test))
#' trt01.test    <- rbinom(10 * n.obs, 1, prob = trt.prob.test)
#'
#' trt.test      <- 2 * trt01.test - 1
#'
#' # simulate response
#' delta.test <- 2 * (0.5 + x.test[,2] - x.test[,3] - x.test[,11] + x.test[,1] * x.test[,12])
#' xbeta.test <- x.test[,1] + x.test[,11] - 2 * x.test[,12]^2 + x.test[,13]
#' xbeta.test <- xbeta.test + delta.test * trt.test
#'
#' y.test <- drop(xbeta.test) + rnorm(10 * n.obs, sd = 2)
#'
#' valmod <- validate.subgroup(subgrp.model, B = 10,
#'                             method = "training_test",
#'                             train.fraction = 0.75)
#' valmod
#'
#' bene.score.test <- subgrp.model$predict(x.test)
#'
#' mean(y.test[bene.score.test > 0 & trt01.test == 1]) -
#'        mean(y.test[bene.score.test > 0 & trt01.test == 0])
#' mean(y.test[bene.score.test <= 0 & trt01.test == 0]) -
#'        mean(y.test[bene.score.test <= 0 & trt01.test == 1])
#'
#' quantile(valmod$boot.results[[1]][,1], c(0.025, 0.975))
#' quantile(valmod$boot.results[[1]][,2], c(0.025, 0.975))
#' @export
validate.subgroup <- function(model,
                              B              = 50L,
                              method         = c("training_test_replication",
                                                 "boot_bias_correction"),
                              train.fraction = 0.5)
{
    method <- match.arg(method)


    if (class(model)[1] != "subgroup_fitted")
    {
        stop("model should be a fitted object returned by the 'fit.subgroup' function")
    }


    B <- as.integer(B[1])
    if (B <= 1) stop("B must be a strictly positive integer")

    train.fraction <- as.numeric(train.fraction[1])

    if (train.fraction >= 1 | train.fraction <= 0)
    {
        stop("train.fraction must be between 0 and 1")
    }

    if (is.null(model$call)) stop("retcall argument must be set to TRUE for fitted model
                                  to use validate.subgroup()")


    # no need to store passed arguments for the bootstrap or training/testing replications
    model$call$retcall <- FALSE

    # save data objects because they
    # will be written over by resampled versions later
    x   <- model$call$x
    trt <- model$call$trt
    y   <- model$call$y

    n.obs <- NROW(x)

    # create objects to store results
    boot.list      <- vector(mode = "list", length = length(model$subgroup.trt.effects))
    boot.list[[1]] <- array(NA, dim = c(B, length(model$subgroup.trt.effects[[1]])))
    boot.list[[2]] <- boot.list[[3]] <- array(NA, dim = c(B, dim(model$subgroup.trt.effects[[2]])))

    dimnames(boot.list[[2]]) <- dimnames(boot.list[[2]]) <-
        c(list(NULL), dimnames(model$subgroup.trt.effects$avg.outcomes))

    for (b in 1:B)
    {
        if (method == "training_test_replication")
        {
            # randomly split/partition data into training and testing sets
            train.samp.size <- floor(n.obs * train.fraction)
            samp.idx        <- sample.int(n.obs, train.samp.size, replace = FALSE)
            model$call$x    <- x[samp.idx,]
            model$call$trt  <- trt[samp.idx]

            x.test          <- x[-samp.idx,]

            # need to handle differently if outcome is a matrix
            if (is.matrix(y))
            {
                model$call$y <- y[samp.idx,]
                y.test       <- y[-samp.idx,]
            } else
            {
                model$call$y <- y[samp.idx]
                y.test       <- y[-samp.idx]
            }
            trt.test <- trt[-samp.idx]

            # fit subgroup model on training data
            mod.b    <- do.call(fit.subgroup, model$call)

            # compute benefit scores on testing data
            benefit.scores.test <- mod.b$predict(x.test)

            # estimate subgroup treatment effects on test data
            sbgrp.trt.eff.test  <- subgroup.effects(benefit.scores.test,
                                                    y.test, trt.test,
                                                    model$call$cutpoint)

            # save results
            boot.list[[1]][b,]  <- sbgrp.trt.eff.test[[1]]
            boot.list[[2]][b,,] <- sbgrp.trt.eff.test[[2]]
            boot.list[[3]][b,,] <- sbgrp.trt.eff.test[[3]]

        } else if (method == "boot_bias_correction")
        {   # bootstrap bias correction

            # take a bootstrap sample with replacement
            samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
            model$call$x   <- x[samp.idx,]

            if (is.matrix(y))
            {
                model$call$y   <- y[samp.idx,]
            } else
            {
                model$call$y   <- y[samp.idx]
            }
            model$call$trt <- trt[samp.idx]

            # fit subgroup model on resampled data
            mod.b    <- do.call(fit.subgroup, model$call)

            # calculate benefit scores and resulting
            # subgroup treatment effects on the original data
            benefit.scores.orig <- mod.b$predict(x)

            sbgrp.trt.eff.orig  <- subgroup.effects(benefit.scores.orig,
                                                    y, trt,
                                                    model$call$cutpoint)

            # subtract estimated bias for current bootstrap iteration
            boot.list[[1]][b,]  <- model$subgroup.trt.effects[[1]] -
                (mod.b$subgroup.trt.effects[[1]] - sbgrp.trt.eff.orig[[1]]) # bias estimate portion
            boot.list[[2]][b,,] <- model$subgroup.trt.effects[[2]] -
                (mod.b$subgroup.trt.effects[[2]] - sbgrp.trt.eff.orig[[2]]) # bias estimate portion
            boot.list[[3]][b,,] <- model$subgroup.trt.effects[[3]] -
                (mod.b$subgroup.trt.effects[[3]] - sbgrp.trt.eff.orig[[3]]) # bias estimate portion

        } else
        {   # bootstrap

            # bootstrap is not available because it
            # results in overly optimistic results
            samp.idx       <- sample.int(n.obs, n.obs, replace = TRUE)
            model$call$x   <- x[samp.idx,]
            model$call$y   <- y[samp.idx]
            model$call$trt <- trt[samp.idx]

            mod.b               <- do.call(fit.subgroup, model$call)
            boot.list[[1]][b,]  <- mod.b$subgroup.trt.effects[[1]] # subgroup-specific trt effects
            boot.list[[2]][b,,] <- mod.b$subgroup.trt.effects[[2]] # mean of outcome for 2x2 table (trt received vs recommended)
            boot.list[[3]][b,,] <- mod.b$subgroup.trt.effects[[3]] # sample sizes for 2x2 table

        }
    }

    # ugly way to handle cases where
    # some subgroups have no members
    for (l in 1:length(boot.list))
    {
        boot.list[[l]][is.nan(boot.list[[l]])] <- 0
    }

    # compute averages and standard
    # deviations across iterations
    summary.stats    <- list(colMeans(boot.list[[1]]),
                             apply(boot.list[[2]], c(2, 3), mean),
                             apply(boot.list[[3]], c(2, 3), mean))

    summary.stats.se <- list(apply(boot.list[[1]], 2, sd),
                             apply(boot.list[[2]], c(2, 3), sd),
                             apply(boot.list[[3]], c(2, 3), sd))

    names(summary.stats) <- names(boot.list) <- names(model$subgroup.trt.effects)

    names(summary.stats$subgroup.effects) <- names(model$subgroup.trt.effects$subgroup.effects)
    dimnames(summary.stats$sample.sizes)  <- dimnames(model$subgroup.trt.effects$sample.sizes)

    names(summary.stats.se[[1]])          <- names(model$subgroup.trt.effects$subgroup.effects)
    dimnames(summary.stats.se[[3]])       <- dimnames(model$subgroup.trt.effects$sample.sizes)

    names(summary.stats.se) <- paste("SE", names(summary.stats), sep = ".")


    ret <- list(avg.results  = summary.stats,    # means
                se.results   = summary.stats.se, # std errors
                boot.results = boot.list,        # this is a list of results for each iter
                family       = model$family,     # model family
                loss         = model$loss,       # model loss
                method       = model$method,     # subgroup method (weighting vs a-learning)
                val.method   = method)
    class(ret) <- "subgroup_validated"
    ret
}
