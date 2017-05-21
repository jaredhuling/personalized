#' Validating fitted subgroup identification models
#'
#' @description Validates subgroup treatment effects for fitted
#'  subgroup identification model class of Chen, et al (2017)
#'
#' @param model fitted model object returned by \code{fit.subgrp()} function
#' @param method validation method
#' @param B number of bootstrap replications or refitting replications
#' @param train.fraction fraction (between 0 and 1) of samples to be used for training in
#' training/test replication. Only used for \code{method = "training_test_replication"}
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
#' subgrp.model <- fit.subgrp(x = x, y = y,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            family = "gaussian",
#'                            loss   = "sq_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#' subgrp.model$subgroup.trt.effects
#'
#' x.test <- matrix(rnorm(10 * n.obs * n.vars, sd = 3), 10 * n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat.test   <- 0.5 + 0.5 * x.test[,21] - 0.5 * x.test[,41]
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
#' valmod <- validate.subgrp(subgrp.model, B = 10,
#'                           method = "training_test",
#'                           train.fraction = 0.75)
#' valmod$avg.results
#'
#' bene.score.test <- subgrp.model$predict(x.test)
#'
#' mean(y.test[bene.score.test > 0 & trt01.test == 1]) - mean(y.test[bene.score.test > 0 & trt01.test == 0])
#' mean(y.test[bene.score.test <= 0 & trt01.test == 0]) - mean(y.test[bene.score.test <= 0 & trt01.test == 1])
#'
#' quantile(valmod$boot.results[[1]][,1], c(0.025, 0.975))
#' quantile(valmod$boot.results[[1]][,2], c(0.025, 0.975))
#' @export
validate.subgrp <- function(model,
                            B              = 50L,
                            method         = c("training_test_replication",
                                               "boot_bias_correction"),
                            train.fraction = 0.5)
{
    method <- match.arg(method)


    if (class(model)[1] != "subgroup_fitted")
    {
        stop("model should be a fitted object returned by the 'fit.subgrp' function")
    }


    B <- as.integer(B[1])
    if (B <= 1) stop("B must be a strictly positive integer")

    train.fraction <- as.numeric(train.fraction[1])

    if (train.fraction >= 1 | train.fraction <= 0)
    {
        stop("train.fraction must be between 0 and 1")
    }

    if (is.null(model$call)) stop("retcall argument must be set to TRUE for fitted model")


    model$call$retcall <- FALSE

    x   <- model$call$x
    trt <- model$call$trt
    y   <- model$call$y

    n.obs <- NROW(x)

    boot.list <- vector(mode = "list", length = length(model$subgroup.trt.effects))
    boot.list[[1]] <- array(NA, dim = c(B, length(model$subgroup.trt.effects[[1]])))
    boot.list[[2]] <- boot.list[[3]] <- array(NA, dim = c(B, dim(model$subgroup.trt.effects[[2]])))

    dimnames(boot.list[[2]]) <- dimnames(boot.list[[2]]) <-
        c(list(NULL), dimnames(model$subgroup.trt.effects$avg.outcomes))

    for (b in 1:B)
    {
        if (method == "training_test_replication")
        {
            train.samp.size <- floor(n.obs * train.fraction)
            samp.idx        <- sample.int(n.obs, train.samp.size, replace = FALSE)
            model$call$x    <- x[samp.idx,]
            model$call$y    <- y[samp.idx]
            model$call$trt  <- trt[samp.idx]

            x.test   <- x[-samp.idx,]
            y.test   <- y[-samp.idx]
            trt.test <- trt[-samp.idx]

            mod.b    <- do.call(fit.subgrp, model$call)

            benefit.scores.test <- mod.b$predict(x.test)

            sbgrp.trt.eff.test  <- subgrp.benefit(benefit.scores.test,
                                                  y.test, trt.test,
                                                  model$call$cutpoint)

            boot.list[[1]][b,]  <- sbgrp.trt.eff.test[[1]]
            boot.list[[2]][b,,] <- sbgrp.trt.eff.test[[2]]
            boot.list[[3]][b,,] <- sbgrp.trt.eff.test[[3]]
        } else if (method == "boot_bias_correction")
        {   # bootstrap bias correction
            samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
            model$call$x   <- x[samp.idx,]
            model$call$y   <- y[samp.idx]
            model$call$trt <- trt[samp.idx]

            mod.b    <- do.call(fit.subgrp, model$call)

            benefit.scores.orig <- mod.b$predict(x)

            sbgrp.trt.eff.orig  <- subgrp.benefit(benefit.scores.orig,
                                                  y, trt,
                                                  model$call$cutpoint)

            ## subtract estimated bias for current bootstrap iteration
            boot.list[[1]][b,]  <- model$subgroup.trt.effects[[1]] -
                (mod.b$subgroup.trt.effects[[1]] - sbgrp.trt.eff.orig[[1]]) # bias portion
            boot.list[[2]][b,,] <- model$subgroup.trt.effects[[2]] -
                (mod.b$subgroup.trt.effects[[2]] - sbgrp.trt.eff.orig[[2]]) # bias portion
            boot.list[[3]][b,,] <- model$subgroup.trt.effects[[3]] -
                (mod.b$subgroup.trt.effects[[3]] - sbgrp.trt.eff.orig[[3]]) # bias portion
        } else
        {   # bootstrap
            samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
            model$call$x   <- x[samp.idx,]
            model$call$y   <- y[samp.idx]
            model$call$trt <- trt[samp.idx]

            mod.b    <- do.call(fit.subgrp, model$call)
            boot.list[[1]][b,]  <- mod.b$subgroup.trt.effects[[1]]
            boot.list[[2]][b,,] <- mod.b$subgroup.trt.effects[[2]]
            boot.list[[3]][b,,] <- mod.b$subgroup.trt.effects[[3]]
        }
    }

    for (l in 1:length(boot.list))
    {
        boot.list[[l]][is.nan(boot.list[[l]])] <- 0
    }

    summary.stats <- list(colMeans(boot.list[[1]]),
                          apply(boot.list[[2]], c(2, 3), mean),
                          apply(boot.list[[3]], c(2, 3), mean))

    summary.stats.se <- list(apply(boot.list[[1]], 2, sd),
                             apply(boot.list[[2]], c(2, 3), sd),
                             apply(boot.list[[3]], c(2, 3), sd))

    names(summary.stats) <- names(boot.list) <- names(model$subgroup.trt.effects)

    ret <- list(avg.results  = summary.stats,
                se.results   = summary.stats.se,
                boot.results = boot.list)
    class(ret) <- c("subgroup_validated", class(ret))
    ret
}
