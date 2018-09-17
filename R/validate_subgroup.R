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
#' @param benefit.score.quantiles a vector of quantiles (between 0 and 1) of the benefit score values
#' for which to return bootstrapped information about the subgroups. ie if one of the quantile values is 0.5, the
#' median value of the benefit scores will be used as a cutoff to determine subgroups and summary statistics
#' will be returned about these subgroups
#' @param parallel Should the loop over replications be parallelized? If \code{FALSE}, then no, if \code{TRUE}, then yes.
#' If user sets \code{parallel = TRUE} and the fitted \code{fit.subgroup()} object uses the parallel version of
#' an internal model, say for \code{cv.glmnet()}, then the internal parallelization will be overridden so as
#' not to create a conflict of parallelism.
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models,
#' \code{\link[personalized]{plot.subgroup_validated}} for plotting of validation results, and
#' \code{\link[personalized]{print.subgroup_validated}} for arguments for printing options for \code{validate.subgroup()}.
#' @references Chen, S., Tian, L., Cai, T. and Yu, M. (2017), A general statistical framework for subgroup identification
#' and comparative treatment scoring. Biometrics. doi:10.1111/biom.12676
#'
#' Harrell, F. E., Lee, K. L., and Mark, D. B. (1996). Tutorial in biostatistics multivariable prognostic models: issues in developing models,
#' evaluating assumptions and adequacy, and measuring and reducing errors. Statistics in medicine, 15, 361-387.
#' doi:10.1002/(SICI)1097-0258(19960229)15:4<361::AID-SIM168>3.0.CO;2-4
#' @details Estimates of various quantities conditional on subgroups and treatment statuses are provided and displayed
#' via the \code{\link[personalized]{print.subgroup_validated}} function:
#' \enumerate{
#' \item{"Conditional expected outcomes"}{ The first results shown when printing
#' a \code{subgroup_validated} object are estimates of the expected outcomes conditional on
#' the estimated subgroups (i.e. which subgroup is 'recommended' by the model) and conditional
#' on treatment/intervention status. If there are two total treatment options, this results in a 2x2 table
#' of expected conditional outcomes. }
#' \item{"Treatment effects conditional on subgroups"}{ The second results shown when printing
#' a \code{subgroup_validated} object are estimates of the expected outcomes conditional on
#' the estimated subgroups. If the treatment takes levels \eqn{j \in \{1, \dots, K\}}, a total of \eqn{K}
#' conditional treatment effects will be shown. For example, of the outcome is continuous, the
#' \eqn{j}th conditional treatment effect is defined as \eqn{E(Y|Trt = j, Subgroup=j) - E(Y|Trt = j, Subgroup =/= j)}},
#' where \eqn{Subgroup=j} if treatment \eqn{j} is recommended, i.e. treatment \eqn{j} results in the largest/best
#' expected potential outcomes given the fitted model.
#' \item{"Overall treatment effect conditional on subgroups "}{ The third quantity displayed shows the overall improvement
#' in outcomes resulting from all treatment recommendations. This is essentially an average over all of the
#' conditional treatment effects weighted by the proportion of the population recommended each respective
#' treatment level.}
#' }
#' @return An object of class \code{"subgroup_validated"}
#' \item{avg.results}{Estimates of average conditional treatment effects when
#' subgroups are determined based on the provided cutoff value for the benefit score. For example,
#' if \code{cutoff = 0} and there is a treatment and control only, then the treatment is
#' recommended if the benefit score is greater than 0.}
#' \item{se.results}{Standard errors of the estimates from \code{avg.estimates}}
#' \item{boot.results}{Contains the individual results for each replication. \code{avg.results} is comprised
#' of averages of the values from \code{boot.results}}
#' \item{avg.quantile.results}{Estimates of average conditional treatment effects when
#' subgroups are determined based on different quntile cutoff values for the benefit score. For example,
#' if \code{benefit.score.quantiles = 0.75} and there is a treatment and control only, then the treatment is
#' recommended if the benefit score is greater than the 75th upper quantile of all benefit scores. If multiple quantile
#' values are provided, e.g. \code{benefit.score.quantiles = c(0.15, 0.5, 0.85)}, then results will be provided
#' for all quantile levels.}
#' \item{se.quantile.results}{Standard errors corresponding to \code{avg.quantile.results}}
#' \item{boot.results.quantiles}{Contains the individual results for each replication. \code{avg.quantile.results} is comprised
#' of averages of the values from \code{boot.results.quantiles}}
#' \item{family}{Family of the outcome. For example, \code{"gaussian"} for continuous outcomes}
#' \item{method}{Method used for subgroup identification model. Weighting or A-learning}
#' \item{n.trts}{The number of treatment levels}
#' \item{comparison.trts}{All treatment levels other than the reference level}
#' \item{reference.trt}{The reference level for the treatment. This should usually be the control group/level}
#' \item{larger.outcome.better}{If larger outcomes are preferred for this model}
#' \item{cutpoint}{Benefit score cutoff value used for determining subgroups}
#' \item{val.method}{Method used for validation}
#' \item{iterations}{Number of replications used in the validation process}
#' @importFrom stats predict sd
#' @import foreach
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
#' valmod <- validate.subgroup(subgrp.model, B = 3,
#'                             method = "training_test",
#'                             train.fraction = 0.75)
#' valmod
#'
#' print(valmod, which.quant = c(4, 5))
#'
#' @export
validate.subgroup <- function(model,
                              B              = 50L,
                              method         = c("training_test_replication",
                                                 "boot_bias_correction"),
                              train.fraction = 0.5,
                              benefit.score.quantiles = c(0.1666667, 0.3333333, 0.5000000, 0.6666667, 0.8333333),
                              parallel       = FALSE)
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
    x               <- model$call$x
    trt             <- model$call$trt
    y               <- model$call$y
    match.id        <- model$call$match.id
    propensity.func <- model$call$propensity.func

    if (is.null(model$pi.x)) stop("Please re-fit model and run again")

    pi.x            <- model$pi.x

    ## override any internal parallelization
    ## if there is a conflict of parallelization
    if (parallel & ("parallel" %in% names(model$call)) )
    {
        model$call$parallel <- FALSE
    }

    if (!all(benefit.score.quantiles < 1) || !all(benefit.score.quantiles > 0) )
    {
        stop("all quantile values must be strictly between 0 and 1")
    }

    if (is.null(benefit.score.quantiles))
    {
        benefit.score.quantiles <- seq(0, 1, length.out = 7)[-c(1,7)]
    }

    if (!is.vector(benefit.score.quantiles)) stop("benefit.score.quantiles must be a vector")

    n.quantiles <- length(benefit.score.quantiles)

    if (n.quantiles < 1)
    {
        benefit.score.quantiles <- seq(0, 1, length.out = 7)[-c(1,7)]
        n.quantiles <- length(benefit.score.quantiles)
    }

    n.obs <- NROW(x)

    # create objects to store results
    boot.list      <- vector(mode = "list", length = length(model$subgroup.trt.effects))
    boot.list[[1]] <- array(NA, dim = c(B, length(model$subgroup.trt.effects[[1]])))
    boot.list[[2]] <- boot.list[[3]] <- array(NA, dim = c(B, dim(model$subgroup.trt.effects[[2]])))
    # Add extra element to boot.list to collect coefficient information on each run
    boot.list[[4]] <- matrix(rep(NA,B*length(model$coefficients)),
                             ncol=B,
                             dimnames = list(row.names(model$coefficients),paste0("B",1:B)))
    boot.list[[5]] <- numeric(B)

    dimnames(boot.list[[2]]) <- dimnames(boot.list[[3]]) <-
        c(list(NULL), dimnames(model$subgroup.trt.effects$avg.outcomes))

    boot.list.quantiles <- rep(list(boot.list), n.quantiles)

    sbgrp.trt.eff.mod.orig <- vector(mode = "list", length = n.quantiles)
    for (q in 1:n.quantiles)
    {
        cutpt     <- round(100 * benefit.score.quantiles[q])
        cutpt.txt <- paste0("quant", cutpt)

        # estimate subgroup treatment effects on original data with original model
        sbgrp.trt.eff.mod.orig[[q]]  <- subgroup.effects(model$benefit.scores,
                                                         y,
                                                         trt,
                                                         model$pi.x,
                                                         cutpt.txt,
                                                         model$larger.outcome.better,
                                                         model$reference.trt)
    }

    if (parallel)
    {
        ## no longer needed:
        #comb <- function(x, ...)
        #{
        #    lapply(seq_along(x),
        #           function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
        #}

        outlist <- foreach(i = seq(B) ) %dopar% #, .combine = "comb", .multicombine = TRUE, .init=list(list(), list(), list(), list(), list()))
                           {
                               if (method == "training_test_replication")
                               {
                                   # randomly split/partition data into training and testing sets
                                   if (is.null(match.id))
                                   {
                                       train.samp.size <- floor(n.obs * train.fraction)
                                       samp.idx        <- sample.int(n.obs, train.samp.size, replace = FALSE)
                                   } else
                                   {
                                       # Draw at the cluster level (train.fraction interpreted as fraction of clusters)
                                       n.levels    <- length(levels(match.id))
                                       samp.levels <- sample.int(n.levels, n.levels * train.fraction, replace = FALSE)
                                       samp.idx    <- which(match.id %in% levels(match.id)[samp.levels])
                                       model$call$match.id <- match.id[samp.idx]
                                   }
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
                                   trt.test  <- trt[-samp.idx]
                                   pi.x.test <- pi.x[-samp.idx]

                                   # fit subgroup model on training data
                                   mod.b    <- do.call(fit.subgroup, model$call)

                                   # compute benefit scores on testing data
                                   benefit.scores.test <- mod.b$predict(x.test)

                                   # estimate subgroup treatment effects on test data
                                   sbgrp.trt.eff.test  <- subgroup.effects(benefit.scores.test,
                                                                           y.test, trt.test,
                                                                           pi.x.test,
                                                                           model$call$cutpoint,
                                                                           model$larger.outcome.better,
                                                                           model$reference.trt)

                                   # save results

                                   res  <- rep(list(list(sbgrp.trt.eff.test[[1]],
                                                    sbgrp.trt.eff.test[[2]],
                                                    sbgrp.trt.eff.test[[3]],
                                                    as.vector(mod.b$coefficients),
                                                    sbgrp.trt.eff.test[[4]])), 1 + n.quantiles)

                                   for (q in 1:n.quantiles)
                                   {

                                       cutpt     <- round(100 * benefit.score.quantiles[q])
                                       cutpt.txt <- paste0("quant", cutpt)

                                       # estimate subgroup treatment effects on test data
                                       sbgrp.trt.eff.test.q  <- subgroup.effects(benefit.scores.test,
                                                                                 y.test, trt.test,
                                                                                 pi.x.test,
                                                                                 cutpt.txt,
                                                                                 model$larger.outcome.better,
                                                                                 model$reference.trt)

                                       # save results
                                       res[[1 + q]][[1]] <- sbgrp.trt.eff.test.q[[1]]
                                       res[[1 + q]][[2]] <- sbgrp.trt.eff.test.q[[2]]
                                       res[[1 + q]][[3]] <- sbgrp.trt.eff.test.q[[3]]
                                       res[[1 + q]][[4]] <- as.vector(mod.b$coefficients)
                                       res[[1 + q]][[5]] <- sbgrp.trt.eff.test.q[[4]]

                                   }

                               } else if (method == "boot_bias_correction")
                               {   # bootstrap bias correction

                                   # take a bootstrap sample with replacement
                                   if (is.null(match.id))
                                   {
                                       samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
                                   } else
                                   {
                                       # Draw at the cluster level
                                       samp.levels <- sample(levels(match.id), replace = TRUE)
                                       samp.lookup <- lapply(samp.levels, function(z) {which(match.id == z)})
                                       samp.idx    <- unlist(samp.lookup)
                                       # Remap matching IDs so that each cluster draw is assigned a unique matching ID
                                       samp.lengths           <- lapply(samp.lookup,length)
                                       model$call$match.id <- unlist(lapply(1:length(samp.lengths),function(z) {rep(z,samp.lengths[[z]])}))
                                   }


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

                                   pi.x.b   <- mod.b$pi.x

                                   # calculate benefit scores and resulting
                                   # subgroup treatment effects on the original data
                                   benefit.scores.orig <- mod.b$predict(x)

                                   sbgrp.trt.eff.orig  <- subgroup.effects(benefit.scores.orig,
                                                                           y, trt, model$pi.x,
                                                                           model$call$cutpoint,
                                                                           model$larger.outcome.better,
                                                                           model$reference.trt)

                                   # subtract estimated bias for current bootstrap iteration
                                   res  <- rep(list(list(model$subgroup.trt.effects[[1]] -
                                                             (mod.b$subgroup.trt.effects[[1]] - sbgrp.trt.eff.orig[[1]]), # bias estimate portion
                                                         model$subgroup.trt.effects[[2]] -
                                                             (mod.b$subgroup.trt.effects[[2]] - sbgrp.trt.eff.orig[[2]]), # bias estimate portion
                                                         model$subgroup.trt.effects[[3]] -
                                                             (mod.b$subgroup.trt.effects[[3]] - sbgrp.trt.eff.orig[[3]]), # bias estimate portion
                                                         as.vector(mod.b$coefficients),
                                                         model$subgroup.trt.effects[[4]] -
                                                             (mod.b$subgroup.trt.effects[[4]] - sbgrp.trt.eff.orig[[4]]))), 1 + n.quantiles)
                                   # bias estimate portion


                                   for (q in 1:n.quantiles)
                                   {

                                       cutpt     <- round(100 * benefit.score.quantiles[q])
                                       cutpt.txt <- paste0("quant", cutpt)

                                       # estimate subgroup treatment effects on original data with resampled model
                                       sbgrp.trt.eff.orig.q  <- subgroup.effects(benefit.scores.orig,
                                                                                 y,
                                                                                 trt,
                                                                                 model$pi.x,
                                                                                 cutpt.txt,
                                                                                 model$larger.outcome.better,
                                                                                 model$reference.trt)

                                       # estimate subgroup treatment effects on resampled data with resampled model
                                       sbgrp.trt.eff.mod  <- subgroup.effects(mod.b$benefit.scores,
                                                                              y[samp.idx],
                                                                              trt[samp.idx],
                                                                              pi.x.b,
                                                                              cutpt.txt,
                                                                              model$larger.outcome.better,
                                                                              model$reference.trt)


                                       # save results
                                       # subtract estimated bias for current bootstrap iteration
                                       res[[1 + q]][[1]]  <- sbgrp.trt.eff.mod.orig[[q]][[1]] -
                                           (sbgrp.trt.eff.mod[[1]] - sbgrp.trt.eff.orig.q[[1]]) # bias estimate portion
                                       res[[1 + q]][[2]]  <- sbgrp.trt.eff.mod.orig[[q]][[2]] -
                                           (sbgrp.trt.eff.mod[[2]] - sbgrp.trt.eff.orig.q[[2]]) # bias estimate portion
                                       res[[1 + q]][[3]]  <- sbgrp.trt.eff.mod.orig[[q]][[3]] -
                                           (sbgrp.trt.eff.mod[[3]] - sbgrp.trt.eff.orig.q[[3]]) # bias estimate portion
                                       res[[1 + q]][[4]]  <- as.vector(mod.b$coefficients)
                                       res[[1 + q]][[5]]  <- sbgrp.trt.eff.mod.orig[[q]][[4]] -
                                           (sbgrp.trt.eff.mod[[4]] - sbgrp.trt.eff.orig.q[[4]]) # bias estimate portion
                                   }

                               }

                               res
                           } # end parallel foreach loop

        ## collect results as they are collected for non-parallel loop
        #boot.list[[1]] <- unname(Reduce('rbind', outlist[[1]]))
        for (b in 1:B)
        {
            boot.list[[1]][b,]  <- outlist[[b]][[1]][[1]]
            boot.list[[2]][b,,] <- outlist[[b]][[1]][[2]]
            boot.list[[3]][b,,] <- outlist[[b]][[1]][[3]]
            boot.list[[4]][,b]  <- outlist[[b]][[1]][[4]]
            boot.list[[5]][b]   <- outlist[[b]][[1]][[5]]

            for (q in 1:n.quantiles)
            {
                boot.list.quantiles[[q]][[1]][b,]  <- outlist[[b]][[1 + q]][[1]]
                boot.list.quantiles[[q]][[2]][b,,] <- outlist[[b]][[1 + q]][[2]]
                boot.list.quantiles[[q]][[3]][b,,] <- outlist[[b]][[1 + q]][[3]]
                boot.list.quantiles[[q]][[4]][,b]  <- outlist[[b]][[1 + q]][[4]]
                boot.list.quantiles[[q]][[5]][b]   <- outlist[[b]][[1 + q]][[5]]
            }
        }

        #boot.list[[4]] <- unname(Reduce('cbind', outlist[[4]]))
        #boot.list[[5]] <- unname(Reduce('c', outlist[[5]]))
        rm(outlist)
    } else
    {
        for (b in 1:B)
        {
            if (method == "training_test_replication")
            {
                # randomly split/partition data into training and testing sets
                if (is.null(match.id))
                {
                    train.samp.size <- floor(n.obs * train.fraction)
                    samp.idx        <- sample.int(n.obs, train.samp.size, replace = FALSE)
                } else
                {
                    # Draw at the cluster level (train.fraction interpreted as fraction of clusters)
                    n.levels    <- length(levels(match.id))
                    samp.levels <- sample.int(n.levels, n.levels * train.fraction, replace = FALSE)
                    samp.idx    <- which(match.id %in% levels(match.id)[samp.levels])
                    model$call$match.id <- match.id[samp.idx]
                }
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
                trt.test  <- trt[-samp.idx]
                pi.x.test <- pi.x[-samp.idx]

                # fit subgroup model on training data
                mod.b    <- do.call(fit.subgroup, model$call)

                # compute benefit scores on testing data
                benefit.scores.test <- mod.b$predict(x.test)

                # estimate subgroup treatment effects on test data
                sbgrp.trt.eff.test  <- subgroup.effects(benefit.scores.test,
                                                        y.test, trt.test,
                                                        pi.x.test,
                                                        model$call$cutpoint,
                                                        model$larger.outcome.better,
                                                        model$reference.trt)

                # save results
                boot.list[[1]][b,]  <- sbgrp.trt.eff.test[[1]]
                boot.list[[2]][b,,] <- sbgrp.trt.eff.test[[2]]
                boot.list[[3]][b,,] <- sbgrp.trt.eff.test[[3]]
                boot.list[[4]][,b]  <- as.vector(mod.b$coefficients)
                boot.list[[5]][b]   <- sbgrp.trt.eff.test[[4]]

                for (q in 1:n.quantiles)
                {

                    cutpt <- round(100 * benefit.score.quantiles[q])
                    cutpt.txt <- paste0("quant", cutpt)

                    # estimate subgroup treatment effects on test data
                    sbgrp.trt.eff.test.q  <- subgroup.effects(benefit.scores.test,
                                                              y.test, trt.test,
                                                              pi.x.test,
                                                              cutpt.txt,
                                                              model$larger.outcome.better,
                                                              model$reference.trt)

                    # save results
                    boot.list.quantiles[[q]][[1]][b,]  <- sbgrp.trt.eff.test.q[[1]]
                    boot.list.quantiles[[q]][[2]][b,,] <- sbgrp.trt.eff.test.q[[2]]
                    boot.list.quantiles[[q]][[3]][b,,] <- sbgrp.trt.eff.test.q[[3]]
                    boot.list.quantiles[[q]][[4]][,b]  <- as.vector(mod.b$coefficients)
                    boot.list.quantiles[[q]][[5]][b]   <- sbgrp.trt.eff.test.q[[4]]

                }

            } else if (method == "boot_bias_correction")
            {   # bootstrap bias correction

                # take a bootstrap sample with replacement
                if (is.null(match.id))
                {
                    samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
                } else
                {
                    # Draw at the cluster level
                    samp.levels <- sample(levels(match.id), replace = TRUE)
                    samp.lookup <- lapply(samp.levels, function(z) {which(match.id == z)})
                    samp.idx    <- unlist(samp.lookup)
                    # Remap matching IDs so that each cluster draw is assigned a unique matching ID
                    samp.lengths           <- lapply(samp.lookup,length)
                    model$call$match.id <- unlist(lapply(1:length(samp.lengths),function(z){rep(z,samp.lengths[[z]])}))
                }
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
                pi.x.b   <- mod.b$pi.x

                # calculate benefit scores and resulting
                # subgroup treatment effects on the original data
                benefit.scores.orig <- mod.b$predict(x)

                sbgrp.trt.eff.orig  <- subgroup.effects(benefit.scores.orig,
                                                        y, trt, model$pi.x,
                                                        model$call$cutpoint,
                                                        model$larger.outcome.better,
                                                        model$reference.trt)

                # subtract estimated bias for current bootstrap iteration
                boot.list[[1]][b,]  <- model$subgroup.trt.effects[[1]] -
                    (mod.b$subgroup.trt.effects[[1]] - sbgrp.trt.eff.orig[[1]]) # bias estimate portion
                boot.list[[2]][b,,] <- model$subgroup.trt.effects[[2]] -
                    (mod.b$subgroup.trt.effects[[2]] - sbgrp.trt.eff.orig[[2]]) # bias estimate portion
                boot.list[[3]][b,,] <- model$subgroup.trt.effects[[3]] -
                    (mod.b$subgroup.trt.effects[[3]] - sbgrp.trt.eff.orig[[3]]) # bias estimate portion
                boot.list[[4]][,b]  <- as.vector(mod.b$coefficients)
                boot.list[[5]][b]   <- model$subgroup.trt.effects[[4]] -
                    (mod.b$subgroup.trt.effects[[4]] - sbgrp.trt.eff.orig[[4]]) # bias estimate portion


                for (q in 1:n.quantiles)
                {

                    cutpt <- round(100 * benefit.score.quantiles[q])
                    cutpt.txt <- paste0("quant", cutpt)

                    # estimate subgroup treatment effects on original data with resampled model
                    sbgrp.trt.eff.orig  <- subgroup.effects(benefit.scores.orig,
                                                            y, trt, model$pi.x,
                                                            cutpt.txt,
                                                            model$larger.outcome.better,
                                                            model$reference.trt)

                    # estimate subgroup treatment effects on resampled data with resampled model
                    sbgrp.trt.eff.mod  <- subgroup.effects(mod.b$benefit.scores,
                                                           y[samp.idx],
                                                           trt[samp.idx],
                                                           pi.x.b,
                                                           cutpt.txt,
                                                           model$larger.outcome.better,
                                                           model$reference.trt)

                    # save results
                    # subtract estimated bias for current bootstrap iteration
                    boot.list.quantiles[[q]][[1]][b,]  <- sbgrp.trt.eff.mod.orig[[q]][[1]] -
                        (sbgrp.trt.eff.mod[[1]] - sbgrp.trt.eff.orig[[1]]) # bias estimate portion
                    boot.list.quantiles[[q]][[2]][b,,] <- sbgrp.trt.eff.mod.orig[[q]][[2]] -
                        (sbgrp.trt.eff.mod[[2]] - sbgrp.trt.eff.orig[[2]]) # bias estimate portion
                    boot.list.quantiles[[q]][[3]][b,,] <- sbgrp.trt.eff.mod.orig[[q]][[3]] -
                        (sbgrp.trt.eff.mod[[3]] - sbgrp.trt.eff.orig[[3]]) # bias estimate portion
                    boot.list.quantiles[[q]][[4]][,b]  <- as.vector(mod.b$coefficients)
                    boot.list.quantiles[[q]][[5]][b]   <- sbgrp.trt.eff.mod.orig[[q]][[4]] -
                        (sbgrp.trt.eff.mod[[4]] - sbgrp.trt.eff.orig[[4]]) # bias estimate portion

                }

            }
        } ## end resampling loop
    }

    # ugly way to handle cases where
    # some subgroups have no members
    #for (l in 1:length(boot.list))
    #{
    #    boot.list[[l]][is.nan(boot.list[[l]])] <- 0
    #}

    # compute averages and standard
    # deviations across iterations
    summary.stats    <- list(colMeans(boot.list[[1]], na.rm = TRUE),
                             apply(boot.list[[2]], c(2, 3), function(x) mean(x, na.rm = TRUE)),
                             apply(boot.list[[3]], c(2, 3), function(x) mean(x, na.rm = TRUE)),
                             mean(boot.list[[5]], na.rm = TRUE))

    summary.stats.se <- list(apply(boot.list[[1]], 2, function(x) sd(x, na.rm = TRUE)),
                             apply(boot.list[[2]], c(2, 3), function(x) sd(x, na.rm = TRUE)),
                             apply(boot.list[[3]], c(2, 3), function(x) sd(x, na.rm = TRUE)),
                             sd(boot.list[[5]], na.rm = TRUE))

    names(summary.stats) <- names(model$subgroup.trt.effects)
    names(boot.list) <- 1:length(names(boot.list))
    names(boot.list)[c(1:3, 5)] <- names(summary.stats)
    names(boot.list)[4] <- "coefficients"

    colnames(boot.list[[1]]) <- names(model$subgroup.trt.effects$subgroup.effects)

    summary.stats.quantile <- summary.stats.quantile.se <- vector(mode = "list", length = n.quantiles)
    for (q in 1:n.quantiles)
    {
        summary.stats.quantile[[q]]    <- list(colMeans(boot.list.quantiles[[q]][[1]], na.rm = TRUE),
                                               apply(boot.list.quantiles[[q]][[2]], c(2, 3), function(x) mean(x, na.rm = TRUE)),
                                               apply(boot.list.quantiles[[q]][[3]], c(2, 3), function(x) mean(x, na.rm = TRUE)),
                                               mean(boot.list.quantiles[[q]][[5]], na.rm = TRUE))

        summary.stats.quantile.se[[q]] <- list(apply(boot.list.quantiles[[q]][[1]], 2, function(x) sd(x, na.rm = TRUE)),
                                               apply(boot.list.quantiles[[q]][[2]], c(2, 3), function(x) sd(x, na.rm = TRUE)),
                                               apply(boot.list.quantiles[[q]][[3]], c(2, 3), function(x) sd(x, na.rm = TRUE)),
                                               sd(boot.list.quantiles[[q]][[5]], na.rm = TRUE))

        names(summary.stats.quantile[[q]]) <- names(model$subgroup.trt.effects)

        names(summary.stats.quantile[[q]]$subgroup.effects) <- names(model$subgroup.trt.effects$subgroup.effects)
        dimnames(summary.stats.quantile[[q]]$sample.sizes)  <- dimnames(model$subgroup.trt.effects$sample.sizes)

        names(summary.stats.quantile.se[[q]][[1]])          <- names(model$subgroup.trt.effects$subgroup.effects)
        dimnames(summary.stats.quantile.se[[q]][[3]])       <- dimnames(model$subgroup.trt.effects$sample.sizes)

        names(summary.stats.quantile.se[[q]]) <- paste("SE", names(summary.stats), sep = ".")

        names(boot.list.quantiles[[q]]) <- 1:length(names(boot.list))
        names(boot.list.quantiles[[q]])[c(1:3, 5)] <- names(summary.stats)
        names(boot.list.quantiles[[q]])[4] <- "coefficients"

        colnames(boot.list.quantiles[[q]][[1]]) <- names(model$subgroup.trt.effects$subgroup.effects)
    }

    names(summary.stats.quantile) <- names(summary.stats.quantile.se) <-
        names(boot.list.quantiles) <- paste0("Quant_", round(100 * benefit.score.quantiles))


    names(summary.stats$subgroup.effects) <- names(model$subgroup.trt.effects$subgroup.effects)
    dimnames(summary.stats$sample.sizes)  <- dimnames(model$subgroup.trt.effects$sample.sizes)

    names(summary.stats.se[[1]])          <- names(model$subgroup.trt.effects$subgroup.effects)
    dimnames(summary.stats.se[[3]])       <- dimnames(model$subgroup.trt.effects$sample.sizes)

    names(summary.stats.se) <- paste("SE", names(summary.stats), sep = ".")

    names(summary.stats.quantile) <- names(summary.stats.quantile.se) <- paste0("Quantile_", round(100 * benefit.score.quantiles))

    ret <- list(avg.results  = summary.stats,    # means
                se.results   = summary.stats.se, # std errors
                boot.results = boot.list,        # this is a list of results for each iter
                avg.quantile.results = summary.stats.quantile,
                se.quantile.results  = summary.stats.quantile.se,
                boot.results.quantiles = boot.list.quantiles,
                family       = model$family,     # model family
                loss         = model$loss,       # model loss
                method       = model$method,     # subgroup method (weighting vs a-learning)
                n.trts       = model$n.trts,
                comparison.trts       = model$comparison.trts,
                reference.trt         = model$reference.trt,
                larger.outcome.better = model$larger.outcome.better,
                cutpoint              = model$cutpoint,
                val.method   = method,
                iterations   = B)
    class(ret) <- "subgroup_validated"
    ret
}
