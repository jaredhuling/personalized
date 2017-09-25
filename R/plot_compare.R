

#' Plot a comparison results for fitted or validated subgroup identification models
#'
#' @description Plots comparison of results for estimated subgroup treatment effects
#' @param ... the fitted (model or validation) objects to be plotted. Must be either
#' objects returned from \code{fit.subgroup()} or \code{validate.subgroup()}
#' @param type type of plot. \code{"density"} results in a density plot for the results
#' across all observations (if \code{x} is from \code{fit.subgroup()}) or if \code{x} is from \code{validate.subgroup()}
#' across iterations of either the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"boxplot"} results in boxplots across all observations/iterations of either
#' the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"interaction"} creates an
#' interaction plot for the different subgroups (crossing lines here means a meaningful subgroup)
#' @param avg.line boolean value of whether or not to plot a line for the average
#' value in addition to the density (only valid for \code{type = "density"})
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models and
#' \code{\link[personalized]{validate.subgroup}} for function which creates validation results.
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
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "sq_loss_lasso",
#'                            nfolds = 5)              # option for cv.glmnet
#'
#'
#' subgrp.modelg <- fit.subgroup(x = x, y = y,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "sq_loss_lasso_gam")
#'
#' plotCompare(subgrp.model, subgrp.modelg)
#'
#' @export
plotCompare <- function(...,
                        type = c("boxplot", "density", "interaction"),
                        avg.line = TRUE)
{

    list.obj  <- list(...)

    n.obj <- length(list.obj)

    if (n.obj == 0)
        stop("no fitted or validated model objects specified to be plotted")

    ok <- sapply(list.obj, function(lo) (class(lo)[1] == "subgroup_fitted") |
                                        (class(lo)[1] == "subgroup_validated"))

    ## get names of supplied objects
    ## and print them if they are not appropriate objects
    obj.names <- sapply(match.call(expand.dots = TRUE)[-1L], deparse)
    if (!all(ok))
    {
        n <- sum(!ok)
        stop(sprintf(ngettext(n, "object %s not either subgroup_fitted or subgroup_validated object",
                              "objects %s not either subgroup_fitted or subgroup_validated objects"), paste(sQuote(obj.names[!ok]),
                                                             collapse = ", ")), domain = NA)
    }


    type <- match.arg(type)

    avg.line <- as.logical(avg.line[1])

    dat.list <- avg.list <- vector(mode = "list", length = n.obj)

    types.vec <- character(n.obj)


    for (l in 1:n.obj)
    {

        obj.type <- class(list.obj[[l]])[1]

        if (obj.type == "subgroup_fitted")
        {
            avg.res  <- list.obj[[l]]$subgroup.trt.effects
        } else
        {
            avg.res  <- list.obj[[l]]$avg.results
        }

        types.vec[l] <- obj.type

        avg.res.2.plot <- data.frame(Recommended = rep(colnames(avg.res$avg.outcomes),
                                                       each = ncol(avg.res$avg.outcomes)),
                                     Received    = rep(rownames(avg.res$avg.outcomes),
                                                       ncol(avg.res$avg.outcomes)),
                                     Value       = as.vector(avg.res$avg.outcomes),
                                     Model       = obj.names[l])

        avg.list[[l]] <- avg.res.2.plot

        if (type != "interaction")
        {

            if (obj.type == "subgroup_fitted")
            {
                if (is.null(list.obj[[l]]$call))
                    stop("retcall argument must be set to TRUE for fitted model object to be plotted with
                     non-interaction plot")

                benefit.scores <- list.obj[[l]]$benefit.scores
                trt.rec        <- list.obj[[l]]$recommended.trts
                B <- NROW(benefit.scores)

                res.2.plot <- array(NA, dim = c(B, 3))
                colnames(res.2.plot) <- c("Recommended", "Received", "Value")
                res.2.plot <- data.frame(res.2.plot)

                cutpoint <- list.obj[[l]]$call$cutpoint
                lb       <- list.obj[[l]]$call$larger.outcome.better

                #res.2.plot[, 1] <- ifelse(trt.rec == 1, "Recommended Trt", "Recommended Ctrl")
                #res.2.plot[, 2] <- ifelse(x$call$trt == 1, "Received Trt", "Received Ctrl")
                res.2.plot[, 1] <- paste("Recommended", trt.rec)
                res.2.plot[, 2] <- paste("Received", list.obj[[l]]$call$trt)

                if (class(list.obj[[l]]$call$y) == "Surv")
                {
                    res.2.plot[, 3] <- log(list.obj[[l]]$call$y[,1])
                } else
                {
                    res.2.plot[, 3] <- list.obj[[l]]$call$y
                }
            } else
            {
                boot.res <- list.obj[[l]]$boot.results$avg.outcomes
                boot.dims <- dim(boot.res)


                n.entries <- prod(boot.dims[2:3])
                B <- boot.dims[1]

                res.2.plot <- array(NA, dim = c(B * n.entries, 3))
                colnames(res.2.plot) <- c("Recommended", "Received", "Value")
                res.2.plot <- data.frame(res.2.plot)

                for (b in 1:B)
                {
                    cur.idx <- c(((b - 1) * n.entries + 1):(b * n.entries))
                    res.2.plot[cur.idx, 1] <- rep(colnames(boot.res[b,,]),
                                                  each = ncol(boot.res[b,,]))
                    res.2.plot[cur.idx, 2] <- rep(rownames(boot.res[b,,]),
                                                  ncol(boot.res[b,,]))
                    res.2.plot[cur.idx, 3] <- as.vector(boot.res[b,,])
                }

            }

            res.2.plot$Model <- obj.names[l]
            dat.list[[l]] <- res.2.plot
        }

    }

    avg.res.2.plot <- Reduce(rbind, avg.list)
    res.2.plot     <- Reduce(rbind, dat.list)

    Recommended <- Received <- Value <- Model <- NULL


    if (type == "density")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Value, fill = Received)) +
            geom_density(alpha = 0.65) +
            geom_rug(aes(colour = Received), alpha = 0.85) +
            coord_flip() +
            facet_grid(Recommended ~ Model) +
            theme(legend.position = "bottom") +
            xlab("Outcome")
        if (avg.line)
        {
            pl.obj <- pl.obj + geom_vline(data = avg.res.2.plot,
                                          aes(xintercept = Value),
                                          size = 1.25) +
                geom_vline(data = avg.res.2.plot,
                           aes(xintercept = Value, colour = Received))
        }
    } else if (type == "boxplot")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Received, y = Value)) +
            geom_boxplot(aes(fill = Received)) +
            geom_rug(aes(colour = Received), alpha = 0.85) +
            facet_grid(Recommended ~ Model) +
            theme(legend.position = "bottom") +
            ylab("Outcome")
    } else
    {
        pl.obj <- ggplot(avg.res.2.plot,
                         aes(x = Recommended, y = Value, group = Received)) +
            geom_line(aes(colour = Received), size = 1.25) +
            geom_point(aes(colour = Received), size = 2) +
            facet_grid( ~ Model) +
            theme(legend.position = "bottom") +
            scale_x_discrete(expand = c(0.25, 0.25)) +
            ylab("Average Outcome")
    }
    pl.obj
}
