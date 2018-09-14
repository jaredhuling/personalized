

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
#' interaction plot for the different subgroups (crossing lines here means a meaningful subgroup).
#' \code{"conditional"} plots smoothed (via a GAM smoother) means of the outcomes as a function of the estimated benefit score
#' separately for the treated and untreated groups.
#' @param avg.line boolean value of whether or not to plot a line for the average
#' value in addition to the density (only valid for \code{type = "density"})
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models and
#' \code{\link[personalized]{validate.subgroup}} for function which creates validation results.
#'
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 100
#' n.vars <- 15
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,1] - 0.5 * x[,4]
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
#' subgrp.model.o <- fit.subgroup(x = x, y = y,
#'                            trt = trt01,
#'                            propensity.func = prop.func,
#'                            loss   = "owl_logistic_flip_loss_lasso",
#'                            nfolds = 5)
#'
#' plotCompare(subgrp.model, subgrp.model.o)
#'
#' @export
plotCompare <- function(...,
                        type = c("boxplot", "density", "interaction", "conditional"),
                        avg.line = TRUE)
{

    list.obj  <- list(...)

    n.obj <- length(list.obj)

    bs <- Outcome <- NULL

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

    is_fitted_obj    <- sapply(list.obj, function(lo) (class(lo)[1] == "subgroup_fitted"))
    is_validated_obj <- sapply(list.obj, function(lo) (class(lo)[1] == "subgroup_validated"))

    if (type == "conditional" & (!all(is_fitted_obj) & !all(is_validated_obj)))
    {
        stop("type == 'conditional' only allowed if all objects are subgroup_fitted or subgroup_validated,
             not a mix of the two")
    }

    all_binom <- all(sapply(list.obj, function(lo) lo$family == "binomial"))


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
                                     Received    = gsub("^Received ", "", rep(rownames(avg.res$avg.outcomes),
                                                       ncol(avg.res$avg.outcomes))),
                                     Value       = as.vector(avg.res$avg.outcomes),
                                     Model       = obj.names[l])

        avg.res.2.plot$Received <- as.factor(avg.res.2.plot$Received)

        avg.res.2.plot.dens <- avg.res.2.plot

        avg.res.2.plot$Recommended <- gsub("^Recommended ", "", avg.res.2.plot$Recommended)

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

                if (type != "conditional")
                {
                    res.2.plot <- array(NA, dim = c(B, 3))
                    colnames(res.2.plot) <- c("Recommended", "Received", "Value")
                    res.2.plot <- data.frame(res.2.plot)

                    cutpoint <- list.obj[[l]]$call$cutpoint
                    lb       <- list.obj[[l]]$call$larger.outcome.better

                    #res.2.plot[, 1] <- ifelse(trt.rec == 1, "Recommended Trt", "Recommended Ctrl")
                    #res.2.plot[, 2] <- ifelse(x$call$trt == 1, "Received Trt", "Received Ctrl")
                    res.2.plot[, 1] <- paste("Recommended", trt.rec)
                    res.2.plot[, 2] <- list.obj[[l]]$call$trt #paste("Received", list.obj[[l]]$call$trt)

                    if (class(list.obj[[l]]$call$y) == "Surv")
                    {
                        res.2.plot[, 3] <- log(list.obj[[l]]$call$y[,1])
                    } else
                    {
                        res.2.plot[, 3] <- list.obj[[l]]$call$y
                    }
                } else
                {
                    res.2.plot <- array(NA, dim = c(B, 3))
                    colnames(res.2.plot) <- c("bs", "Received", "Outcome")
                    res.2.plot <- data.frame(res.2.plot)

                    cutpoint <- list.obj[[l]]$call$cutpoint
                    lb       <- list.obj[[l]]$call$larger.outcome.better

                    #res.2.plot[, 1] <- ifelse(trt.rec == 1, "Recommended Trt", "Recommended Ctrl")
                    #res.2.plot[, 2] <- ifelse(x$call$trt == 1, "Received Trt", "Received Ctrl")
                    res.2.plot[, 1] <- benefit.scores
                    res.2.plot[, 2] <- list.obj[[l]]$call$trt #paste("Received", list.obj[[l]]$call$trt)

                    if (class(list.obj[[l]]$call$y) == "Surv")
                    {
                        res.2.plot[, 3] <- log(list.obj[[l]]$call$y[,1])
                    } else
                    {
                        res.2.plot[, 3] <- list.obj[[l]]$call$y
                    }
                }
            } else
            {
                boot.res <- list.obj[[l]]$boot.results$avg.outcomes
                boot.dims <- dim(boot.res)


                n.entries <- prod(boot.dims[2:3])
                B <- boot.dims[1]

                if (type != "conditional")
                {
                    res.2.plot <- array(NA, dim = c(B * n.entries, 3))
                    colnames(res.2.plot) <- c("Recommended", "Received", "Value")
                    res.2.plot <- data.frame(res.2.plot)

                    for (b in 1:B)
                    {
                        cur.idx <- c(((b - 1) * n.entries + 1):(b * n.entries))
                        res.2.plot[cur.idx, 1] <- rep(colnames(boot.res[b,,]),
                                                      each = ncol(boot.res[b,,]))
                        res.2.plot[cur.idx, 2] <- gsub("^Received ", "", rep(rownames(boot.res[b,,]),
                                                      ncol(boot.res[b,,])))
                        res.2.plot[cur.idx, 3] <- as.vector(boot.res[b,,])
                    }
                } else
                {
                    n.quantiles    <- length(list.obj[[l]]$boot.results.quantiles)
                    quantile.names <- paste("Cutoff:", names(list.obj[[l]]$boot.results.quantiles))

                    res.2.plot <- array(NA, dim = c(B * n.entries * n.quantiles, 4))
                    colnames(res.2.plot) <- c("Recommended", "Received", "Value", "Quantile")
                    res.2.plot <- data.frame(res.2.plot)

                    ct <- 0
                    for (q in 1:n.quantiles)
                    {
                        for (b in 1:B)
                        {
                            res.cur.mat <- list.obj[[l]]$boot.results.quantiles[[q]]$avg.outcomes[b,,]
                            cur.idx <- c(((b - 1) * n.entries + 1):(b * n.entries)) + ct
                            res.2.plot[cur.idx, 1] <- rep(colnames(res.cur.mat),
                                                          each = ncol(res.cur.mat))
                            res.2.plot[cur.idx, 2] <- gsub("^Received ", "", rep(rownames(res.cur.mat),
                                                          ncol(res.cur.mat)))
                            res.2.plot[cur.idx, 3] <- as.vector(res.cur.mat)
                            res.2.plot[cur.idx, 4] <- quantile.names[q]

                        }
                        ct <- ct + B * n.entries
                    }
                }

            }
            res.2.plot$Received <- as.factor(res.2.plot$Received)

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
            geom_density(alpha = 0.65, na.rm = TRUE) +
            geom_rug(aes(colour = Received), alpha = 0.85, sides = "l", na.rm = TRUE) +
            coord_flip() +
            facet_grid(Recommended ~ Model) +
            theme(legend.position = "bottom") +
            xlab("Outcome")
        if (avg.line)
        {
            pl.obj <- pl.obj + geom_vline(data = avg.res.2.plot.dens,
                                          aes(xintercept = Value),
                                          size = 1.25) +
                geom_vline(data = avg.res.2.plot.dens,
                           aes(xintercept = Value, colour = Received))
        }
    } else if (type == "conditional")
    {
        if (all(is_fitted_obj))
        {
            if (all_binom)
            {
                pl.obj <- ggplot(res.2.plot,
                                 aes(x = bs, y = Outcome,
                                     group = factor(Received),
                                     color = factor(Received) )) +
                    geom_point(na.rm = TRUE) +
                    geom_smooth(method = "gam", method.args = list(family = "binomial"), na.rm = TRUE) +
                    facet_grid( ~ Model) +
                    theme(legend.position = "bottom") +
                    scale_color_discrete(name = "Received") +
                    ggtitle("Individual Observations by Treatment Group")
            } else
            {
                pl.obj <- ggplot(res.2.plot,
                                 aes(x = bs, y = Outcome,
                                     group = factor(Received),
                                     color = factor(Received) )) +
                    geom_point(na.rm = TRUE) +
                    geom_smooth(method = "gam", na.rm = TRUE) +
                    facet_grid( ~ Model) +
                    theme(legend.position = "bottom") +
                    scale_color_discrete(name = "Received") +
                    ggtitle("Individual Observations by Treatment Group")
            }
        } else
        {
            pl.obj <- ggplot(res.2.plot,
                             aes(x = Received, y = Value)) +
                geom_boxplot(aes(fill = Received), na.rm = TRUE) +
                geom_rug(aes(colour = Received), alpha = 0.85, sides = "l", na.rm = TRUE) +
                facet_grid(Quantile ~ Recommended + Model) +
                theme(legend.position = "none") +
                ylab("Average Outcome")
        }
    } else if (type == "boxplot")
    {
        if (all_binom)
        {
            res.2.plot$Value <- as.factor(res.2.plot$Value)
            pl.obj <- ggplot(res.2.plot,
                             aes(x = Received, fill = factor(Value) )) +
                geom_bar(position = "fill", na.rm = TRUE) +
                facet_grid(Recommended ~ Model) +
                theme(legend.position = "bottom") +
                ylab("Outcome")
                guides(fill = guide_legend(title = "Observed Response"))
        } else
        {
            pl.obj <- ggplot(res.2.plot,
                             aes(x = Received, y = Value)) +
                geom_boxplot(aes(fill = Received), na.rm = TRUE) +
                geom_rug(aes(colour = Received), alpha = 0.85, sides = "l", na.rm = TRUE) +
                facet_grid(Recommended ~ Model) +
                theme(legend.position = "none") +
                ylab("Outcome")
        }
    } else
    {
        pl.obj <- ggplot(avg.res.2.plot,
                         aes(x = Recommended, y = Value, group = Received)) +
            geom_line(aes(colour = Received), size = 1.25, na.rm = TRUE) +
            geom_point(aes(colour = Received), size = 2, na.rm = TRUE) +
            facet_grid( ~ Model) +
            theme(legend.position = "bottom") +
            scale_x_discrete(expand = c(0.25, 0.25)) +
            ylab("Average Outcome")
    }
    pl.obj
}
