#' Plotting validation results for fitted subgroup identification models
#'
#' @description Plots validation results for estimated subgroup treatment effects
#'
#' @param x fitted object returned by \code{validate.subgrp()} or \code{fit.subgrp()} function
#' @param type type of plot. \code{"density"} results in a density plot for the results
#' across all observations (if \code{x} is from \code{fit.subgrp()}) or if \code{x} is from \code{validate.subgrp()}
#' across iterations of either the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"boxplot"} results in boxplots across all observations/iterations of either
#' the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"interaction"} creates an
#' interaction plot for the different subgroups (crossing lines here means a meaningful subgroup)
#' @param avg.line boolean value of whether or not to plot a line for the average
#' value in addition to the density (only valid for \code{type = "density"})
#' @param ... not used
#' @seealso \code{\link[personalized]{validate.subgrp}} for function which creates validation results
#' and \code{\link[personalized]{fit.subgrp}} for function which fits subgroup identification models.
#' @rdname plot
#' @import ggplot2
#'
#' @examples
#'
#' valmod <- validate.subgrp(subgrp.model, B = 10,
#'                           method = "training_test",
#'                           train.fraction = 0.75)
#' valmod$avg.results
#'
#' plot(valmod)
#'
#' plot(valmod, type = "boxplot")
#'
#' plot(valmod, type = "interaction")
#'
#' @export
plot.subgroup_validated <- function(x,
                                    type = c("density", "boxplot", "interaction"),
                                    avg.line = TRUE,
                                    ...)
{
    type <- match.arg(type)

    avg.line <- as.logical(avg.line[1])

    boot.res <- x$boot.results$avg.outcomes
    avg.res  <- x$avg.results
    B <- dim(boot.res)[1]

    res.2.plot <- array(NA, dim = c(B * 4, 3))
    colnames(res.2.plot) <- c("Recommended", "Received", "Value")
    res.2.plot <- data.frame(res.2.plot)

    avg.res.2.plot <- data.frame(Recommended = c("Recommended Trt", "Recommended Trt",
                                                 "Recommended Ctrl", "Recommended Ctrl"),
                                 Received    = c("Received Trt", "Received Ctrl",
                                                 "Received Trt", "Received Ctrl"),
                                 Value       = as.vector(avg.res$avg.outcomes))

    Recommended <- Received <- Value <- NULL

    for (b in 1:B)
    {
        cur.idx <- c(((b - 1) * 4 + 1):(b * 4))
        res.2.plot[cur.idx, 1] <- c("Recommended Trt", "Recommended Trt",
                                    "Recommended Ctrl", "Recommended Ctrl")
        res.2.plot[cur.idx, 2] <- c("Received Trt", "Received Ctrl",
                                    "Received Trt", "Received Ctrl")
        res.2.plot[cur.idx, 3] <- as.vector(boot.res[b,,])
    }

    if (type == "density")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Value, fill = Received)) +
                      geom_density(alpha = 0.65) +
                      geom_rug(aes(colour = Received), alpha = 0.85) +
                      coord_flip() +
                      facet_grid( ~ Recommended) +
                      theme(legend.position = "bottom")
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
            facet_grid( ~ Recommended) +
            theme(legend.position = "bottom")
    } else
    {
        pl.obj <- ggplot(avg.res.2.plot,
                         aes(x = Recommended, y = Value, group = Received)) +
            geom_line(aes(colour = Received), size = 1.25) +
            theme(legend.position = "bottom") +
            scale_x_discrete(expand = c(0.25, 0.25))
    }
    pl.obj
}

