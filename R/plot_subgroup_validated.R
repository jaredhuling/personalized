#' Plotting validation results for fitted subgroup identification models
#'
#' @description Plots validation results for estimated subgroup treatment effects
#'
#' @param x fitted object returned by \code{validate.subgroup()} or \code{fit.subgroup()} function
#' @param type type of plot. \code{"density"} results in a density plot for the results
#' across all observations (if \code{x} is from \code{fit.subgroup()}) or if \code{x} is from \code{validate.subgroup()}
#' across iterations of either the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"boxplot"} results in boxplots across all observations/iterations of either
#' the bootstrap or training/test re-fitting. For the latter
#' case the test results will be plotted. \code{"interaction"} creates an
#' interaction plot for the different subgroups (crossing lines here means a meaningful subgroup). For the interaction plot,
#' the intervals around each point represent +1 one SE
#' \code{"conditional"} For subgroup_fitted objects, plots smoothed (via a GAM smoother) means of the outcomes as a function of the estimated benefit score
#' separately for the treated and untreated groups. For subgroup_validated objects, boxplots of summary statistics
#' within subgroups will be plotted as subgroups are defined by different cutoffs of the benefit scores.
#' These cutoffs can be specified via the \code{benefit.score.quantiles} argument of
#' \code{\link[personalized]{validate.subgroup}}.
#' @param avg.line boolean value of whether or not to plot a line for the average
#' value in addition to the density (only valid for \code{type = "density"})
#' @param ... not used
#' @seealso \code{\link[personalized]{validate.subgroup}} for function which creates validation results
#' and \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname plot
#' @import plotly
#' @importFrom ggplot2 ggplot aes geom_density geom_rug coord_flip facet_grid theme xlab
#' @importFrom ggplot2 ylab ggtitle geom_vline geom_boxplot geom_line geom_point geom_smooth geom_errorbar
#' @importFrom ggplot2 scale_x_discrete scale_color_discrete geom_histogram geom_rect geom_hline xlim geom_bar
#'
#' @examples
#'
#' valmod <- validate.subgroup(subgrp.model, B = 3,
#'                           method = "training_test",
#'                           benefit.score.quantiles = c(0.25, 0.5, 0.75),
#'                           train.fraction = 0.75)
#'
#' plot(valmod)
#'
#'
#' plot(valmod, type = "interaction")
#'
#' # see how summary statistics of subgroups change
#' # when the subgroups are defined based on different cutoffs
#' # (25th quantile of bene score, 50th, and 75th)
#' plot(valmod, type = "conditional")
#'
#' # visualize the frequency of particular variables
#' # of being selected across the resampling iterations with
#' # 'type = "stability"'
#' # not run:
#' # plot(valmod, type = "stability")
#'
#' @export
plot.subgroup_validated <- function(x,
                                    type = c("boxplot", "density", "interaction", "conditional", "stability"),
                                    avg.line = TRUE,
                                    ...)
{
    type <- match.arg(type)

    family   <- x$family

    avg.line <- as.logical(avg.line[1])

    boot.res <- x$boot.results$avg.outcomes
    avg.res  <- x$avg.results

    boot.dims <- dim(boot.res)

    n.entries <- prod(boot.dims[2:3])
    B <- boot.dims[1]

    avg.res.2.plot <- data.frame(Recommended = rep(colnames(avg.res$avg.outcomes),
                                                   each = ncol(avg.res$avg.outcomes)),
                                 Received    = gsub("^Received ", "", rep(rownames(avg.res$avg.outcomes),
                                                   ncol(avg.res$avg.outcomes))),
                                 Value       = as.vector(avg.res$avg.outcomes),
                                 SE          = as.vector(x$se.results$SE.avg.outcomes))



    avg.res.2.plot.dens <- avg.res.2.plot

    if (type == "interaction")
    {
        avg.res.2.plot$Recommended <- gsub("^Recommended ", "", avg.res.2.plot$Recommended)
    }

    avg.res.2.plot$Received <- as.factor(avg.res.2.plot$Received)
    avg.res.2.plot.dens$Received <- as.factor(avg.res.2.plot.dens$Received)

    avg.res.2.plot$Recommended <- as.factor(avg.res.2.plot$Recommended)
    avg.res.2.plot.dens$Recommended <- as.factor(avg.res.2.plot.dens$Recommended)

    ## reorder factors how they were ordered originally
    avg.res.2.plot$Received <- factor(avg.res.2.plot$Received,
                                      levels = levels(avg.res.2.plot$Received)[match(x$trts, sort(x$trts))])

    avg.res.2.plot$Recommended <- factor(avg.res.2.plot$Recommended,
                                         levels = levels(avg.res.2.plot$Recommended)[match(x$trts, sort(x$trts))])

    avg.res.2.plot.dens$Received <- factor(avg.res.2.plot.dens$Received,
                                           levels = levels(avg.res.2.plot.dens$Received)[match(x$trts, sort(x$trts))])

    Recommended <- Received <- Value <- bs <- Quantile <- Outcome <- SE <- NULL

    if (type == "conditional")
    {
        n.quantiles    <- length(x$boot.results.quantiles)
        quantile.names <- paste("Cutoff:", names(x$boot.results.quantiles))

        res.2.plot <- array(NA, dim = c(B * n.entries * n.quantiles, 4))
        colnames(res.2.plot) <- c("Recommended", "Received", "Value", "Quantile")
        res.2.plot <- data.frame(res.2.plot)

        ct <- 0
        for (q in 1:n.quantiles)
        {
            for (b in 1:B)
            {
                res.cur.mat <- x$boot.results.quantiles[[q]]$avg.outcomes[b,,]
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
    } else
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
    }

    res.2.plot$Received <- as.factor(res.2.plot$Received)
    res.2.plot$Recommended <- as.factor(res.2.plot$Recommended)

    res.2.plot$Received <- factor(res.2.plot$Received,
                                  levels = levels(res.2.plot$Received)[match(x$trts, sort(x$trts))])

    res.2.plot$Recommended <- factor(res.2.plot$Recommended,
                                     levels = levels(res.2.plot$Recommended)[match(x$trts, sort(x$trts))])


    title.text <- NULL
    if (x$val.method == "training_test_replication")
    {
        title.text <- "Average Test Set Outcome Across Replications Among Subgroups"
    } else
    {
        title.text <- "Average Bias-Corrected Outcome Across Replications Among Subgroups"
    }

    ylab.text <- "Average Outcome"

    if (family == "cox")
    {
        ylab.text <- "Average Restricted Mean Survival"
    }

    if (type == "density")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Value, fill = Received)) +
            geom_density(alpha = 0.65, na.rm = TRUE) +
            geom_rug(aes(colour = Received), alpha = 0.85, na.rm = TRUE, sides = "l") +
            coord_flip() +
            facet_grid( ~ Recommended) +
            theme(legend.position = "bottom") +
            xlab(ylab.text) +
            ggtitle(title.text)
        if (avg.line)
        {
            pl.obj <- pl.obj + geom_vline(data = avg.res.2.plot.dens,
                                          aes(xintercept = Value),
                                          size = 1.25) +
                geom_vline(data = avg.res.2.plot.dens,
                           aes(xintercept = Value, colour = Received))
        }
    } else if (type == "boxplot")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Received, y = Value)) +
            geom_boxplot(aes(fill = Received), na.rm = TRUE) +
            geom_rug(aes(colour = Received), alpha = 0.85, na.rm = TRUE, sides = "l") +
            facet_grid( ~ Recommended) +
            theme(legend.position = "none") +
            ylab(ylab.text) +
            ggtitle(title.text)
    } else if (type == "conditional")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Received, y = Value)) +
            geom_boxplot(aes(fill = Received), na.rm = TRUE) +
            geom_rug(aes(colour = Received), alpha = 0.85, na.rm = TRUE, sides = "l") +
            facet_grid(Recommended ~ Quantile) +
            theme(legend.position = "none") +
            ylab(ylab.text) +
            ggtitle(title.text)
    } else if (type == "stability")
    {
      # Acquire coefficients for each bootstrap iteration (exclude Intercept and Trt terms)
        coefmat <- x$boot.results$coefficients

        if (!is.null(rownames(coefmat)) &&
            "(Intercept)" == rownames(coefmat)[1])
        {
            coefmat <- coefmat[-c(1,2),]
        } else
        {
            coefmat <- coefmat[-c(1),]
        }

        if (is.null(rownames(coefmat)))
        {
            rownames(coefmat) <- paste0("V", 1:NROW(coefmat))
        }

        d <- as.data.frame(coefmat)

        # Variables to be created in this code block
        pct.selected <- signs <- is.consistent <- summary.stats <- min.stat  <- med.stat <- max.stat <-
            bar.type <- name <- plot.idx <- Variable <- Selection <- Median <- Range <- p.primary <- p.secondary <- NULL

        # Compute percentage of times each variable was selected
        d$pct.selected <- apply(d,1,function(x){sum(x!=0)}/ncol(d)*100)

        # Remove instances where variables were never selected in any bootstrap iteration
        d <- subset(d, pct.selected != 0)

        # Compute percentage of time variable has consistent sign.
        # A variable is deemed consistent if it has the same sign at least 95% of the times it was selected.
        signs <- apply(d[,grep("B",colnames(d), value = TRUE)], 1, function(x){sign(x)[x!=0]})

        if (is.matrix(signs))
        {
            signsmat <- signs
            signs <- vector(mode = "list", length = NCOL(signsmat))
            names(signs) <- colnames(signsmat)
            for (j in 1:NCOL(signsmat))
            {
                signs[[j]] <- signsmat[,j]
                names(signs[[j]]) <- rownames(signsmat[,j])
            }
        }

        d$is.consistent <- sapply(signs,function(x){any(table(x) / length(x) >= .95)})

        # Calculate min, median, and max
        summary.stats <- apply(d[,grep("B",colnames(d), value=TRUE)],1,function(x){summary(x[x!=0])})
        d$min.stat <- summary.stats["Min.",]
        d$med.stat <- summary.stats["Median",]
        d$max.stat <- summary.stats["Max.",]

        # Create label for bar type (Positive/Negative Tendency or Mixed)
        d$bar.type <- factor(ifelse(d$is.consistent, ifelse(d$med.stat > 0,"Positive Tendency","Negative Tendency"),"Mixed"),
                             levels=c("Negative Tendency", "Mixed", "Positive Tendency"))

        # Order by most frequently selected and bar type
        d <- d[order(d$bar.type,-d$pct.selected, -abs(d$med.stat) ),]

        # Add variable name and plot index to data for plotting purposes
        d$name <- rownames(d)
        d$plot.idx <- 1:nrow(d)

        # Remove individual bootstrap values from plotting data frame
        d <- d[,!(names(d) %in% grep("B",names(d),value=TRUE))]

        # Create tooltip statistics
        d$Variable <- d$name
        d$Selection <- paste0(d$pct.selected,"%")
        d$Median <- round(d$med.stat,4)
        d$Range <- paste0("[",round(d$min.stat,4),",",round(d$max.stat,4),"]")

        # Primary Plot - Range with median points
        p.primary <- ggplot(d, aes(xmin = plot.idx-0.5, xmax=plot.idx+0.5, ymin = min.stat, ymax = max.stat, x=plot.idx, y = med.stat, fill = bar.type,
                                   tooltip1 = Variable, tooltip2 = Selection, tooltip3 = Median, tooltip4 = Range )) +
            geom_rect(color="black", stat="identity") +
            geom_point(size=2, shape=21, color="black", fill="azure1", stat="identity") +
            geom_hline(yintercept = 0) +
            geom_vline(xintercept = c(which.min(d$bar.type=="Negative Tendency") - 0.5, which.max(d$bar.type=="Positive Tendency") - 0.5), linetype = "dashed") +
            xlim(0,nrow(d)+1)

        # Secondary Plot - Distribution of selection probability
        p.secondary <- ggplot(d, aes(x = plot.idx, y = pct.selected, fill = bar.type,
                                     tooltip1 = Variable, tooltip2 = Selection)) +
            geom_bar(stat="identity") +
            geom_vline(xintercept = c(which.min(d$bar.type=="Negative Tendency") - 0.5, which.max(d$bar.type=="Positive Tendency") - 0.5), linetype = "dashed")

        # Combine plots and create plotly object
        pl.obj <-
            subplot(ggplotly(p.primary, tooltip = paste0("tooltip",1:4)),
                    ggplotly(p.secondary, tooltip = paste0("tooltip",1:2)),
                    nrows=2,
                    shareX = TRUE,
                    titleX = TRUE,
                    titleY = TRUE
            ) %>%
            layout(title="Variable Selection Across Bootstrap Iterations",
                   showlegend=FALSE,
                   xaxis =  list(title = "Plot Index"),
                   yaxis =  list(title = "Coefficient Value"),
                   yaxis2 = list(title = "Percent of Times Selected")
            )
    } else
    {
        pl.obj <- ggplot(avg.res.2.plot,
                         aes(x = Recommended, y = Value, group = Received)) +
            geom_line(aes(colour = Received), size = 1.25, na.rm = TRUE) +
            geom_point(aes(colour = Received), size = 2, na.rm = TRUE) +
            geom_errorbar(aes(ymin   = Value - SE,
                              ymax   = Value + SE,
                              colour = Received),
                          width = 0.2) +
            theme(legend.position = "bottom") +
            scale_x_discrete(expand = c(0.25, 0.25)) +
            ylab(ylab.text) +
            ggtitle(title.text)
    }
    # Return plot
    pl.obj

}
