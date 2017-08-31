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
#' interaction plot for the different subgroups (crossing lines here means a meaningful subgroup)
#' @param avg.line boolean value of whether or not to plot a line for the average
#' value in addition to the density (only valid for \code{type = "density"})
#' @param ... not used
#' @seealso \code{\link[personalized]{validate.subgroup}} for function which creates validation results
#' and \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname plot
#' @import ggplot2
#'
#' @examples
#'
#' valmod <- validate.subgroup(subgrp.model, B = 5,
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
                                    type = c("boxplot", "density", "interaction", "stability"),
                                    avg.line = TRUE,
                                    ...)
{
    type <- match.arg(type)

    family   <- x$family

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
        ylab.text <- "Average Restricted Mean"
    }

    if (type == "density")
    {
        pl.obj <- ggplot(res.2.plot,
                         aes(x = Value, fill = Received)) +
                      geom_density(alpha = 0.65) +
                      geom_rug(aes(colour = Received), alpha = 0.85) +
                      coord_flip() +
                      facet_grid( ~ Recommended) +
                      theme(legend.position = "bottom") +
            xlab(ylab.text) +
            ggtitle(title.text)
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
            theme(legend.position = "bottom") +
            ylab(ylab.text) +
            ggtitle(title.text)
    } else if (type == "stability")
    {
      # Acquire coefficients for each bootstrap iteration (exclude intercept and Trt terms)
      d <- as.data.frame(validation$boot.results[[4]][-c(1,2),])
      
      # Compute percentage of times each variable was selected
      d$pct.selected <- apply(d,1,function(x){sum(x!=0)}/ncol(d)*100)
      
      # Calculate quartiles
      d$q0 <- apply(d[,grep("B",colnames(d), value=T)],1,function(x){quantile(x[x!=0],0)})
      #d$q1 <- apply(d[,grep("B",colnames(d), value=T)],1,function(x){quantile(x[x!=0],0.25)})
      d$q2 <- apply(d[,grep("B",colnames(d), value=T)],1,function(x){quantile(x[x!=0],0.50)})
      #d$q3 <- apply(d[,grep("B",colnames(d), value=T)],1,function(x){quantile(x[x!=0],0.75)})
      d$q4 <- apply(d[,grep("B",colnames(d), value=T)],1,function(x){quantile(x[x!=0],1)})
      
      # Compute color for bars (Both positive/negative, strictly positive, or strictly negative)
      d$col <- ifelse(d$q4 > 0 & d$q0 < 0, "green", ifelse(d$q4 > 0,"blue","red"))
      
      # Order by color and then descending by group median
      d <- d[order(d$col,-d$q2),]
      
      # Remove instances where variables were never selected
      d <- subset(d, pct.selected != 0)
      
      # Capture names for tooltip
      d$name <- rownames(d)
      
      # Assign variable index
      d$plot.idx <- 1:nrow(d)
      
      # Trim to keep only necessary columns
      d <- d[,!(names(d) %in% grep("B",names(d),value=T))]
      # Construct Plot
      pl.obj <- ggplot(data = d) +
      geom_rect(mapping = aes(text=name, xmin = plot.idx-0.5, xmax=plot.idx+0.5, ymin = q0, ymax = q4, fill = col), color="black", stat="identity") +
      geom_point(mapping = aes(text=name, x = plot.idx, y = q2, size=pct.selected), shape=21, color="black", fill="purple", stat="identity") +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = c(which.min(d$col=="blue") - 0.5, which.max(d$col=="red") - 0.5), linetype = "dashed") +
      xlab("Variable Index") +
      ylab("Coefficient Value") +
      ggtitle("Validation Stability") +
      xlim(0,nrow(d)+1) +
      theme(plot.title = element_text(hjust = 0.5))
    } else
    {
        pl.obj <- ggplot(avg.res.2.plot,
                         aes(x = Recommended, y = Value, group = Received)) +
            geom_line(aes(colour = Received), size = 1.25) +
            geom_point(aes(colour = Received), size = 2) +
            theme(legend.position = "bottom") +
            scale_x_discrete(expand = c(0.25, 0.25)) +
            ylab(ylab.text) +
            ggtitle(title.text)
    }
    pl.obj
}
