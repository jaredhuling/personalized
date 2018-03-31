#' Printing results for fitted subgroup identification models
#'
#' @description Prints results for estimated subgroup treatment effects
#'
#' @param x a fitted object from either \code{fit.subgroup}, \code{validate.subgroup}, or \code{summarize.subgroups()}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @seealso \code{\link[personalized]{validate.subgroup}} for function which creates validation results
#' and \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname print
#' @importFrom stats quantile
#' @export
print.subgroup_fitted <- function(x, digits = max(getOption('digits')-3, 3), ...)
{
    cat("family: ", x$family, "\n")
    cat("loss:   ", x$loss, "\n")
    cat("method: ", x$method, "\n")

    if (!is.null(x$augment.func))
    {
        func.name <- as.character(substitute(x$augment.func))
        func.name <- func.name[length(func.name)]
        cat("augmentation \nfunction:",
            func.name,
            "\n")
    }
    if (!is.null(x$propensity.func))
    {
        func.name <- as.character(substitute(x$propensity.func))
        func.name <- func.name[length(func.name)]
        cat("propensity \nfunction:",
            func.name,
            "\n")
    }

    cat("\n")

    Cf <- matrix(paste0(round(x$subgroup.trt.effects$avg.outcomes, digits),
                        " (n = ", x$subgroup.trt.effects$sample.sizes, ")"), ncol = ncol(x$subgroup.trt.effects$avg.outcomes))
    dimnames(Cf) <- dimnames(x$subgroup.trt.effects$avg.outcomes)

    cat("Average Outcomes:\n")
    print.default(Cf, quote = FALSE, right = TRUE, na.print = "NA",
                  ...)

    cat("\n")

    Cf2 <- paste0(round(x$subgroup.trt.effects$subgroup.effects, digits),
                               " (n = ", colSums(x$subgroup.trt.effects$sample.sizes), ")")
    names(Cf2) <- names(x$subgroup.trt.effects$subgroup.effects)
    print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                  ...)

    ncol.bs <- NCOL(x$benefit.scores)

    if (is.null(ncol.bs) || ncol.bs == 1)
    {
        cat("\nBenefit score quantiles: \n")
        print(quantile(x$benefit.scores), digits = digits)
    } else
    {
        for (cc in 1:ncol.bs)
        {
            cat("\nBenefit score", cc, "quantiles: \n")
            print(quantile(x$benefit.scores[,cc]), digits = digits)
        }
    }
}

#' @param sample.pct boolean variable of whether to print the percent of the test sample within each subgroup. If false
#' the sample size itself, not the percent is printed. This may not be informative if the test sample size is much different
#' from the total sample size
#' @param which.quant when \code{validate.subgroup()} is called with a vector of quantile values specified for \code{benefit.score.quantiles},
#' i.e. \code{benefit.score.quantiles = c(0.25, 0.5, 0.75)}, the argument \code{which.quant} can be a vector of indexes specifying which
#' quantile cutoff value validation results to display, i.e. \code{which.quant = c(1,3)} in the above example results in the display of
#' validation results for subgroups defined by cutoff values of the benefit score defined by the 25th abnd 75th quantiles of the benefit score
#' @rdname print
#' @export
print.subgroup_validated <- function(x, digits = max(getOption('digits')-3, 3), sample.pct = FALSE, which.quant = NULL, ...)
{

    if (!is.null(which.quant))
    {
        if (!is.vector(which.quant) | !is.numeric(which.quant) | any(which.quant <= 0))
        {
            stop("which.quant must be a vector of positive integers")
        }

        nq <- length(which.quant)

        which.quant <- as.integer(which.quant)

        if (nq > length(x$boot.results.quantiles) | max(which.quant) > length(x$boot.results.quantiles))
        {
            stop("trying to request a non-existent quantile. Try a smaller index for which.quant")
        }

        ct <- 0
        for (q in which.quant)
        {
            ct <- ct + 1
            cat("family: ", x$family, "\n")
            cat("loss:   ", x$loss, "\n")
            cat("method: ", x$method, "\n\n")
            cat("validation method: ", x$val.method, "\n")
            cat("Cutoff:            ", names(x$boot.results.quantiles)[q], "\n")

            if (is.null(x$iterations))
            {
                iters <- length(x$boot.results$overall.subgroup.effect)
            } else
            {
                iters <- x$iterations
            }
            cat("iterations:        ", iters, "\n\n")

            if (x$val.method == "training_test_replication")
            {
                valtext <- "Test Set Outcomes:"
            } else
            {
                valtext <- "Bootstrap Bias-Corrected Outcomes:"
            }

            cat(paste0("Average ", valtext, "\n"))

            if (sample.pct)
            {
                Cf <- matrix(paste0(round(x$avg.quantile.results[[q]]$avg.outcomes, digits),
                                    " (SE = ", round(x$se.quantile.results[[q]]$SE.avg.outcomes, digits),
                                    ", ", round(100 * x$avg.quantile.results[[q]]$sample.sizes /
                                                    sum(x$avg.quantile.results[[q]]$sample.sizes),
                                                digits), "%)"),
                             ncol = ncol(x$avg.quantile.results[[q]]$avg.outcomes))
            } else
            {
                Cf <- matrix(paste0(round(x$avg.quantile.results[[q]]$avg.outcomes, digits),
                                    " (SE = ", round(x$se.quantile.results[[q]]$SE.avg.outcomes, digits),
                                    ", n = ", round(x$avg.quantile.results[[q]]$sample.sizes, digits), ")"),
                             ncol = ncol(x$avg.quantile.results[[q]]$avg.outcomes))
            }
            dimnames(Cf) <- dimnames(x$avg.quantile.results[[q]]$avg.outcomes)

            print.default(Cf, quote = FALSE, right = TRUE, na.print = "NA",
                          ...)

            cat("\n")

            if (sample.pct)
            {
                Cf2 <- paste0(round(x$avg.quantile.results[[q]]$subgroup.effects, digits),
                              " (SE = ", round(x$se.quantile.results[[q]]$SE.subgroup.effects, digits),
                              ", ", round(100 * colSums(x$avg.quantile.results[[q]]$sample.sizes) /
                                              sum(x$avg.quantile.results[[q]]$sample.sizes), digits), "%)")
            } else
            {
                Cf2 <- paste0(round(x$avg.quantile.results[[q]]$subgroup.effects, digits),
                              " (SE = ", round(x$se.quantile.results[[q]]$SE.subgroup.effects, digits),
                              ", n = ", round(colSums(x$avg.quantile.results[[q]]$sample.sizes), digits), ")")
            }
            names(Cf2) <- names(x$avg.quantile.results[[q]]$subgroup.effects)
            print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                          ...)

            cat("\n")

            overall <- paste0(round(x$avg.quantile.results[[q]]$overall.subgroup.effect, digits),
                              " (SE = ", round(x$se.quantile.results[[q]]$SE.overall.subgroup.effect, digits), ")")
            names(overall) <- "Overall Subgroup Effect"

            print.default(overall, quote = FALSE, right = TRUE, na.print = "NA",
                          ...)

            if (ct < nq)
            {
                cat("\n<===============================================>\n\n")
            }

        }
    } else
    {
        cat("family: ", x$family, "\n")
        cat("loss:   ", x$loss, "\n")
        cat("method: ", x$method, "\n\n")
        cat("validation method: ", x$val.method, "\n")

        if (is.null(x$iterations))
        {
            iters <- length(x$boot.results$overall.subgroup.effect)
        } else
        {
            iters <- x$iterations
        }
        cat("iterations: ", iters, "\n\n")

        if (x$val.method == "training_test_replication")
        {
            valtext <- "Test Set Outcomes:"
        } else
        {
            valtext <- "Bootstrap Bias-Corrected Outcomes:"
        }

        cat(paste0("Average ", valtext, "\n"))

        if (sample.pct)
        {
            Cf <- matrix(paste0(round(x$avg.results$avg.outcomes, digits),
                                " (SE = ", round(x$se.results$SE.avg.outcomes, digits),
                                ", ", round(100 * x$avg.results$sample.sizes / sum(x$avg.results$sample.sizes),
                                            digits), "%)"),
                         ncol = ncol(x$avg.results$avg.outcomes))
        } else
        {
            Cf <- matrix(paste0(round(x$avg.results$avg.outcomes, digits),
                                " (SE = ", round(x$se.results$SE.avg.outcomes, digits),
                                ", n = ", round(x$avg.results$sample.sizes, digits), ")"),
                         ncol = ncol(x$avg.results$avg.outcomes))
        }
        dimnames(Cf) <- dimnames(x$avg.results$avg.outcomes)

        print.default(Cf, quote = FALSE, right = TRUE, na.print = "NA",
                      ...)

        cat("\n")

        if (sample.pct)
        {
            Cf2 <- paste0(round(x$avg.results$subgroup.effects, digits),
                          " (SE = ", round(x$se.results$SE.subgroup.effects, digits),
                          ", ", round(100 * colSums(x$avg.results$sample.sizes) /
                                          sum(x$avg.results$sample.sizes), digits), "%)")
        } else
        {
            Cf2 <- paste0(round(x$avg.results$subgroup.effects, digits),
                          " (SE = ", round(x$se.results$SE.subgroup.effects, digits),
                          ", n = ", round(colSums(x$avg.results$sample.sizes), digits), ")")
        }
        names(Cf2) <- names(x$avg.results$subgroup.effects)
        print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                      ...)

        cat("\n")

        overall <- paste0(round(x$avg.results$overall.subgroup.effect, digits),
                          " (SE = ", round(x$se.results$SE.overall.subgroup.effect, digits), ")")
        names(overall) <- "Overall Subgroup Effect"

        print.default(overall, quote = FALSE, right = TRUE, na.print = "NA",
                      ...)
    }


}
