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
    cat("family:   ", x$family, "\n")
    cat("loss:     ", x$loss, "\n")
    cat("method:   ", x$method, "\n")
    cat("cutpoint: ", x$cutpoint, "\n")

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
        cat("propensity \nfunction: ",
            func.name,
            "\n")
    }
    cat("\n")
    if (x$n.trts == 2)
    {
        if (x$larger.outcome.better)
        {
            cat("benefit score: f(x), \nTrt recom =",
                paste0(x$comparison.trts, "*I(f(x)>c)+",
                       x$reference.trt, "*I(f(x)<=c)"), "where c is 'cutpoint'\n")
        } else
        {
            cat("benefit score: f(x), \nTrt recom =",
                paste0(x$comparison.trts, "*I(f(x)<c)+",
                       x$reference.trt, "*I(f(x)>=c)"), "where c is 'cutpoint'\n")
        }

    } else if (x$n.trts == NCOL(x$benefit.scores)) ## multinomial model case
    {
        if (x$larger.outcome.better)
        {
            bene.text <- paste(paste0("f_", (x$trts), "(x): ",
                                      paste(x$trts)),
                               collapse = ",  ")
            cat("benefit scores:", bene.text, "\n")
            trt.rec.text <- paste0("maxval = max(", paste(paste0("f_", (x$trts), "(x)"), collapse = ", "), ")")
            cat(trt.rec.text, "\n")
            cat("which.max(maxval) = The trt level which maximizes maxval\n")
            cat("Trt recom = which.max(maxval)\n")
        } else
        {
            bene.text <- paste(paste0("f_", (x$trts), "(x): ",
                                      paste(x$trts)),
                               collapse = ",  ")
            cat("benefit scores:", bene.text, "\n")
            trt.rec.text <- paste0("maxval = min(", paste(paste0("f_", (x$trts), "(x)"), collapse = ", "), ")")
            cat(trt.rec.text, "\n")
            cat("which.min(minval) = The trt level which mininizes minval\n")
            cat("Trt recom = which.min(minval)\n")
        }
    } else
    {
        if (x$larger.outcome.better)
        {
            bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                      paste(x$comparison.trts, "vs", x$reference.trt)),
                               collapse = ",  ")
            bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
            cat("benefit score:", bene.text, bene.txt2, "\n")
            trt.rec.text <- paste0("maxval = max(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
            cat(trt.rec.text, "\n")
            cat("which.max(maxval) = The trt level which maximizes maxval\n")
            cat("Trt recom = which.max(maxval)*I(maxval > c) +", paste0(x$reference.trt, "*I(maxval <= c) where c is 'cutpoint'\n") )
        } else
        {
            bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                      paste(x$comparison.trts, "vs", x$reference.trt)),
                               collapse = ",  ")
            bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
            cat("benefit score:", bene.text, bene.txt2, "\n")
            trt.rec.text <- paste0("minval = min(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
            cat(trt.rec.text, "\n")
            cat("which.min(minval) = The trt level which mininizes minval\n")
            cat("Trt recom = which.min(minval)*I(minval < c) +", paste0(x$reference.trt, "*I(minval >= c) where c is 'cutpoint'\n") )
        }
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
    cat("Treatment effects conditional on subgroups:\n")
    print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                  ...)

    ncol.bs <- NCOL(x$benefit.scores)

    cat("\nNOTE: The above average outcomes are biased estimates of\n      the expected outcomes conditional on subgroups. \n      Use 'validate.subgroup()' to obtain unbiased estimates.\n")

    cat("\n---------------------------------------------------\n")

    if (is.null(ncol.bs) || ncol.bs == 1)
    {
        cname <- paste0(x$comparison.trts, " vs ", x$reference.trt)
        cat("\nBenefit score", paste0("quantiles (f(X) for ", cname, "): \n") )
        print(quantile(x$benefit.scores), digits = digits)
    } else if (ncol.bs == length(x$trts))
    {
        for (cc in 1:ncol.bs)
        {
            cname <- paste0(x$trts[cc])
            cat("\nBenefit score", cc, paste0("quantiles (f(X) for ", cname, "): \n") )
            print(quantile(x$benefit.scores[,cc]), digits = digits)
        }
    } else
    {
        for (cc in 1:ncol.bs)
        {
            cname <- paste0(x$comparison.trts[cc], " vs ", x$reference.trt)
            cat("\nBenefit score", cc, paste0("quantiles (f(X) for ", cname, "): \n") )
            print(quantile(x$benefit.scores[,cc]), digits = digits)
        }
    }

    cat("\n---------------------------------------------------\n")
    cat("\n")
    print.individual_treatment_effects(x$individual.trt.effects, digits = digits, ...)
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
            cat("cutpoint:          ", names(x$boot.results.quantiles)[q], "\n")

            if (is.null(x$iterations))
            {
                iters <- length(x$boot.results$overall.subgroup.effect)
            } else
            {
                iters <- x$iterations
            }
            cat("replications:      ", iters, "\n\n")

            if (x$val.method == "training_test_replication")
            {
                valtext <- "Test Set Outcomes:"
            } else
            {
                valtext <- "Bootstrap Bias-Corrected Outcomes:"
            }



            if (x$n.trts == 2)
            {
                if (x$larger.outcome.better)
                {
                    cat("benefit score: f(x), \nTrt recom =",
                        paste0(x$comparison.trts, "*I(f(x)>c)+",
                               x$reference.trt, "*I(f(x)<=c)"), "where c is 'cutpoint'\n\n")
                } else
                {
                    cat("benefit score: f(x), \nTrt recom =",
                        paste0(x$comparison.trts, "*I(f(x)<c)+",
                               x$reference.trt, "*I(f(x)>=c)"), "where c is 'cutpoint'\n\n")
                }

            } else
            {
                if (x$larger.outcome.better)
                {
                    bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                              paste(x$comparison.trts, "vs", x$reference.trt)),
                                       collapse = ",  ")
                    bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
                    cat("benefit score:", bene.text, bene.txt2, "\n")
                    trt.rec.text <- paste0("maxval = max(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
                    cat(trt.rec.text, "\n")
                    cat("which.max(maxval) = The trt level which maximizes maxval\n")
                    cat("Trt recom = which.max(maxval)*I(maxval > c) +", paste0(x$reference.trt, "*I(maxval <= c) where c is 'cutpoint'\n") )
                } else
                {
                    bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                              paste(x$comparison.trts, "vs", x$reference.trt)),
                                       collapse = ",  ")
                    bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
                    cat("benefit score:", bene.text, bene.txt2, "\n")
                    trt.rec.text <- paste0("minval = min(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
                    cat(trt.rec.text, "\n")
                    cat("which.min(minval) = The trt level which minimizes minval\n")
                    cat("Trt recom = which.min(minval)*I(minval < c) +", paste0(x$reference.trt, "*I(minval >= c) where c is 'cutpoint'\n") )
                }
                cat("\n")
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

            cat("Treatment effects conditional on subgroups:\n")
            print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                          ...)

            cat("\n")

            overall <- paste0(round(x$avg.quantile.results[[q]]$overall.subgroup.effect, digits),
                              " (SE = ", round(x$se.quantile.results[[q]]$SE.overall.subgroup.effect, digits), ")")
            #names(overall) <- "Overall treatment effect conditional on subgroups (E[Y|Trt Received = Trt recommended])"
            cat("Est of E[Y|Trt received = Trt recom] - E[Y|Trt received =/= Trt recom]:")
            names(overall) <- ""

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
        cat("cutpoint:          ", x$cutpoint, "\n")

        if (is.null(x$iterations))
        {
            iters <- length(x$boot.results$overall.subgroup.effect)
        } else
        {
            iters <- x$iterations
        }
        cat("replications:      ", iters, "\n\n")

        if (x$val.method == "training_test_replication")
        {
            valtext <- "Test Set Outcomes:"
        } else
        {
            valtext <- "Bootstrap Bias-Corrected Outcomes:"
        }

        if (x$n.trts == 2)
        {
            if (x$larger.outcome.better)
            {
                cat("benefit score: f(x), \nTrt recom =",
                    paste0(x$comparison.trts, "*I(f(x)>c)+",
                           x$reference.trt, "*I(f(x)<=c)"), "where c is 'cutpoint'\n\n")
            } else
            {
                cat("benefit score: f(x), \nTrt recom =",
                    paste0(x$comparison.trts, "*I(f(x)<c)+",
                           x$reference.trt, "*I(f(x)>=c)"), "where c is 'cutpoint'\n\n")
            }

        } else
        {
            if (x$larger.outcome.better)
            {
                bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                          paste(x$comparison.trts, "vs", x$reference.trt)),
                                   collapse = ",  ")
                bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
                cat("benefit score:", bene.text, bene.txt2, "\n")
                trt.rec.text <- paste0("maxval = max(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
                cat(trt.rec.text, "\n")
                cat("which.max(maxval) = The trt level which maximizes maxval\n")
                cat("Trt recom = which.max(maxval)*I(maxval > c) +", paste0(x$reference.trt, "*I(maxval <= c) where c is 'cutpoint'\n") )
            } else
            {
                bene.text <- paste(paste0("f_", (x$comparison.trts), "(x): ",
                                          paste(x$comparison.trts, "vs", x$reference.trt)),
                                   collapse = ",  ")
                bene.txt2 <- paste0("\n               f_", x$reference.trt,  "(x): 0")
                cat("benefit score:", bene.text, bene.txt2, "\n")
                trt.rec.text <- paste0("minval = min(", paste(paste0("f_", (x$comparison.trts), "(x)"), collapse = ", "), ")")
                cat(trt.rec.text, "\n")
                cat("which.min(minval) = The trt level which minimizes minval\n")
                cat("Trt recom = which.min(minval)*I(minval < c) +", paste0(x$reference.trt, "*I(minval >= c) where c is 'cutpoint'\n") )
            }
            cat("\n")
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
        cat("Treatment effects conditional on subgroups:\n")
        print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                      ...)

        cat("\n")

        overall <- paste0(round(x$avg.results$overall.subgroup.effect, digits),
                          " (SE = ", round(x$se.results$SE.overall.subgroup.effect, digits), ")")
        #names(overall) <- "Overall treatment effect conditional on subgroups"
        cat("Est of \nE[Y|Trt received = Trt recom] - E[Y|Trt received =/= Trt recom]:")
        names(overall) <- ""

        print.default(overall, quote = FALSE, right = FALSE, na.print = "NA",
                      ...)
    }

}
