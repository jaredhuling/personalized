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
    cat("method: ", x$method, "\n\n")
    Cf <- matrix(paste0(round(x$subgroup.trt.effects$avg.outcomes, digits),
                        " (n = ", x$subgroup.trt.effects$sample.sizes, ")"), ncol = 2)
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

    cat("\nBenefit score quantiles: \n")
    print(quantile(x$benefit.scores), digits = digits)
}

#' @rdname print
#' @export
print.subgroup_validated <- function(x, digits = max(getOption('digits')-3, 3), ...)
{
    cat("family: ", x$family, "\n")
    cat("loss:   ", x$loss, "\n")
    cat("method: ", x$method, "\n\n")
    cat("validation method: ", x$val.method, "\n\n")

    if (x$val.method == "training_test_replication")
    {
        valtext <- "Test Set Outcomes:"
    } else
    {
        valtext <- "Bootstrap Bias-Corrected Outcomes:"
    }

    cat(paste0("Average ", valtext, "\n"))

    Cf <- matrix(paste0(round(x$avg.results$avg.outcomes, digits),
                        " (SE = ", round(x$se.results$SE.avg.outcomes, digits), ")"), ncol = 2)
    dimnames(Cf) <- dimnames(x$avg.results$avg.outcomes)

    print.default(Cf, quote = FALSE, right = TRUE, na.print = "NA",
                  ...)

    cat("\n")

    Cf2 <- paste0(round(x$avg.results$subgroup.effects, digits),
                         " (SE = ", round(x$se.results$SE.subgroup.effects, digits), ")")
    names(Cf2) <- names(x$avg.results$subgroup.effects)
    print.default(Cf2, quote = FALSE, right = TRUE, na.print = "NA",
                  ...)
}
