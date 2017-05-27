
#' Summary of results for fitted subgroup identification models
#'
#' @description Prints summary of results for estimated subgroup treatment effects
#'
#' @param object a fitted object from either \code{fit.subgroup} or \code{validate.subgroup}
#' @param digits minimal number of significant digits to print.
#' @param ... further arguments passed to or from \code{\link[base]{print.default}}.
#' @seealso \code{\link[personalized]{validate.subgroup}} for function which creates validation results
#' and \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname summary
#' @export
summary.subgroup_fitted <- function(object, digits = max(getOption('digits')-3, 3), ...)
{
    print.subgroup_fitted(object, digits, ...)

    ## have special summary info printed for cv.glmnet objects
    if (class(object$model)[1] == "cv.glmnet")
    {
        est.coef <- predict(object$model, type = "coef", s = "lambda.min")
        vnames   <- rownames(est.coef)
        est.coef <- as.vector(est.coef)
        sel.idx  <- which(est.coef != 0)

        cat("\n")

        ## if variables are selected print out how many are selected
        ## and their coefficient estimates
        if (length(sel.idx) > 0)
        {
            sel.varnames <- vnames[sel.idx]
            cat(length(sel.idx), "variables selected by the lasso (cross validation criterion).\n\n")
            coefmat <- matrix(est.coef[sel.idx], ncol = 1)
            rownames(coefmat) <- sel.varnames
            colnames(coefmat) <- "Estimate"
            print.default(coefmat, quote = FALSE, right = TRUE, na.print = "NA", ...)
        } else
        {
            cat("No variables selected by the lasso with cross validation. \n\n")
        }
    } else
    {
        return(summary(object$model))
    }
}
