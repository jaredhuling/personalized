
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
        sel.varnames <- vnames[sel.idx]
        cat(length(sel.idx) - 1,
            "out of",
            length(est.coef) - 1,
            "variables selected in total by the lasso (cross validation criterion).\n\n")

        if (object$n.trts == 2)
        {
            coefmat <- matrix(est.coef[sel.idx], ncol = 1)
            rownames(coefmat) <- sel.varnames
            colnames(coefmat) <- "Estimate"

            print.default(round(coefmat, digits), quote = FALSE, right = TRUE, na.print = "NA", ...)
        } else
        {
            est.coef  <- predict(object$model, type = "coef", s = "lambda.min")[-1,,drop=FALSE]
            vnames    <- rownames(est.coef)

            ## remove the unecessary ".#" artificially added to variable names
            vnames    <- gsub("[.][0-9]*$", "", vnames)

            all.coefs <- unname(as.vector(drop(all.coefs)))

            n.coefs.per.trt <- length(all.coefs) / (object$n.trts - 1)
            for (t in 1:(object$n.trts - 1))
            {
                idx.coefs.cur <- (n.coefs.per.trt * (t - 1) + 1):(n.coefs.per.trt * t)
                coefs.cur     <- all.coefs[idx.coefs.cur]

                sel.idx  <- which(coefs.cur != 0)

                cat("\n")

                ## if variables are selected print out how many are selected
                ## and their coefficient estimates
                sel.varnames.cur <- vnames[sel.idx]
                cat(length(sel.idx) - 1,
                    "out of",
                    length(coefs.cur),
                    "variables selected for delta", t, "by the lasso (cross validation criterion).\n\n")

                coefmat <- matrix(coefs.cur[sel.idx], ncol = 1)

                rownames(coefmat) <- sel.varnames.cur
                colnames(coefmat) <- paste0("Estimates for delta(",
                                            object$comparison.trts[t], " vs ", object$reference.trt, ")" )

                print.default(round(coefmat, digits), quote = FALSE, right = TRUE, na.print = "NA", ...)
            }
        }

    } else
    {
        return(summary(object$model))
    }
}
