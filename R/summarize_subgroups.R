
#' Summarizing covariates within estimated subgroups
#'
#' @description Summarizes covariate values within the estimated subgroups
#'
#' @param object a fitted object from \code{fit.subgroup()}
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @importFrom stats t.test chisq.test
#' @export
summarize.subgroups <- function(object)
{

    if (is.null(object$call)) stop("retcall argument must be set to TRUE for fitted model
                                  to use summarize.subgroups()")


    # save data objects because they
    # will be written over by resampled versions later
    x      <- object$call$x
    bene.scores <- object$benefit.scores
    if (object$larger.outcome.better)
    {
        subgroup <- 1 * (bene.scores > 0)
    } else
    {
        subgroup <- 1 * (bene.scores < 0)
    }
    y      <- object$call$y

    n.obs  <- NROW(x)
    n.vars <- NCOL(x)
    vnames <- object$var.names

    # find which variables are binary
    var.levels <- numeric(n.vars)
    for (v in 1:n.vars)
    {
        var.levels[v] <- length(unique(x[,v]))
    }

    contin.vars <- vnames[var.levels > 2]
    binary.vars <- vnames[var.levels == 2]

    compare.mat <- array(0, dim = c(n.vars, 6L))

    for (v in 1:n.vars)
    {
        ## means within each subgroup
        compare.mat[,1] <- colMeans(x[subgroup == 1, ])
        compare.mat[,2] <- colMeans(x[subgroup == 0, ])
        compare.mat[,3] <- compare.mat[,1] - compare.mat[,2]



        if (var.levels[v] > 2)
        {
            ## standard errors within each subgroup
            compare.mat[v,5] <- sd(x[subgroup == 1, v]) / sqrt(sum(subgroup == 1))
            compare.mat[v,6] <- sd(x[subgroup == 0, v]) / sqrt(sum(subgroup == 0))

            ## run t.test for contin vars
            tt <- t.test(x[subgroup == 1, v], x[subgroup == 0, v])
            compare.mat[v,4] <- tt$p.value
        } else
        {
            ## run chi squared test for binary vars
            ct <- chisq.test(subgroup, x[, v])
            compare.mat[v,4] <- ct$p.value
        }
    }
    rownames(compare.mat) <- vnames

    compare.mat <- data.frame(compare.mat)
    colnames(compare.mat) <- c("mean (recomm. trt)", "mean (recomm. ctrl)", "diff",
                               "p.value", "SE (recomm. trt)", "SE (recomm. ctrl)")
    class(compare.mat) <- c("subgroup_summary", "data.frame")
    compare.mat
}

#' Printing summary results for fitted subgroup identification models
#'
#' @description Prints summary results for estimated subgroup treatment effects
#'
#' @param p.value a p-value threshold for mean differences below which covariates will be displayed. For example,
#' setting \code{p.value = 0.05} will display all covariates that have a significant difference between subgroups
#'  with p-value less than 0.05. Defaults to 1, which displays all covariates
#' @seealso \code{\link[personalized]{summarize.subgroups}} for function which summarizes subgroup covariate values
#' @rdname print
#' @export
print.subgroup_summary <- function(x, p.value = 1, digits = max(getOption('digits')-3, 3), ...)
{
    compare.mat <- x[x$p.value <= p.value,]
    print.data.frame(compare.mat, digits = digits, quote = FALSE, right = TRUE, na.print = "NA", ...)
}

