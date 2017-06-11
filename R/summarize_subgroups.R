#' Summarizing covariates within estimated subgroups
#'
#' @description Summarizes covariate values within the estimated subgroups
#' @param x a fitted object from \code{fit.subgroup()} or a matrix of covariate values
#' @param ... optional arguments to \code{summarize.subgroups} methods
#' @importFrom stats t.test chisq.test
#' @export
summarize.subgroups <- function(x, ...) UseMethod("summarize.subgroups")


#' @param subgroup vector of indicators of same length as the number of rows in x if x is a matrix.
#' A value of 1 in the ith position of \code{subgroup} indicates patient i is in the subgroup
#' of patients recommended the treatment and a value of 0 in the ith position of \code{subgroup} indicates patient i is in the subgroup
#' of patients recommended the control.
#' If x is a fitted object returned by \code{fit.subgroup()}, \code{subgroup} is not needed.
#' @rdname summarize.subgroups
#' @export
summarize.subgroups.default <- function(x, subgroup, ...)
{
    vnames <- colnames(x)

    n.obs  <- NROW(x)
    n.vars <- NCOL(x)

    if (is.null(vnames))
    {
        vnames <- paste0("V", 1:n.vars)
    }


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
    colnames(compare.mat) <- c("avg (recom trt)", "avg (recom ctrl)", "diff",
                               "p.value", "SE (recom trt)", "SE (recom ctrl)")
    class(compare.mat) <- c("subgroup_summary", "data.frame")
    compare.mat
}


#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models.
#' @rdname summarize.subgroups
#' @export
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 1000
#' n.vars <- 50
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
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
#'                              trt = trt01,
#'                              propensity.func = prop.func,
#'                              loss   = "sq_loss_lasso",
#'                              nfolds = 5)    # option for cv.glmnet
#'
#' comp <- summarize.subgroups(subgrp.model)
#' print(comp, p.value = 0.01)
#'
#' # or we can simply supply the matrix x and the subgroups
#' comp2 <- summarize.subgroups(x, subgroup = 1 * (subgrp.model$benefit.scores > 0))
summarize.subgroups.subgroup_fitted <- function(x, ...)
{

    object <- x
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

    vnames <- object$var.names

    colnames(x) <- vnames

    summarize.subgroups.default(x = x, subgroup = subgroup)
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

