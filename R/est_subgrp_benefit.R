
#' Computes means within various subgroups
#'
#' @description Computes means within various subgroups to estimate subgroup treatment effects
#'
#' @param benefit.scores vector of estimated benefit scores
#' @param y The response vector
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group
#' @param larger.outcome.better boolean value of whether a larger outcome is better. Set to \code{TRUE}
#' if a larger outcome is better and set to \code{FALSE} if a smaller outcome is better. Defaults to \code{TRUE}.
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models which generate
#' benefit scores.
#' @export
subgrp.benefit <- function(benefit.scores, y, trt, cutpoint = 0, larger.outcome.better = TRUE)
{

    benefit.scores <- drop(benefit.scores)
    y   <- drop(y)
    trt <- drop(trt)

    if (length(benefit.scores) != length(y))   stop("length of benefit.scores and y do not match")
    if (length(benefit.scores) != length(trt)) stop("length of benefit.scores and trt do not match")

    cutpoint <- as.numeric(cutpoint[1])

    # meaning of larger vs smaller benefit score
    # is different depending on whether larger means
    # better or not for the outcome
    if (larger.outcome.better)
    {
        recommended.trt <- 1 * (benefit.scores > cutpoint)
    } else
    {
        recommended.trt <- 1 * (benefit.scores < cutpoint)
    }

    # compute mean of outcome within
    # group of patients who both
    # received and were recommended the treatment group
    idx.11  <- (recommended.trt == 1) & (trt == 1)
    mean.11 <- mean(y[idx.11])

    # compute mean of outcome within
    # group of patients who
    # received control and were recommended the treatment group
    idx.10  <- (recommended.trt == 1) & (trt == 0)
    mean.10 <- mean(y[idx.10])

    # compute mean of outcome within
    # group of patients who
    # received treatment and were recommended the control group
    idx.01  <- (recommended.trt == 0) & (trt == 1)
    mean.01 <- mean(y[idx.01])

    # compute mean of outcome within
    # group of patients who both
    # received and were recommended the control group
    idx.00  <- (recommended.trt == 0) & (trt == 0)
    mean.00 <- mean(y[idx.00])

    res.mat <- matrix(0, ncol = 2, nrow = 2)
    colnames(res.mat) <- c("Recommended Trt", "Recommended Ctrl")
    rownames(res.mat) <- c("Received Trt",    "Received Ctrl")

    sample.size.mat <- res.mat

    res.mat[1,1] <- mean.11
    res.mat[1,2] <- mean.01
    res.mat[2,1] <- mean.10
    res.mat[2,2] <- mean.00

    sample.size.mat[1,1] <- sum(idx.11)
    sample.size.mat[1,2] <- sum(idx.01)
    sample.size.mat[2,1] <- sum(idx.10)
    sample.size.mat[2,2] <- sum(idx.00)

    subgroup.effects <- c(mean.11 - mean.10,
                          mean.00 - mean.01)

    names(subgroup.effects) <- c("Trt  Effect Among Recommended Trt",
                                 "Ctrl Effect Among Recommended Ctrl")

    list(subgroup.effects = subgroup.effects,
         avg.outcomes     = res.mat,
         sample.sizes     = sample.size.mat)
}
