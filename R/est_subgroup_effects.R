
#' Computes treatment effects within various subgroups
#'
#' @description Computes treatment effects within various subgroups to estimate subgroup treatment effects
#'
#' @param benefit.scores vector of estimated benefit scores
#' @param y The response vector
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param pi.x The propensity score for each observation
#' @param cutpoint numeric value for patients with benefit scores above which
#' (or below which if \code{larger.outcome.better = FALSE})
#' will be recommended to be in the treatment group. Can also set \code{cutpoint = "median"}, which will
#' use the median value of the benefit scores as the cutpoint or can set specific quantile values via \code{"quantx"}
#' where \code{"x"} is a number between 0 and 100 representing the quantile value; e.g. \code{cutpoint = "quant75"}
#' will use the 75th perent upper quantile of the benefit scores as the quantile.
#' @param larger.outcome.better boolean value of whether a larger outcome is better. Set to \code{TRUE}
#' if a larger outcome is better and set to \code{FALSE} if a smaller outcome is better. Defaults to \code{TRUE}.
#' @param reference.trt index of which treatment is the reference (in the case of multiple treatments).
#' This should be known already, as for a \code{trt} with K-levels, there will be K-1 benefit scores (1 per column)
#' of \code{benefit.scores}, where each column is a comparison of each K-1 treatments with the reference treatment.
#' The default is the last level of \code{trt} if it is a factor.
#' @seealso \code{\link[personalized]{fit.subgroup}} for function which fits subgroup identification models which generate
#' benefit scores.
#' @importFrom stats median weighted.mean
#' @export
subgroup.effects <- function(benefit.scores, y, trt,
                             pi.x,
                             cutpoint = 0,
                             larger.outcome.better = TRUE,
                             reference.trt = NULL)
{

    benefit.scores <- drop(benefit.scores)
    y   <- drop(y)

    family <- "standard"
    if (class(y) == "Surv")
    {
        family <- "cox"
    }

    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
        n.trts      <- length(unique.trts)
    } else
    {
        unique.trts <- sort(unique(trt))
        n.trts      <- length(unique.trts)
    }

    if (n.trts - 1 > NCOL(benefit.scores))
    {
        stop("number of unique treatments is larger than the number of treatments represented by the benefit scores provided")
    }

    if (NROW(benefit.scores) != NROW(y))     stop("length of benefit.scores and y do not match")
    if (NROW(benefit.scores) != length(trt)) stop("length of benefit.scores and trt do not match")

    cutpoint <- convert.cutpoint(cutpoint, benefit.scores)

    if (is.null(reference.trt) | n.trts == 2) # two trt options defaults to first as reference
    {
        reference.idx   <- 1L
        comparison.idx  <- (1:n.trts)[-reference.idx]
        comparison.trts <- unique.trts[-reference.idx]
        reference.trt   <- unique.trts[reference.idx]
    } else
    {
        reference.idx   <- which(unique.trts == reference.trt)
        comparison.idx  <- (1:n.trts)[-reference.idx]
        comparison.trts <- unique.trts[-reference.idx]
    }

    if (n.trts > 2)
    {
        # meaning of larger vs smaller benefit score
        # is different depending on whether larger means
        # better or not for the outcome
        if (larger.outcome.better)
        {
            best.comp.idx   <- apply(benefit.scores, 1, which.max)
            recommended.trt <- 1 * (benefit.scores > cutpoint)
            rec.ref         <- rowSums(recommended.trt) == 0

            recommended.trt <- ifelse(rec.ref, reference.trt, comparison.trts[best.comp.idx])
        } else
        {
            best.comp.idx   <- apply(benefit.scores, 1, which.min)
            recommended.trt <- 1 * (benefit.scores < cutpoint)
            rec.ref         <- rowSums(recommended.trt) == 0

            recommended.trt <- ifelse(rec.ref, reference.trt, comparison.trts[best.comp.idx])
        }
    } else
    {
        # meaning of larger vs smaller benefit score
        # is different depending on whether larger means
        # better or not for the outcome
        if (larger.outcome.better)
        {
            recommended.trt <- ifelse(benefit.scores > cutpoint, comparison.trts, reference.trt)
        } else
        {
            recommended.trt <- ifelse(benefit.scores < cutpoint, comparison.trts, reference.trt)
        }
    }

    wts <- create.weights(pi.x   = pi.x,
                          trt    = trt,
                          method = "weighting")

    ## old way before multiple trtments::
    ##
    ## compute mean of outcome within
    ## group of patients who both
    ## received and were recommended the treatment group
    #idx.11  <- (recommended.trt == 1) & (trt == 1)

    ## compute mean of outcome within
    ## group of patients who
    ## received control and were recommended the treatment group
    #idx.10  <- (recommended.trt == 1) & (trt == 0)

    ## compute mean of outcome within
    ## group of patients who
    ## received treatment and were recommended the control group
    #idx.01  <- (recommended.trt == 0) & (trt == 1)

    ## compute mean of outcome within
    ## group of patients who both
    ## received and were recommended the control group
    #idx.00  <- (recommended.trt == 0) & (trt == 0)


    idx.list <- rep(list(vector(mode = "list", length = n.trts)), n.trts)

    ## loop over all combinations of
    ## recommendation and observed trt status
    ## and find who has what combination
    ## of recommended trt and received trt
    for (t.recom in 1:n.trts)
    {
        for (t.receiv in 1:n.trts)
        {
            idx.list[[t.recom]][[t.receiv]] <- (recommended.trt == unique.trts[t.recom]) &
                                                           (trt == unique.trts[t.receiv])
        }
    }


    # mean.11 <- mean.10 <- mean.01 <- mean.00 <- numeric(1L)
    # if (family == "cox")
    # {
    #     survf.11 <- survfit(y[idx.11] ~ 1)
    #     mean.11  <- summary(survf.11)$table[5]
    #
    #     survf.10 <- survfit(y[idx.10] ~ 1)
    #     mean.10  <- summary(survf.10)$table[5]
    #
    #     survf.01 <- survfit(y[idx.01] ~ 1)
    #     mean.01  <- summary(survf.01)$table[5]
    #
    #     survf.00 <- survfit(y[idx.00] ~ 1)
    #     mean.00  <- summary(survf.00)$table[5]
    # } else
    # {
    #     mean.11 <- mean(y[idx.11])
    #     mean.10 <- mean(y[idx.10])
    #     mean.01 <- mean(y[idx.01])
    #     mean.00 <- mean(y[idx.00])
    # }

    res.mat <- matrix(0, ncol = n.trts, nrow = n.trts)
    colnames(res.mat) <- paste("Recommended", unique.trts)
    rownames(res.mat) <- paste("Received", unique.trts)

    #res.mat <- matrix(0, ncol = 2, nrow = 2)
    #colnames(res.mat) <- c("Recommended Trt", "Recommended Ctrl")
    #rownames(res.mat) <- c("Received Trt",    "Received Ctrl")

    sample.size.mat <- res.mat

    subgroup.effects <- numeric(n.trts)

    if (family == "cox")
    {
        for (t.recom in 1:n.trts)
        {
            for (t.receiv in 1:n.trts)
            {
                idx.cur <- idx.list[[t.recom]][[t.receiv]]

                if (sum(idx.cur))
                {
                    survf <- survfit(y[idx.cur] ~ 1, weights = wts[idx.cur])
                    restricted.mean <- summary(survf)$table[5]

                    res.mat[t.receiv, t.recom] <- restricted.mean
                    sample.size.mat[t.receiv, t.recom] <- sum(idx.cur)

                    if (t.recom == t.receiv)
                    {
                        #idx.disagree <- Reduce("|", idx.list[[t.recom]][-t.receiv])
                        idx.disagree <- (recommended.trt == unique.trts[t.recom]) &
                            (trt != unique.trts[t.recom])
                        if (sum(idx.disagree))
                        {
                            survf <- survfit(y[idx.disagree] ~ 1, weights = wts[idx.disagree])
                            restricted.mean <- summary(survf)$table[5]

                            subgroup.effects[t.recom] <- res.mat[t.receiv, t.recom] - restricted.mean
                        } else
                        {
                            subgroup.effects[t.recom] <- NaN
                        }
                    }
                } else
                {
                    res.mat[t.receiv, t.recom] <- NaN
                    subgroup.effects[t.recom]  <- NaN
                }

            }
        }
    } else
    {
        for (t.recom in 1:n.trts)
        {
            for (t.receiv in 1:n.trts)
            {
                idx.cur <- idx.list[[t.recom]][[t.receiv]]
                res.mat[t.receiv, t.recom] <- weighted.mean(y[idx.cur], w = wts[idx.cur])
                sample.size.mat[t.receiv, t.recom] <- sum(idx.cur)

                if (t.recom == t.receiv)
                {
                    #idx.disagree <- Reduce("|", idx.list[[t.recom]][-t.receiv])
                    idx.disagree <- (recommended.trt == unique.trts[t.recom]) &
                        (trt != unique.trts[t.recom])
                    subgroup.effects[t.recom] <- res.mat[t.receiv, t.recom] - weighted.mean(y[idx.disagree], w = wts[idx.disagree])
                }

            }
        }
    }

    idx.agree <- recommended.trt == trt

    if (family == "cox")
    {
        if (sum(idx.agree))
        {
            survf.agree <- survfit(y[idx.agree] ~ 1, weights = wts[idx.agree])
            restricted.mean.agree <- summary(survf.agree)$table[5]
        } else
        {
            restricted.mean.agree <- NaN
        }

        if (sum(!idx.agree))
        {
            survf.disagree <- survfit(y[!idx.agree] ~ 1, weights = wts[!idx.agree])
            restricted.mean.disagree <- summary(survf.disagree)$table[5]
        } else
        {
            restricted.mean.disagree <- NaN
        }

        overall.subgroup.effect <- restricted.mean.agree - restricted.mean.disagree
    } else
    {
        overall.subgroup.effect <- weighted.mean(y[idx.agree], w = wts[idx.agree]) -
            weighted.mean(y[!idx.agree], w = wts[!idx.agree])
    }

    #res.mat[1,1] <- mean.11
    #res.mat[1,2] <- mean.01
    #res.mat[2,1] <- mean.10
    #res.mat[2,2] <- mean.00

    #sample.size.mat[1,1] <- sum(idx.11)
    #sample.size.mat[1,2] <- sum(idx.01)
    #sample.size.mat[2,1] <- sum(idx.10)
    #sample.size.mat[2,2] <- sum(idx.00)

    #subgroup.effects <- c(mean.11 - mean.10,
    #                      mean.00 - mean.01)

    #names(subgroup.effects) <- c("Trt  Effect Among Recommended Trt",
    #                             "Ctrl Effect Among Recommended Ctrl")

    names(subgroup.effects) <- paste(unique.trts, "effect among recommended", unique.trts)

    list(subgroup.effects = subgroup.effects, # subgroup-specific effects
                                              # (trt effect among those recommended trt and
                                              #  ctrl effect among rec ctrl)
         avg.outcomes     = res.mat,          # means within 2x2 table (trt status vs trt rec)
         sample.sizes     = sample.size.mat,  # sample sizes for 2x2 table
         overall.subgroup.effect = overall.subgroup.effect) # overall difference in means among those
                                                            # whose recommendation agrees with what they received
                                                            # vs those whose recommendation differs from what they received
}
