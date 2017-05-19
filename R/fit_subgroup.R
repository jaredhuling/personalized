
fit.subgrp <- function(x,
                       y,
                       trt,
                       pi.x,
                       family     = c("gaussian", "binomial", "cox"),
                       loss       = c("sq_loss_lasso", "logistic_loss_lasso", "cox_loss_lasso"),
                       method     = c("weighting", "a_learning"),
                       nfolds     = 10L,
                       foldid     = NULL,
                       ...)
{

    family <- match.arg(family)
    loss   <- match.arg(loss)
    method <- match.arg(method)

    this.call <- match.call()

    unique.trts <- sort(unique(trt))
    if (length(unique.trts) != 2) stop("trt must have 2 distinct levels")
    if (any(unique.trts) != c(0, 1)) stop("trt should be coded as 0 and 1")

    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("pi.x should be between 0 and 1")

    trt2 <- 2 * trt - 1

    if (method == "weighting")
    {
        x.tilde <- trt2 * cbind(1, x)
        wts     <- pi.x * (trt == 1) + (1 - pi.x) * (trt == 0)
    } else
    {
        x.tilde <- (trt - pi.x) * cbind(1, x)
        wts     <- rep(1, nrow(x))
    }

    fit_fun      <- paste0("fit_", loss)
    fitted.model <- do.call(fit_fun, list(x.tilde, y, wts, family, ...))

    fitted.model$loss <- loss
    fitted.model$method <- method
}
