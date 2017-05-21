
#' @export
validate.subgrp <- function(model,
                            B              = 100L,
                            method         = c("boot_bias_correction",
                                               "training_test_replication"),
                            train.fraction = 0.5)
{
    method <- match.arg(method)

    B <- as.integer(B[1])
    if (B <= 1) stop("B must be a strictly positive integer")

    train.fraction <- as.numeric(train.fraction[1])

    if (train.fraction >= 1 | train.fraction <= 0)
    {
        stop("train.fraction must be between 0 and 1")
    }

    if (is.null(model$call)) stop("retcall argument must be set to TRUE for fitted model")


    model$call$retcall <- FALSE

    x   <- model$call$x
    trt <- model$call$trt
    y   <- model$call$y

    n.obs <- NROW(x)

    boot.list <- vector(mode = "list", length = B)
    for (b in 1:B)
    {
        samp.idx <- sample.int(n.obs, n.obs, replace = TRUE)
        model$call$x   <- x[samp.idx,]
        model$call$y   <- y[samp.idx]
        model$call$trt <- trt[samp.idx]

        mod.b    <- do.call(fit.subgrp, model$call)
        boot.list[[b]] <- mod.b$subgroup.trt.effects
    }

    summary.stats <- boot.list[[1]]

    for (b in 2:B)
    {
        for (l in 1:length(summary.stats))
        {
            summary.stats[[l]] <- summary.stats[[l]] + boot.list[[b]][[l]]
        }
    }

    for (l in 1:length(summary.stats))
    {
        summary.stats[[l]] <- summary.stats[[l]] / B
    }

    list(avg.results = summary.stats, boot.results = boot.list)
}
