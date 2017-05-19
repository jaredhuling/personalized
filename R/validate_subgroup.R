
validate.subgrp <- function(model,
                            B              = 100L,
                            method         = c("boot_bias_correction",
                                               "training_test_replication"),
                            train.fraction = 0.5)
{
    method <- match.arg(method)

    B <- as.integer(B[1])
    if (B < 1) stop("B must be a strictly positive integer")

    train.fraction <- as.numeric(train.fraction[1])

    if (train.fraction >= 1 | train.fraction <= 0)
    {
        stop("train.fraction must be between 0 and 1")
    }
}
