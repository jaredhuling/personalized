


create.augmentation.function <- function(family, crossfit = TRUE, nfolds.crossfit = 10, cv.glmnet.args = NULL)
{
    if (family == "binomial")
    {
        tm <- "auc"
    } else
    {
        tm <- "mse"
    }

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    cv.glmnet.args
}





glmnet_aug_kfold_crossfit <- function(x, y, trt, wts = NULL, K = 10,
                                      predtype = c("link", "response"),
                                      family = c("gaussian", "binomial", "poisson", "multinomial", "cox", "mgaussian"),
                                      interactions = TRUE, cv.glmnet.args = NULL)
{

    predtype <- match.arg(predtype)
    family   <- match.arg(family)


    if (family == "binomial")
    {
        tm <- "auc"
    } else
    {
        tm <- "mse"
    }

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "weights", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    if (is.null(wts))
    {
        wts <- rep(1, NROW(x))
    }

    if (interactions)
    {
        ## full model for nonzeroness
        df_all <- data.frame(x, trt = trt)
        df_1   <- data.frame(x, trt = 1)
        df_0   <- data.frame(x, trt = -1)

        mm_all <- model.matrix(~x*trt-1, data = df_all)
        mm_1   <- model.matrix(~x*trt-1, data = df_1)
        mm_0   <- model.matrix(~x*trt-1, data = df_0)
    } else
    {
        mm_all <- x
    }

    n <- NROW(mm_all)

    predvec <- numeric(n)

    foldid = sample(rep(seq(K), length = n))

    for (i in seq(K))
    {
        which <- foldid == i

        ## full model for nonzeroness
        # glmfit_zero_main  <- cv.glmnet(y            = y[!which],
        #                                x            = mm_all[!which,,drop=FALSE],
        #                                weights      = wts[!which],
        #                                family       = family,
        #                                parallel     = FALSE,
        #                                type.measure = type.measure)

        glmfit_zero_main <- do.call(cv.glmnet, c(list(y = y[!which], x = mm_all[!which,,drop=FALSE],
                                                      weights = wts[!which], family = family), cv.glmnet.args))

        if (interactions)
        {
            ## get predictions for trt = 1 & -1
            pred1_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_1[which,,drop=FALSE], s = "lambda.min", type = predtype)))
            pred0_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_0[which,,drop=FALSE], s = "lambda.min", type = predtype)))

            predvec[which] <- 0.5 * (pred1_zerr + pred0_zerr)
        } else
        {
            ## get predictions for trt = 1 & -1
            pred_zerr <- unname(drop(predict(glmfit_zero_main, newx = mm_all[which,,drop=FALSE], s = "lambda.min", type = predtype)))

            predvec[which] <- pred_zerr
        }

    }

    predvec
}



glmnet_propensity_kfold_crossfit <- function(x, trt, K = 10, cv.glmnet.args = NULL)
{

    n <- NROW(x)

    tm <- "auc"

    if (is.null(cv.glmnet.args))
    {

        cv.glmnet.args <- list(type.measure = tm, nfolds = 10)
    }

    cv.glmnet.args[c("x", "y", "family", "parallel")] <- NULL
    cv.glmnet.args$parallel <- FALSE


    if (!("type.measure" %in% names(cv.glmnet.args) ))
    {
        cv.glmnet.args$type.measure <- tm
    }

    propensvec <- numeric(n)

    foldid = sample(rep(seq(K), length = n))

    for (i in seq(K))
    {
        which <- foldid == i

        ## propensity score model fit on K-1 folds
        # glmfit_propens  <- cv.glmnet(y = trt[!which], x = x[!which,,drop=FALSE],
        #                              family = "binomial", #parallel = TRUE,
        #                              type.measure = type.measure)

        glmfit_propens <- do.call(cv.glmnet, c(list(y = trt[!which], x = x[!which,,drop=FALSE],
                                                    family = "binomial"), cv.glmnet.args))

        ## get propensity scores for the held out fold
        propensvec[which] <- unname(drop(predict(glmfit_propens, newx = x[which,,drop=FALSE],
                                                 s = "lambda.1se", type = "response")))


    }

    ## propensity scores will never be outside of 0 or 1 and
    ## shouldn't have missing values, but this code is a safety
    ## check just in case
    propensvec[is.na(propensvec)] <- mean(propensvec[!is.na(propensvec)])
    propensvec[propensvec <= 0] <- 1e-5
    propensvec[propensvec >= 1] <- 1 - 1e-5

    propensvec
}
