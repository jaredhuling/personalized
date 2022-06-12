# Define common predictions function types
get.pred.func <- function(fit.name, model, env = parent.frame())
{
    n.trts    <- env$n.trts
    vnames    <- env$vnames
    sel.idx   <- env$sel.idx
    best.iter <- env$best.iter
    family    <- env$family
    # GAM models
    if (grepl("_gam$",fit.name))
    {
        if (grepl("_cox", fit.name))
        {
            pred.func <- function(x, type = c("link", "class"))
            {
                type <- match.arg(type)
                df.pred <- data.frame(cbind(1, x[,sel.idx[-1] - 1]))
                colnames(df.pred) <- vnames
                df.pred$trt_1n1 <- 1
                -drop(predict(model, newdata = df.pred, type = "link"))
            }
        } else
        {
            pred.func <- function(x, type = c("link", "class"))
            {
                type <- match.arg(type)
                df.pred <- data.frame(cbind(1, x[,sel.idx[-1] - 1]))
                colnames(df.pred) <- vnames
                df.pred$trt_1n1 <- 1
                drop(predict(model, newdata = df.pred, type = "link"))
            }
        }
        # xgboost models
    } else if (grepl("_xgboost$", fit.name))
    {
        pred.func <- function(x, type = c("link", "class"))
        {
            type <- match.arg(type)
            df <- xgb.DMatrix(x)

            if (type == "link")
            {
                outputmargin <- TRUE
            } else
            {
                outputmargin <- FALSE
            }

            drop(predict(model, newdata = df, outputmargin = outputmargin))
        }
        # non-GAM/xgboost LASSO models (loss ends in _lasso)
    } else if (grepl("_lasso$",fit.name))
    {
        if (grepl("_cox", fit.name))
        {
            pred.func <- function(x, type = c("link", "class"))
            {
                type <- match.arg(type)
                if (n.trts == 2)
                {
                    -drop(predict(model, newx = cbind(1, x),
                                  type = "link",
                                  s = "lambda.min",
                                  newoffset = rep(0, NROW(x)) ))
                } else
                {
                    ## need to handle cases with multiple treatments specially
                    ## because we don't want to sum up over all the estimated deltas.
                    ## for K-trtments we estimate K-1 delta functions and thus need
                    ## to extract each one individually.
                    all.coefs <- as.vector(predict(model, type = "coefficients", s = "lambda.min"))
                    n.coefs.per.trt <- length(all.coefs) / (n.trts - 1)

                    n.preds  <- NROW(x)
                    pred.mat <- array(NA, dim = c(n.preds, n.trts - 1))
                    for (t in 1:(n.trts - 1))
                    {
                        idx.coefs.cur <- (n.coefs.per.trt * (t - 1) + 1):(n.coefs.per.trt * t)
                        coefs.cur     <- all.coefs[idx.coefs.cur]

                        pred.mat[,t]  <- drop(cbind(1, x) %*% coefs.cur)
                    }
                    -pred.mat

                }
            }
        } else
        {
            pred.func <- function(x, type = c("link", "class"))
            {
                type <- match.arg(type)
                if (n.trts == 2)
                {
                    drop(predict(model, newx = cbind(1, x),
                                 type = "link",
                                 s = "lambda.min",
                                 newoffset = rep(0, NROW(x)) ))
                } else
                {
                    ## need to handle cases with multiple treatments specially
                    ## because we don't want to sum up over all the estimated deltas.

                    if (family == "multinomial")
                    {
                        drop(predict(model, cbind(1, x),
                                     type = type,
                                     s = "lambda.min"))
                    } else
                    {
                        ## for K-trtments we estimate K-1 delta functions and thus need
                        ## to extract each one individually.
                        all.coefs <- as.vector(predict(model, type = "coefficients", s = "lambda.min"))[-1]
                        n.coefs.per.trt <- length(all.coefs) / (n.trts - 1)

                        n.preds  <- NROW(x)
                        pred.mat <- array(NA, dim = c(n.preds, n.trts - 1))
                        for (t in 1:(n.trts - 1))
                        {
                            idx.coefs.cur <- (n.coefs.per.trt * (t - 1) + 1):(n.coefs.per.trt * t)
                            coefs.cur     <- all.coefs[idx.coefs.cur]

                            pred.mat[,t]  <- drop(cbind(1, x) %*% coefs.cur)
                        }
                        pred.mat
                    }
                }
            }
        }
    } else if (grepl("hinge_loss$", fit.name))
    {
        pred.func <- function(x, type = c("link", "class"))
        {
            drop(predict(model, newx = cbind(1, x), type = "linear.predictor"))
        }
    } else
    {
        stop(paste0("No prediction method found for loss: ", fit.name))
    }
    return(pred.func)
} # End get.pred.func

# Define common coefficient return methods
get.coef.func <- function(fit.name, env = parent.frame())
{
    n.trts <- env$n.trts
    # GAM or LASSO_GAM models (using cv.glmnet())
    if ( grepl("_lasso$", fit.name) )
    {
        coef.func <- function(mod)
        {
            coef(mod, s = "lambda.min")
        }
        # LOSS_GAM models (using gam() )
    } else if ( grepl("_loss_gam$",fit.name) & !grepl("lasso_gam$", fit.name))
    {
        coef.func <- function(mod)
        {
            coef(mod)
        }
    } else
    {
        coef.func <- function(mod)
        {
            return(NULL)
        }
    }
    return(coef.func)
} # End get.coef.func

#' @import glmnet
#' @importFrom stats coef
fit_sq_loss_lasso <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, intercept = FALSE, ...)
{
    # this function must return a fitted model
    # in addition to a function which takes in
    # a design matrix and outputs estimated benefit scores

    ###################################################################
    ##
    ## IMPORTANT NOTE: the name of this function *must*
    ##                 begin with "fit_" and end with
    ##                 the text string to associated with
    ##                 this function in the options for the
    ##                 'loss' argument of the fit.subgrp()
    ##                 function
    ##
    ###################################################################

    list.dots <- list(...)
    dot.names <- names(list.dots)

    n.unique.vars <- ncol(x) / (n.trts - 1)
    zero.pen.idx  <- ((1:(n.trts - 1) ) - 1) * n.unique.vars + 1

    list.dots$intercept <- intercept

    if ("penalty.factor" %in% dot.names)
    {
        ## ensure treatment is not penalized
        list.dots$penalty.factor[zero.pen.idx] <- 0
    } else
    {
        list.dots$penalty.factor <- rep(1, ncol(x))
        list.dots$penalty.factor[zero.pen.idx] <- 0
    }

    ## Establish nfolds for cv.glmnet()
    if ("nfolds" %in% dot.names)
    {
        nfolds <- list.dots$nfolds
        if (nfolds < 3)
        {
            stop("nfolds must be bigger than 3; nfolds = 10 recommended")
        }
    } else
    {
        nfolds <- 10
    }
    list.dots$nfolds <- nfolds

    nsel <- 0
    ct   <- 0
    ntry <- 4

    while(nsel == 0 & ct <= ntry)
    {
        ct <- ct + 1
        ## Establish foldid for cv.glmnet()
        ## if match.id was supplied, foldid will be structured around the clusters
        if (!is.null(match.id))
        {
            if ("foldid" %in% dot.names)
            {
                warning("User-supplied foldid will be ignored since match.id was detected.
                        Folds will be randomly assigned to clusters according to match.id.")
            }
            # Assign a fold ID for each cluster level
            df.folds <- data.frame(match.id = sample(levels(match.id)),
                                   fold.id = 1:length(levels(match.id)) %% nfolds)
            # Obtain vector of fold IDs with respect to the data
            foldid <- sapply(match.id, function(z) {df.folds[which(z == df.folds$match.id),"fold.id"]}) + 1
        } else
        {
            if ("foldid" %in% dot.names)
            {
                foldid <- list.dots$foldid
            } else
            {
                foldid <- sample(rep(seq(nfolds), length = nrow(x)))
            }
        }
        list.dots$foldid <- foldid

        # fit a model with a lasso
        # penalty and desired loss
        model <- do.call(cv.glmnet, c(list(x = x, y = y, weights = wts, family = family), list.dots))

        # this is needed for OWL losses, as glmnet
        # no longer allows constant columns (ie an intercept)
        # to have estimated coefficients
        if (intercept)
        {
            if (family != "multinomial")
            {
                model$glmnet.fit$beta[1,] <- unname(model$glmnet.fit$a0)
                model$glmnet.fit$a0       <- rep(0, length(model$glmnet.fit$a0))
            } else
            {
                for (cl in 1:nrow(model$glmnet.fit$a0))
                {
                    model$glmnet.fit$beta[[cl]][1,] <- unname(model$glmnet.fit$a0[cl,])
                    model$glmnet.fit$a0[cl,]        <- rep(0, length(model$glmnet.fit$a0[cl,]))
                }
            }
        }

        coefs <- get.coef.func("fit_sq_loss_lasso")(model)

        if (is.list(coefs))
        {
            nsel  <- sum(sapply(coefs, function(cfs) sum(cfs != 0))) - (n.trts - 1)
        } else
        {
            nsel  <- sum(coefs != 0) - (n.trts - 1)
        }
    }



    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_sq_loss_lasso", model),
         model        = model,
         coefficients = coefs)
}


fit_logistic_loss_lasso <- fit_sq_loss_lasso
fit_poisson_loss_lasso  <- fit_sq_loss_lasso

#' @import survival
fit_cox_loss_lasso <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, ...)
{

    list.dots <- list(...)
    dot.names <- names(list.dots)

    n.unique.vars <- ncol(x) / (n.trts - 1)
    zero.pen.idx  <- ((1:(n.trts - 1) ) - 1) * n.unique.vars + 1

    if ("penalty.factor" %in% dot.names)
    {
        ## ensure treatment is not penalized
        list.dots$penalty.factor[zero.pen.idx] <- 0
    } else
    {
        list.dots$penalty.factor <- rep(1, ncol(x))
        list.dots$penalty.factor[zero.pen.idx] <- 0
    }

    ## Establish nfolds for cv.glmnet()
    if ("nfolds" %in% dot.names)
    {
        nfolds <- list.dots$nfolds
        if (nfolds < 3)
        {
            stop("nfolds must be bigger than 3; nfolds = 10 recommended")
        }
    } else
    {
        nfolds <- 10
    }
    list.dots$nfolds <- nfolds

    nsel <- 0
    ct   <- 0
    ntry <- 4

    while(nsel == 0 & ct <= ntry)
    {
        ct <- ct + 1
        ## Establish foldid for cv.glmnet()
        ## if match.id was supplied, foldid will be structured around the clusters
        if (!is.null(match.id))
        {
            if ("foldid" %in% dot.names)
            {
                warning("User-supplied foldid will be ignored since match.id was detected.
                        Folds will be randomly assigned to clusters according to match.id.")
            }
            # Assign a fold ID for each cluster level
            df.folds <- data.frame(match.id = sample(levels(match.id)),
                                   fold.id = 1:length(levels(match.id)) %% nfolds)
            # Obtain vector of fold IDs with respect to the data
            foldid <- sapply(match.id, function(z) {df.folds[which(z == df.folds$match.id),"fold.id"]}) +1
        } else
        {
            if ("foldid" %in% dot.names)
            {
                foldid <- list.dots$foldid
            } else
            {
                foldid <- sample(rep(seq(nfolds), length = nrow(x)))
            }
        }
        list.dots$foldid <- foldid

        # fit a model with a lasso
        # penalty and desired loss
        model <- do.call(cv.glmnet, c(list(x = x, y = y, weights = wts, family = "cox"), list.dots))

        coefs <- get.coef.func("fit_cox_loss_lasso")(model)

        if (is.list(coefs))
        {
            nsel  <- sum(sapply(coefs, function(cfs) sum(cfs != 0))) - (n.trts - 1)
        } else
        {
            nsel  <- sum(coefs != 0) - (n.trts - 1)
        }
    }



    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_cox_loss_lasso", model),
         model        = model,
         coefficients = coefs)
}


#' @import mgcv
#' @importFrom stats as.formula binomial gaussian
fit_sq_loss_lasso_gam <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, intercept = FALSE, ...)
{
    # this function must return a fitted model
    # in addition to a function which takes in
    # a design matrix and outputs estimated benefit scores

    ###################################################################
    ##
    ## IMPORTANT NOTE: the name of this function *must*
    ##                 begin with "fit_" and end with
    ##                 the text string to associated with
    ##                 this function in the options for the
    ##                 'loss' argument of the fit.subgrp()
    ##                 function
    ##
    ###################################################################

    # need to inspect the dots to extract
    # the arguments supplied to cv.glmnet
    # and those supplied to gam
    list.dots <- list(...)
    dot.names <- names(list.dots)

    if ("penalty.factor" %in% dot.names)
    {
        ## ensure treatment is not penalized
        list.dots$penalty.factor[1] <- 0
    } else
    {
        list.dots$penalty.factor <- c(0, rep(1, ncol(x) - 1))
    }

    list.dots$intercept <- intercept


    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
    } else
    {
        unique.trts <- sort(unique(trt))
    }


    if (n.trts == 2)
    {
        trt.y <- trt

        ## can also use 'trt.multiplier'
        trt_1n1 <- ifelse(trt == unique.trts[2], 1, -1)
    } else
    {
        stop("gam loss not yet available for multiple treatments scenarios.")
    }

    ## Establish nfolds for cv.glmnet()
    if ("nfolds" %in% dot.names)
    {
        nfolds <- list.dots$nfolds
        if (nfolds < 3)
        {
            stop("nfolds must be bigger than 3; nfolds = 10 recommended")
        }
    } else
    {
        nfolds <- 10
    }
    list.dots$nfolds <- nfolds

    ## Establish foldid for cv.glmnet()
    ## if match.id was supplied, foldid will be structured around the clusters
    if (!is.null(match.id))
    {
        if ("foldid" %in% dot.names)
        {
            warning("User-supplied foldid will be ignored since match.id was detected.
                     Folds will be randomly assigned to clusters according to match.id.")
        }
        # Assign a fold ID for each cluster level
        df.folds <- data.frame(match.id = sample(levels(match.id)),
                               fold.id = 1:length(levels(match.id)) %% nfolds)
        # Obtain vector of fold IDs with respect to the data
        foldid <- sapply(match.id, function(z) {df.folds[which(z == df.folds$match.id),"fold.id"]}) + 1
    } else
    {
        if ("foldid" %in% dot.names)
        {
            foldid <- list.dots$foldid
        } else
        {
            foldid <- sample(rep(seq(nfolds), length = nrow(x)))
        }
    }
    list.dots$foldid <- foldid

    glmnet.argnames <- union(names(formals(cv.glmnet)), names(formals(glmnet)))
    gam.argnames    <- names(formals(gam))

    # since 'method' is an argument of 'fit.subgrp',
    # let the user change the gam 'method' arg by supplying
    # 'method.gam' arg instead of 'method'
    dot.names[dot.names == "method.gam"] <- "method"
    names(list.dots)[names(list.dots) == "method.gam"] <- "method"

    # find the arguments relevant for each
    # possible ...-supplied function
    dots.idx.glmnet <- match(glmnet.argnames, dot.names)
    dots.idx.gam    <- match(gam.argnames, dot.names)

    dots.idx.glmnet <- dots.idx.glmnet[!is.na(dots.idx.glmnet)]
    dots.idx.gam    <- dots.idx.gam[!is.na(dots.idx.gam)]

    # fit a model with a lasso
    # penalty and desired loss:
    sel.model <- do.call(cv.glmnet, c(list(x = trt.multiplier * x, y = y, weights = wts, family = family),
                                      list.dots[dots.idx.glmnet]))


    vnames <- colnames(x)

    sel.idx <- drop(predict(sel.model, type = "nonzero", s = "lambda.min")[[1]])

    # always include treatment main effect
    sel.idx <- union(1L, sel.idx)

    # names of selected variables
    sel.vnames <- vnames[sel.idx]

    # find which variables are binary
    var.levels <- numeric(length(sel.idx))
    for (v in 1:length(sel.idx))
    {
        var.levels[v] <- length(unique(x[,sel.idx[v]]))
    }

    contin.vars <- sel.vnames[var.levels > 2]
    binary.vars <- sel.vnames[var.levels <= 2]

    # create formula for gam
    contin.formula <- binary.formula <- NULL

    # don't create smoother for binary vars
    if (length(binary.vars) > 0)
    {
        binary.formula <- paste(binary.vars, collapse = "+")
    }

    # create smoother for each continuous var
    if (length(contin.vars) > 0)
    {
        num_unique_values <- apply(x[,contin.vars,drop=FALSE], 2, function(x) length(unique(x)) )

        form.cur <- paste0("s(", contin.vars, ", by = trt_1n1)")

        form.cur[num_unique_values <= 10] <- paste0("s(", contin.vars[num_unique_values <= 10], ", by = trt_1n1, k=",
                                                    num_unique_values[num_unique_values <= 10]-1, ")")

        contin.formula <- paste(form.cur, collapse = "+")
    }

    family.func <- gaussian()



    if (family == "cox")
    {
        rhs.formula <- paste(c(binary.formula, contin.formula), collapse = "+")
        family.func <- cox.ph()
    } else
    {
        rhs.formula <- paste("-1 +", paste(c(binary.formula, contin.formula), collapse = "+"))
        if (family == "binomial")
        {
            family.func <- binomial()
            y <- as.integer(y)
        } else if (family == "poisson")
        {
            family.func <- poisson()
            y <- as.integer(y)
        }
    }
    gam.formula <- as.formula(paste("y ~", rhs.formula))

    # create data frame
    df <- data.frame(y = y, x = x[,sel.idx], trt_1n1 = trt_1n1)


    colnames(df) <- c("y", sel.vnames)

    ## need binary vars to also be multiplied!!
    df[,binary.vars] <- trt.multiplier * df[,binary.vars]

    vnames <- sel.vnames


    oversmoothing_factor <- sqrt(ncol(x) / (length(contin.vars) + 1))

    # fit gam model:
    # only add in dots calls if they exist
    if (length(dots.idx.glmnet) > 0)
    {
        model <- do.call(gam, c(list(formula = gam.formula, data = df,
                                     weights = wts, family = family.func,
                                     gamma = oversmoothing_factor, ## oversmooth since we're in a post-selection scenario
                                     drop.intercept = TRUE),
                                list.dots[dots.idx.gam]))
    } else
    {
        model <- do.call(gam, list(formula = gam.formula, data = df,
                                   weights = wts, family = family.func,
                                   gamma = oversmoothing_factor, ## oversmooth since we're in a post-selection scenario
                                   drop.intercept = TRUE))
    }

    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_sq_loss_lasso_gam", model),
         model        = model,
         coefficients = get.coef.func("fit_sq_loss_lasso_gam")(model))
}

fit_logistic_loss_lasso_gam <- fit_sq_loss_lasso_gam
fit_cox_loss_lasso_gam      <- fit_sq_loss_lasso_gam
fit_poisson_loss_lasso_gam  <- fit_sq_loss_lasso_gam



fit_sq_loss_gam <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, ...)
{
    # this function must return a fitted model
    # in addition to a function which takes in
    # a design matrix and outputs estimated benefit scores

    ###################################################################
    ##
    ## IMPORTANT NOTE: the name of this function *must*
    ##                 begin with "fit_" and end with
    ##                 the text string to associated with
    ##                 this function in the options for the
    ##                 'loss' argument of the fit.subgrp()
    ##                 function
    ##
    ###################################################################


    list.dots <- list(...)

    # since 'method' is an argument of 'fit.subgrp',
    # let the user change the gam 'method' arg by supplying
    # 'method.gam' arg instead of 'method'
    names(list.dots)[names(list.dots) == "method.gam"] <- "method"

    vnames  <- colnames(x)
    sel.idx <- seq_len(ncol(x))

    # names of selected variables
    sel.vnames <- vnames[sel.idx]

    # if (sel.vnames[1] == "1")
    # {
    #     sel.vnames[1]  <- "Trt1"
    #     colnames(x)[1] <- sel.vnames[1]
    # }


    if (is.factor(trt))
    {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        unique.trts <- levels(trt)
    } else
    {
        unique.trts <- sort(unique(trt))
    }


    if (n.trts == 2)
    {
        trt.y <- trt

        trt_1n1 <- ifelse(trt == unique.trts[2], 1, -1)
    } else
    {
        stop("gam loss not yet available for multiple treatments scenarios.")
    }

    # find which variables are binary
    var.levels <- numeric(length(sel.idx))
    for (v in 1:length(sel.idx))
    {
        var.levels[v] <- length(unique(x[,sel.idx[v]]))
    }

    contin.vars <- sel.vnames[var.levels > 2]
    binary.vars <- sel.vnames[var.levels <= 2]

    # create formula for gam
    contin.formula <- binary.formula <- NULL

    # don't create smoother for binary vars
    if (length(binary.vars) > 0)
    {
        binary.formula <- paste(binary.vars, collapse = "+")
    }

    # create smoother for each continuous var
    if (length(contin.vars) > 0)
    {
        num_unique_values <- apply(x[,contin.vars,drop=FALSE], 2, function(x) length(unique(x)) )

        form.cur <- paste0("s(", contin.vars, ", by = trt_1n1)")

        form.cur[num_unique_values <= 10] <- paste0("s(", contin.vars[num_unique_values <= 10], ", by = trt_1n1, k=",
                                                    num_unique_values[num_unique_values <= 10]-1, ")")

        contin.formula <- paste(form.cur, collapse = "+")
    }

    family.func <- gaussian()

    if (family == "cox")
    {
        rhs.formula <- paste(c(binary.formula, contin.formula), collapse = "+")
        family.func <- cox.ph()
    } else
    {
        rhs.formula <- paste("-1 +", paste(c(binary.formula, contin.formula), collapse = "+"))
        if (family == "binomial")
        {
            family.func <- binomial()
            y <- as.integer(y)
        } else if (family == "poisson")
        {
            family.func <- poisson()
            y <- as.integer(y)
        }
    }
    gam.formula <- as.formula(paste("y ~", rhs.formula))

    # create data frame
    df <- data.frame(y = y, x = x[,sel.idx], trt_1n1 = trt_1n1)
    colnames(df) <- c("y", sel.vnames)

    ## need binary vars to also be multiplied!!
    df[,binary.vars] <- trt.multiplier * df[,binary.vars]

    vnames <- sel.vnames

    # fit gam model:
    # only add in dots calls if they exist
    if (length(list.dots) > 0)
    {
        model <- do.call(gam, c(list(formula = gam.formula, data = df,
                                     weights = wts, family = family.func,
                                     drop.intercept = TRUE),
                                list.dots))
    } else
    {
        model <- do.call(gam, list(formula = gam.formula, data = df,
                                   weights = wts, family = family.func,
                                   drop.intercept = TRUE))
    }


    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_sq_loss_gam", model),
         model        = model,
         coefficients = get.coef.func("fit_sq_loss_gam")(model))
}


fit_logistic_loss_gam <- fit_sq_loss_gam
fit_poisson_loss_gam  <- fit_sq_loss_gam
fit_cox_loss_gam      <- fit_sq_loss_gam




fit_owl_hinge_loss <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, ...)
{

    list.dots <- list(...)
    dot.names <- names(list.dots)

    ipop.argnames  <- names(formals(ipop))
    wksvm.argnames <- names(formals(weighted.ksvm))

    # find the arguments relevant for each
    # possible ...-supplied function
    dots.idx.wksvm <- match(wksvm.argnames, dot.names)
    dots.idx.ipop  <- match(ipop.argnames, dot.names)

    dots.idx.wksvm <- dots.idx.wksvm[!is.na(dots.idx.wksvm)]
    dots.idx.ipop  <- dots.idx.ipop[!is.na(dots.idx.ipop)]

    list.dots <- list.dots[c(dots.idx.wksvm, dots.idx.ipop)]
    dot.names <- dot.names[c(dots.idx.wksvm, dots.idx.ipop)]

    ## Establish nfolds for cv.glmnet()
    if ("nfolds" %in% dot.names)
    {
        nfolds <- list.dots$nfolds
        if (nfolds < 2)
        {
            stop("nfolds must be bigger than 2; nfolds = 10 recommended")
        }
    } else
    {
        nfolds <- 10
    }
    list.dots$nfolds <- nfolds

    ## Establish foldid for cv.glmnet()
    ## if match.id was supplied, foldid will be structured around the clusters
    if (!is.null(match.id))
    {
        if ("foldid" %in% dot.names)
        {
            warning("User-supplied foldid will be ignored since match.id was detected.
                    Folds will be randomly assigned to clusters according to match.id.")
        }
        # Assign a fold ID for each cluster level
        df.folds <- data.frame(match.id = sample(levels(match.id)),
                               fold.id = 1:length(levels(match.id)) %% nfolds)
        # Obtain vector of fold IDs with respect to the data
        foldid <- sapply(match.id, function(z) {df.folds[which(z == df.folds$match.id),"fold.id"]}) +1
    } else
    {
        if ("foldid" %in% dot.names)
        {
            foldid <- list.dots$foldid
        } else
        {
            foldid <- sample(rep(seq(nfolds), length = nrow(x)))
        }
    }
    list.dots$foldid <- foldid

    # fit a model with a lasso
    # penalty and desired loss
    model <- do.call(weighted.ksvm, c(list(x = x, y = as.character(y), weights = wts), list.dots))

    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_hinge_loss", model),
         model        = model,
         coefficients = get.coef.func("fit_hinge_loss")(model))

}
