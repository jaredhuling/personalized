
#' @import glmnet
fit_sq_loss_lasso <- function(x, y, wts, family, ...)
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

    # fit a model with a lasso
    # penalty and desired loss
    model <- cv.glmnet(x = x,  y = y,
                       weights   = wts,
                       family    = family,
                       intercept = FALSE, ...)

    # define a function which inputs a design matrix
    # and outputs estimated benefit scores: one score
    # for each row in the design matrix
    pred.func <- function(x)
    {
        drop(predict(model, newx = cbind(1, x), type = "link", s = "lambda.min"))
    }

    list(predict = pred.func,
         model   = model)
}


fit_logistic_loss_lasso <- fit_sq_loss_lasso

fit_cox_loss_lasso <- function(x, y, wts, family, ...)
{
    model <- cv.glmnet(x = x,  y = y,
                       weights   = wts,
                       family    = "cox",
                       ...)

    pred.func <- function(x)
    {
        drop(predict(model, newx = cbind(1, x), type = "link", s = "lambda.min"))
    }

    list(predict = pred.func,
         model   = model)
}



#' @import mgcv
fit_sq_loss_lasso_gam <- function(x, y, wts, family, ...)
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

    # fit a model with a lasso
    # penalty and desired loss
    sel.model <- cv.glmnet(x = x,  y = y,
                           weights   = wts,
                           family    = family,
                           intercept = FALSE, ...)

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
    binary.vars <- sel.vnames[var.levels == 2]

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
        form.cur <- paste0("s(", contin.vars, ")")
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
        }
    }
    gam.formula <- as.formula(paste("y ~", rhs.formula))

    # create data frame
    df <- data.frame(y = y, x = x[,sel.idx])
    colnames(df) <- c("y", sel.vnames)

    # fit gam model:

    model <- gam(gam.formula, data = df, family = family.func, weights = wts)

    # define a function which inputs a design matrix
    # and outputs estimated benefit scores: one score
    # for each row in the design matrix
    pred.func <- function(x)
    {
        df.pred <- data.frame(cbind(1, x[,sel.idx[-1] - 1]))
        colnames(df.pred) <- colnames(df)[-1] # take out 'y' column name
        drop(predict(model, newdata = df.pred, type = "link"))
    }

    list(predict = pred.func,
         model   = model)
}

fit_logistic_loss_lasso_gam <- fit_sq_loss_lasso_gam
fit_cox_loss_lasso_gam      <- fit_sq_loss_lasso_gam

