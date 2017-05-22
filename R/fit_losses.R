
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
    ##                 'loss' argument of the fit.subgrp() function
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
