
#' @import glmnet
fit_sq_loss_lasso <- function(x, y, wts, family, ...)
{
    model <- cv.glmnet(x = x,  y = y,
                       weights   = wts,
                       family    = family,
                       intercept = FALSE, ...)

    pred.func <- function(x)
    {
        drop(predict(model, newx = cbind(1, x), type = "link", s = "lambda.min"))
    }

    list(predict = pred.func,
         model   = model)
}


fit_logistic_loss_lasso <- fit_sq_loss_lasso
fit_cox_loss_lasso      <- fit_sq_loss_lasso

