

fit_sq_loss_lasso <- function(x, y, wts, family, ...)
{
    model <- cv.glmnet(x = x, y = y, weights = wts, family = family)

    pred.func <- function(x)
    {
        drop(predict(model, newx = x, type = "response", s = "lambda.min"))
    }

    list(predict = pred.func,
         model   = model)
}
