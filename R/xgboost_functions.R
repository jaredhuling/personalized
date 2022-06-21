

return_squared_error_loss_xgboost <- function(trt_multiplier)
{
    weighted_squared_error_loss <- function(preds, dtrain)
    {
        labels <- getinfo(dtrain, "label")
        grad <- trt_multiplier * (trt_multiplier * preds - labels)
        hess <- trt_multiplier * trt_multiplier
        return(list(grad = grad, hess = hess))
    }
    weighted_squared_error_loss
}


return_eval_metric_xgboost <- function(trt_multiplier)
{
    evalerror <- function(preds, dtrain) {
        labels <- getinfo(dtrain, "label")
        weightss <- getinfo(dtrain, "weight")
        preds <- trt_multiplier * preds
        mse <- sqrt(weighted.mean((preds - labels) ^ 2, w = weightss))
        return(list(metric = "wtd_rmse", value = mse))
    }

    evalerror
}


## extract unexported xgboost functions
categorize.callbacks <- utils::getFromNamespace("categorize.callbacks", "xgboost")
check.custom.obj <- utils::getFromNamespace("check.custom.obj", "xgboost")
check.custom.eval <- utils::getFromNamespace("check.custom.eval", "xgboost")
check.deprecation <- utils::getFromNamespace("check.deprecation", "xgboost")
check.booster.params <- utils::getFromNamespace("check.booster.params", "xgboost")
generate.cv.folds <- utils::getFromNamespace("generate.cv.folds", "xgboost")
has.callbacks <- utils::getFromNamespace("has.callbacks", "xgboost")
add.cb <- utils::getFromNamespace("add.cb", "xgboost")
cb.print.evaluation <- utils::getFromNamespace("cb.print.evaluation", "xgboost")
cb.evaluation.log <- utils::getFromNamespace("cb.evaluation.log", "xgboost")
cb.early.stop <- utils::getFromNamespace("cb.early.stop", "xgboost")
cb.cv.predict <- utils::getFromNamespace("cb.cv.predict", "xgboost")
xgb.get.DMatrix <- utils::getFromNamespace("xgb.get.DMatrix", "xgboost")
slicexgb <- utils::getFromNamespace("slice", "xgboost")
xgb.Booster.handle <- utils::getFromNamespace("xgb.Booster.handle", "xgboost")
NVL <- utils::getFromNamespace("NVL", "xgboost")
xgb.iter.update <- utils::getFromNamespace("xgb.iter.update", "xgboost")
xgb.iter.eval <- utils::getFromNamespace("xgb.iter.eval", "xgboost")

#' @importFrom xgboost xgb.train
xgb_cv_personalized <- function (params = list(), data, trt.multiplier, nrounds, nfold, label = NULL,
                                 missing = NA, prediction = FALSE, showsd = TRUE, metrics = list(),
                                 obj_func = NULL, feval_func = NULL, stratified = TRUE, folds = NULL,
                                 train_folds = NULL, verbose = TRUE, print_every_n = 1L, early_stopping_rounds = NULL,
                                 maximize = NULL, callbacks = list(), ...)
{
    check.deprecation(...)
    params <- check.booster.params(params, ...)
    for (m in metrics) params <- c(params, list(eval_metric = m))

    obj <- obj_func(trt.multiplier)
    feval <- feval_func(trt.multiplier)
    check.custom.obj()
    check.custom.eval()
    if ((inherits(data, "xgb.DMatrix") && is.null(getinfo(data,
                                                          "label"))) || (!inherits(data, "xgb.DMatrix") && is.null(label))) {
        stop("Labels must be provided for CV either through xgb.DMatrix, or through 'label=' when 'data' is matrix")
    }
    else if (inherits(data, "xgb.DMatrix")) {
        if (!is.null(label))
            warning("xgb.cv: label will be ignored, since data is of type xgb.DMatrix")
        cv_label <- getinfo(data, "label")
    }
    else {
        cv_label <- label
    }
    if (!is.null(folds)) {
        if (!is.list(folds) || length(folds) < 2)
            stop("'folds' must be a list with 2 or more elements that are vectors of indices for each CV-fold")
        nfold <- length(folds)
    }
    else {
        if (nfold <= 1)
            stop("'nfold' must be > 1")
        folds <- generate.cv.folds(nfold, nrow(data), stratified,
                                             cv_label, params)
    }
    params <- c(params, list(silent = 1))
    print_every_n <- max(as.integer(print_every_n), 1L)
    if (!has.callbacks(callbacks, "cb.print.evaluation") && verbose) {
        callbacks <- add.cb(callbacks, cb.print.evaluation(print_every_n,
                                                           showsd = showsd))
    }
    evaluation_log <- list()
    if (!has.callbacks(callbacks, "cb.evaluation.log")) {
        callbacks <- add.cb(callbacks, cb.evaluation.log())
    }
    stop_condition <- FALSE
    if (!is.null(early_stopping_rounds) && !has.callbacks(callbacks,
                                                          "cb.early.stop")) {
        callbacks <- add.cb(callbacks, cb.early.stop(early_stopping_rounds,
                                                     maximize = maximize, verbose = verbose))
    }
    if (prediction && !has.callbacks(callbacks, "cb.cv.predict")) {
        callbacks <- add.cb(callbacks, cb.cv.predict(save_models = FALSE))
    }
    cb <- categorize.callbacks(callbacks)
    dall <- xgb.get.DMatrix(data, label, missing)
    bst_folds <- lapply(seq_along(folds), function(k) {
        dtest <- slicexgb(dall, folds[[k]])
        if (is.null(train_folds))
            dtrain <- slicexgb(dall, unlist(folds[-k]))
        else dtrain <- slicexgb(dall, train_folds[[k]])
        handle <- xgb.Booster.handle(params, list(dtrain, dtest))
        list(dtrain = dtrain, bst = handle, watchlist = list(test = dtest), index = folds[[k]],
             train_index = unlist(folds[-k]))
    })
    rm(dall)
    basket <- list()
    num_class <- max(as.numeric(NVL(params[["num_class"]], 1)),
                     1)
    num_parallel_tree <- max(as.numeric(NVL(params[["num_parallel_tree"]],
                                                      1)), 1)
    begin_iteration <- 1
    end_iteration <- nrounds
    for (iteration in begin_iteration:end_iteration) {
        for (f in cb$pre_iter) f()
        msg <- lapply(bst_folds, function(fd) {
            train_idx <- fd$train_index
            test_idx <- fd$index
            obj_current <- obj_func(trt.multiplier[train_idx])
            xgb.iter.update(fd$bst, fd$dtrain, iteration - 1,
                                      obj_current)
            feval_current <- feval_func(trt.multiplier[test_idx])
            xgb.iter.eval(fd$bst, fd$watchlist, iteration - 1,
                                    feval_current)
        })

        if (getRversion() >= "4.2.0")
        {
            msg <- simplify2array(msg, except = 0L)
        } else
        {
            msg <- simplify2array(msg)
            if (is.null(dim(msg)))
            {
                msg_names <- names(msg)
                msg <- matrix(msg, nrow = 1)
                rownames(msg) <- msg_names[1]
            }
        }


        bst_evaluation <- rowMeans(msg)
        bst_evaluation_err <- sqrt(rowMeans(msg^2) - bst_evaluation^2)
        for (f in cb$post_iter) f()
        if (stop_condition)
            break
    }
    for (f in cb$finalize) f(finalize = TRUE)
    ret <- list(call = match.call(), params = params, callbacks = callbacks,
                evaluation_log = evaluation_log, niter = end_iteration,
                nfeatures = ncol(data), folds = folds)
    ret <- c(ret, basket)
    class(ret) <- "xgb.cv.synchronous"
    invisible(ret)
}





#' @import xgboost
#' @importFrom xgboost slice
fit_sq_loss_xgboost <- function(x, y, trt, n.trts, wts, family, match.id, trt.multiplier, ...)
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
        stop("xgboost loss not yet available for multiple treatments scenarios.")
    }

    list.dots <- list(...)
    dot.names <- names(list.dots)

    if ("params" %in% dot.names)
    {
        params <- list.dots$params

        params$objective   <- NULL
        params$eval_metric <- NULL

    } else
    {
        params <- list(max_depth = 3,
                       eta       = 0.05,
                       nthread   = 2,
                       verbose   = 1,
                       subsample = 0.623,
                       colsample_bytree = 0.8,
                       booster   = "gbtree")
    }

    if ("verbose" %in% names(params))
    {
        verbose <- params$verbose
        params$verbose <- NULL
        if (!("verbose" %in% names(list.dots)))
        {
            list.dots$verbose <- verbose
        }
    }

    list.dots$params <- NULL

    if (!("early_stopping_rounds" %in% dot.names))
    {
        list.dots$early_stopping_rounds <- 50L
        warning("'early_stopping_rounds' not set; defaulting to 50")
    }

    if (!("nrounds" %in% dot.names))
    {
        list.dots$nrounds <- 500L
        warning("'nrounds' not set; defaulting to 500")
    }

    if (!("nfold" %in% dot.names))
    {
        list.dots$nfold <- 5L
        warning("'nfold' not set; defaulting to 5")
    }

    if (!is.null(match.id)) {
        warning("Matched groups are not guaranteed to remain matched in the cross-validation procedure using xgboost models.")
    }

    if ("offset" %in% dot.names)
    {
        dtrain <- xgb.DMatrix(x[,-1],
                              label = y,
                              weight = wts,
                              base_margin = list.dots$offset)
        list.dots$offset <- NULL
    } else
    {
        dtrain <- xgb.DMatrix(x[,-1],
                              label = y,
                              weight = wts)
    }

    # fit model
    cvfit <- do.call(xgb_cv_personalized, c(list.dots, list(params = params, data = dtrain,
                                                            obj_func = return_squared_error_loss_xgboost,
                                                            feval_func = return_eval_metric_xgboost,
                                                            trt.multiplier = trt.multiplier,
                                                            maximize = FALSE)))

    best.iter <- cvfit$best_iteration

    list.dots$nrounds <- list.dots$early_stopping_rounds <- list.dots$nfold <- NULL

    model <- do.call(xgb.train, c(list.dots, list(params = params,
                                                  data = dtrain,
                                                  nrounds = best.iter,
                                                  obj = return_squared_error_loss_xgboost(trt.multiplier))))


    # Return fitted model and extraction methods
    list(predict      = get.pred.func("fit_sq_loss_xgboost", model),
         model        = model,
         coefficients = get.coef.func("fit_sq_loss_xgboost")(model))
}


