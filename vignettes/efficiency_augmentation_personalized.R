## ----loadpkg, message=FALSE, warning=FALSE-------------------------------
library(personalized)

## ----sim_data_1, message = FALSE, warning = FALSE------------------------
library(personalized)

set.seed(1)
n.obs  <- 1000
n.vars <- 50
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)

# simulate non-randomized treatment
xbetat   <- 0.5 + 0.25 * x[,21] - 0.25 * x[,41]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt      <- rbinom(n.obs, 1, prob = trt.prob)

# simulate delta
delta <- (0.5 + x[,2] - 0.5 * x[,3] - 1 * x[,11] + 1 * x[,1] * x[,12] )

# simulate main effects g(X)
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
xbeta <- xbeta + delta * (2 * trt - 1)

# simulate continuous outcomes
y <- drop(xbeta) + rnorm(n.obs)

## ------------------------------------------------------------------------
# arguments can be passed to cv.glmnet via `cv.glmnet.args`
prop.func <- create.propensity.function(crossfit = TRUE,
                                        nfolds.crossfit = 10,
                                        cv.glmnet.args = list(type.measure = "auc", nfolds = 5))

## ------------------------------------------------------------------------
subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 10)    # option for cv.glmnet (for ITR estimation)

summary(subgrp.model)

## ------------------------------------------------------------------------
aug.func <- create.augmentation.function(family = "gaussian",
                                         crossfit = TRUE,
                                         nfolds.crossfit = 10,
                                         cv.glmnet.args = list(type.measure = "mae", nfolds = 5))

## ------------------------------------------------------------------------
subgrp.model.aug <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             augment.func = aug.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 10)    # option for cv.glmnet (for ITR estimation)

summary(subgrp.model.aug)

## ------------------------------------------------------------------------
valmod <- validate.subgroup(subgrp.model, B = 5,
                            method = "training_test",
                            train.fraction = 0.75)
valmod

## ------------------------------------------------------------------------
valmod.aug <- validate.subgroup(subgrp.model.aug, B = 5,
                                method = "training_test",
                                train.fraction = 0.75)
valmod.aug

