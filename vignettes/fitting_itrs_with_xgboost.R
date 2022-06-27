## ----loadpkg, message=FALSE, warning=FALSE, echo=FALSE------------------------
library(personalized)

## ----sim_data_1, message = FALSE, warning = FALSE-----------------------------
library(personalized)

set.seed(1)
n.obs  <- 500
n.vars <- 10
x <- matrix(rnorm(n.obs * n.vars, sd = 1), n.obs, n.vars)

# simulate non-randomized treatment
xbetat   <- 0.5 + 0.25 * x[,1] - 0.25 * x[,5]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt      <- rbinom(n.obs, 1, prob = trt.prob)

# simulate delta (CATE) as a complex function of the covariates
delta <- 2*(0.25 + x[,1] * x[,2] - x[,3] ^ {-2} * (x[,3] > 0.35) + 
                (x[,1] < x[,3]) - (x[,1] < x[,2]))

# simulate main effects g(X)
xbeta <- x[,1] + x[,2] + x[,4] - 0.2 * x[,4]^2 + x[,5] + 0.2 * x[,9] ^ 2
xbeta <- xbeta + delta * (2 * trt - 1) * 0.5

# simulate continuous outcomes
y <- drop(xbeta) + rnorm(n.obs)

## -----------------------------------------------------------------------------
# arguments can be passed to cv.glmnet via `cv.glmnet.args`.
# normally we would set the number of crossfit folds and internal folds to be larger, 
# but have reduced it to make computation time shorter
prop.func <- create.propensity.function(crossfit = TRUE,
                                        nfolds.crossfit = 4,
                                        cv.glmnet.args = list(type.measure = "auc", nfolds = 3))

## -----------------------------------------------------------------------------
aug.func <- create.augmentation.function(family = "gaussian",
                                         crossfit = TRUE,
                                         nfolds.crossfit = 4,
                                         cv.glmnet.args = list(type.measure = "mse", nfolds = 3))

## -----------------------------------------------------------------------------
subgrp.model.linear <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             augment.func = aug.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 3)    # option for cv.glmnet (for ITR estimation)

summary(subgrp.model.linear)

## -----------------------------------------------------------------------------

## xgboost tuning parameters to use:
param <- list(max_depth = 5, eta = 0.01, nthread = 1, 
              booster = "gbtree", subsample = 0.623, colsample_bytree = 1)

subgrp.model.xgb <- fit.subgroup(x = x, y = y,
                        trt = trt,
                        propensity.func = prop.func,
                        augment.func = aug.func,
                        ## specify xgboost via the 'loss' argument
                        loss   = "sq_loss_xgboost",
                        nfold = 3,
                        params = param, verbose = 0,
                        nrounds = 5000, early_stopping_rounds = 50)

subgrp.model.xgb

## ----eval=FALSE---------------------------------------------------------------
#  valmod.lin <- validate.subgroup(subgrp.model.linear, B = 100,
#                              method = "training_test",
#                              train.fraction = 0.75)
#  valmod.lin

## ----eval=FALSE---------------------------------------------------------------
#  valmod.xgb <- validate.subgroup(subgrp.model.xgb, B = 100,
#                                  method = "training_test",
#                                  train.fraction = 0.75)
#  valmod.xgb

## ---- fig.height=10-----------------------------------------------------------

## RMSE (note: this is still on the in-sample data so
## out-of-sample RMSE is preferred to evaluate methods)

sqrt(mean((delta - treatment.effects(subgrp.model.linear)$delta) ^ 2))
sqrt(mean((delta - treatment.effects(subgrp.model.xgb)$delta) ^ 2))

par(mfrow = c(2,1))
plot(delta ~ treatment.effects(subgrp.model.linear)$delta,
     ylab = "True CATE", xlab = "Estimated Linear CATE")
abline(a=0,b=1,col="blue")
plot(delta ~ treatment.effects(subgrp.model.xgb)$delta,
     ylab = "True CATE", xlab = "CATE via xgboost") 
abline(a=0,b=1,col="blue")

