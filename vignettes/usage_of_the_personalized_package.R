## ----sim_data_1, message = FALSE, warning = FALSE------------------------
library(personalized)

set.seed(123)
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

## ----create_propensity---------------------------------------------------
# create function for fitting propensity score model
prop.func <- function(x, trt)
{
 # fit propensity score model
 propens.model <- cv.glmnet(y = trt,
                            x = x, 
                            family = "binomial")
 pi.x <- predict(propens.model, s = "lambda.min",
                 newx = x, type = "response")[,1]
 pi.x
}

## ----plot_overlap--------------------------------------------------------
check.overlap(x, trt, prop.func)

## ----fit_model-----------------------------------------------------------
subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 10)              # option for cv.glmnet

summary(subgrp.model)

## ----plot_model----------------------------------------------------------
plot(subgrp.model)

## ----plot_model_2--------------------------------------------------------
plot(subgrp.model, type = "interaction")

## ----validate_model------------------------------------------------------
validation <- validate.subgroup(subgrp.model, 
                                B = 25L,  # specify the number of replications
                                method = "training_test_replication",
                                train.fraction = 0.75)

validation

## ----plot_validation-----------------------------------------------------
plot(validation)

## ----plot_validation_2---------------------------------------------------
plot(validation, type = "interaction")

## ----plot_validation_compare---------------------------------------------
plotCompare(subgrp.model, validation, type = "interaction")

## ----propens_func_example_glm--------------------------------------------
propensity.func <- function(x, trt)
{
    # save data in a data.frame
    data.fr <- data.frame(trt = trt, x)
    
    # fit propensity score model
    propensity.model <- glm(trt ~ ., family = binomial(), data = data.fr)
    
    # create estimated probabilities
    pi.x <- predict(propensity.model, type = "response")
    return(pi.x)
}

propensity.func(x, trt)[101:105]
trt[101:105]

## ----propens_func_example_const------------------------------------------
propensity.func <- function(x, trt) 0.5

## ----fit_model2----------------------------------------------------------
subgrp.model2 <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso_gam",
                             nfolds = 10)              # option for cv.glmnet

summary(subgrp.model2)

## ----binary_example------------------------------------------------------
# create binary outcomes
y.binary <- 1 * (xbeta + rnorm(n.obs, sd = 2) > 0 )


## ----fit_binary_1--------------------------------------------------------
subgrp.bin <- fit.subgroup(x = x, y = y.binary,
                           trt = trt,
                           propensity.func = prop.func,
                           loss   = "logistic_loss_lasso",
                           nfolds = 10)      # option for cv.glmnet

## ----fit_binary_2, eval = FALSE------------------------------------------
#  subgrp.bin2 <- fit.subgroup(x = x, y = y.binary,
#                              trt = trt,
#                              propensity.func = prop.func,
#                              loss = "logistic_loss_gbm",
#                              shrinkage = 0.025,  # options for gbm
#                              n.trees = 1500,
#                              interaction.depth = 3,
#                              cv.folds = 5)

## ----plotcompare_bin-----------------------------------------------------
subgrp.bin

## ----tte_example---------------------------------------------------------
# create time-to-event outcomes
surv.time <- exp(-20 - xbeta + rnorm(n.obs, sd = 1))
cens.time <- exp(rnorm(n.obs, sd = 3))
y.time.to.event  <- pmin(surv.time, cens.time)
status           <- 1 * (surv.time <= cens.time)

## ----tte_model_example---------------------------------------------------
library(survival)
set.seed(123)
subgrp.cox <- fit.subgroup(x = x, y = Surv(y.time.to.event, status),
                           trt = trt,
                           propensity.func = prop.func,
                           method = "weighting",
                           loss   = "cox_loss_lasso",
                           nfolds = 10)      # option for cv.glmnet

## ----print_tte_model-----------------------------------------------------
summary(subgrp.cox)

## ----efficiency_func_example---------------------------------------------
adjustment.func <- function(x, y)
{
    df.x  <- data.frame(x)
    
    # add all squared terms to model
    form  <- eval(paste(" ~ -1 + ", 
                paste(paste('poly(', colnames(df.x), ', 2)', sep=''), 
                      collapse=" + ")))
    mm    <- model.matrix(as.formula(form), data = df.x)
    cvmod <- cv.glmnet(y = y, x = mm, nfolds = 10)
    predictions <- predict(cvmod, newx = mm, s = "lambda.min")
    predictions
}

## ----efficiency_ex-------------------------------------------------------

subgrp.model.eff <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             augment.func = adjustment.func,
                             nfolds = 10)              # option for cv.glmnet

summary(subgrp.model.eff)

## ----plot_ex_model_1-----------------------------------------------------
plot(subgrp.model)

## ----plot_ex_model_2-----------------------------------------------------
plot(subgrp.model, type = "density")

## ----plot_ex_model_3-----------------------------------------------------
plot(subgrp.model, type = "interaction")

## ----plot_compare_ex2----------------------------------------------------
plotCompare(subgrp.model, subgrp.model.eff)

## ----summarize_sub-------------------------------------------------------
comp <- summarize.subgroups(subgrp.model)

## ----print_summarize-----------------------------------------------------
print(comp, p.value = 0.01)

## ----summarize_sub2------------------------------------------------------
comp2 <- summarize.subgroups(x, subgroup = subgrp.model$benefit.scores > 0)

## ----train_test_ex-------------------------------------------------------

# check that the object is an object returned by fit.subgroup()
class(subgrp.model.eff)

validation.eff <- validate.subgroup(subgrp.model.eff, 
                                 B = 25L,  # specify the number of replications
                                 method = "training_test_replication",
                                 train.fraction = 0.75)

validation.eff

## ----boot_bias_ex--------------------------------------------------------
validation3 <- validate.subgroup(subgrp.model, 
                                 B = 25L,  # specify the number of replications
                                 method = "boot_bias_correction")

validation3

## ----plot_ex_model_1a----------------------------------------------------
plot(validation)

## ----plot_ex_model_1b----------------------------------------------------
plot(validation, type = "density")

## ----plot_compare_ex3, eval = FALSE--------------------------------------
#  plotCompare(validation, validation.eff)

