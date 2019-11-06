## ----loadpkg, message=FALSE, warning=FALSE-------------------------------
library(personalized)

## ----sim_three_trt_data--------------------------------------------------
set.seed(123)
n.obs  <- 250
n.vars <- 15
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)

# simulated non-randomized treatment with multiple levels
# based off of a multinomial logistic model
xbetat_1 <- 0.1 + 0.5 * x[,1]  - 0.25 * x[,5]
xbetat_2 <- 0.1 - 0.5 * x[,11] + 0.25 * x[,15]
trt.1.prob <- exp(xbetat_1) / (1 + exp(xbetat_1) + exp(xbetat_2))
trt.2.prob <- exp(xbetat_2) / (1 + exp(xbetat_1) + exp(xbetat_2))
trt.3.prob <- 1 - (trt.1.prob + trt.2.prob)

prob.mat <- cbind(trt.1.prob, trt.2.prob, trt.3.prob)
trt.mat <- apply(prob.mat, 1, function(rr) rmultinom(1, 1, prob = rr))
trt.num <- apply(trt.mat, 2, function(rr) which(rr == 1))
trt <- as.factor(paste0("Trt_", trt.num))

# simulate response

# effect of treatment 1 relative to treatment 3
delta1 <- 2 * (0.5 + x[,2] - 2 * x[,3]  )
# effect of treatment 2 relative to treatment 3
delta2 <- (0.5 + x[,6] - 2 * x[,5] )

# main covariate effects with nonlinearities
xbeta <- x[,1] + x[,11] - 2 * x[,12]^2 + x[,13] + 
    0.5 * x[,15] ^ 2 + 2 * x[,2] - 3 * x[,5]

# create entire functional form of E(Y|T,X)
xbeta <- xbeta + 
    delta1 * ((trt.num == 1) - (trt.num == 3) ) + 
    delta2 * ((trt.num == 2) - (trt.num == 3) )


# simulate continuous outcomes E(Y|T,X)
y <- xbeta + rnorm(n.obs, sd = 2)


## ----display_mult_trt_vector---------------------------------------------
trt[1:5]
table(trt)

## ----define_multi_propens------------------------------------------------
propensity.multinom.lasso <- function(x, trt)
{
    if (!is.factor(trt)) trt <- as.factor(trt)
    gfit <- cv.glmnet(y = trt, x = x, family = "multinomial")

    # predict returns a matrix of probabilities:
    # one column for each treatment level
    propens <- drop(predict(gfit, newx = x, 
        type = "response", s = "lambda.min"))

    # return the matrix probability of treatment assignments
    probs <- propens[,match(levels(trt), colnames(propens))]

    probs
}

## ----check_overlap_multitreat, fig.cap = "Propensity score overlap plot for multi-category treatment data."----
check.overlap(x = x, trt = trt, propensity.multinom.lasso)

## ----fit_multi_trt_model-------------------------------------------------
set.seed(123)
subgrp.multi <- fit.subgroup(x = x, y = y,
    trt = trt, propensity.func = propensity.multinom.lasso,
    reference.trt = "Trt_3",
    loss   = "sq_loss_lasso")

summary(subgrp.multi)

## ----plot_multi_trt_model, warning=FALSE, message=FALSE, fig.width=5, fig.cap = "Individual outcome observations by treatment group and subgroup."----
pl <- plot(subgrp.multi)
pl + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----validate_multi_trt_model, eval = TRUE-------------------------------
set.seed(123)
validation.multi <- validate.subgroup(subgrp.multi, 
    B = 5,  # specify the number of replications
    method = "training_test_replication",
    train.fraction = 0.5)

print(validation.multi, digits = 2, sample.pct = TRUE)

## ----plotcomparemultivalidated, warning=FALSE, message=FALSE, eval = TRUE, fig.width=5, fig.cap = "Validation results for multi-category treatment data."----
plv <- plot(validation.multi)
plv + theme(axis.text.x = element_text(angle = 90, hjust = 1))

## ----multinom_propens----------------------------------------------------
propensity.func.multinom <- function(x, trt)
{
    df <- data.frame(trt = trt, x)
    mfit <- nnet::multinom(trt ~ . -trt, data = df)
    # predict returns a matrix of probabilities:
    # one column for each treatment level
    propens <- nnet::predict.nnet(mfit, type = "probs")

    if (is.factor(trt))
    {
        values <- levels(trt)[trt]
    } else 
    {
        values <- trt
    }
    # return the probability corresponding to the
    # treatment that was observed
    probs <- propens[cbind(1:nrow(propens), 
        match(values, colnames(propens)))]
    probs
}

## ----multinom_propens2---------------------------------------------------
propensity.func.multinom <- function(x, trt)
{
    require(nnet)
    df <- data.frame(trt = trt, x)
    mfit <- multinom(trt ~ . -trt, data = df)
    # predict returns a matrix of probabilities:
    # one column for each treatment level
    propens <- predict(mfit, type = "probs")

    if (is.factor(trt))
    {
        levels <- levels(trt)
    } else 
    {
        levels <- sort(unique(trt))
    }
    # return the probability corresponding to the
    # treatment that was observed
    probs <- propens[,match(levels, colnames(propens))]
    probs
}

