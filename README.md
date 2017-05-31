





[![Build Status](https://travis-ci.org/jaredhuling/personalized.svg?branch=master)](https://travis-ci.org/jaredhuling/personalized)


## Introduction to the 'personalized' package

The personalized package provides estimation methods for subgroup identification under the framework of [Chen et al (2017)](http://onlinelibrary.wiley.com/doi/10.1111/biom.12676/abstract)

Install using the **devtools** package:


```r
devtools::install_github("jaredhuling/personalized")
```

or by cloning and building using `R CMD INSTALL`

## Quick Usage Reference

Load the package and access help files for the main functions:

```r
library(personalized)
```


```r
?fit.subgroup
?validate.subgroup
```


```r
set.seed(123)
n.obs  <- 250
n.vars <- 100
x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


# simulate non-randomized treatment
xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41] + 0.5 * x[,1] * x[,12]
trt.prob <- exp(xbetat) / (1 + exp(xbetat))
trt01    <- rbinom(n.obs, 1, prob = trt.prob)

trt      <- 2 * trt01 - 1

# simulate response
delta <- (0.5 + x[,2] - x[,3] - x[,11])
xbeta <- x[,1] + x[,11] - 0.5 * x[,12]^2 + x[,13] + 0.5 * x[,15] ^ 2
xbeta <- xbeta + delta * trt

# continuous outcomes
y <- drop(xbeta) + rnorm(n.obs, sd = 2)

# create function for fitting propensity score model
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    
    # return predicted Pr(trt = 1 | X = x)
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}

subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt01,
                             propensity.func = prop.func,
                             method = "a_learning",
                             loss   = "sq_loss_lasso",
                             nfolds = 10)              # option for cv.glmnet

summary(subgrp.model)
```

```
## family:  gaussian 
## loss:    sq_loss_lasso 
## method:  a_learning 
## 
## Average Outcomes:
##                Recommended Trt Recommended Ctrl
## Received Trt   4.4917 (n = 74)  -4.639 (n = 68)
## Received Ctrl -6.1364 (n = 62)  2.1427 (n = 46)
## 
##  Trt  Effect Among Recommended Trt Ctrl Effect Among Recommended Ctrl 
##                  10.6281 (n = 136)                   6.7817 (n = 114) 
## 
## Benefit score quantiles: 
##       0%      25%      50%      75%     100% 
## -15.5465  -4.1472   0.5254   5.4670  17.7227 
## 
## 7 variables selected by the lasso (cross validation criterion).
## 
##        Estimate
## V2   1.24291538
## V3  -1.30705442
## V6   0.32361980
## V7  -0.57852824
## V11 -0.96475524
## V29  0.06694316
## V33 -0.07032150
```


```r
plot(subgrp.model)
```

![](vignettes/plot_model-1.png)<!-- -->


