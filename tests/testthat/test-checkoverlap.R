

context("check.overlap")

test_that("test plot is returned for hist/density/both", {


    set.seed(123)
    n.obs  <- 50
    n.vars <- 5
    x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)


    # simulate non-randomized treatment
    xbetat   <- 0.25 + 0.5 * x[,1] - 0.5 * x[,5]
    trt.prob <- exp(xbetat) / (1 + exp(xbetat))
    trt01    <- rbinom(n.obs, 1, prob = trt.prob)

    # create function for fitting propensity score model
    prop.func <- function(x, trt)
    {
        # fit propensity score model
        propens.model <- cv.glmnet(y = trt,
                                   x = x, family = "binomial")
        pi.x <- predict(propens.model, s = "lambda.min",
                        newx = x, type = "response")[,1]
        pi.x
    }

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "hist")

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "densitys")

    expect_is(pl, "ggplot")

    pl <- check.overlap(x = x,
                        trt = trt01,
                        propensity.func = prop.func,
                        type = "both")

    expect_is(pl, "ggplot")
})
