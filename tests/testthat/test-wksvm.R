


context("weighted.ksvm")

test_that("weighted.ksvm fitting", {


    set.seed(123)

    x <- matrix(rnorm(200 * 2), ncol = 2)

    y <- 2 * (sin(x[,2]) ^ 2 * exp(-x[,2]) > rnorm(200, sd = 0.1)) - 1

    weights <- runif(100, max = 1.5, min = 0.5)

    wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(0.1, 0.5, 1),
                        weights = weights[1:100])

    expect_is(wk, "wksvm")

    pr <- predict(wk, newx = x[1:100,])

    wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = 1,
                        weights = weights[1:100])

    expect_is(wk, "wksvm")

    expect_error(weighted.ksvm(x = x[1:101,], y = y[1:100], C = c(10),
                        weights = weights[1:100]))

    expect_error(weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(0.1),
                               weights = weights[1:101]))


    foldid <- sample(rep(seq(5), length = NROW(y)))

    wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(1, 3),
                        foldid = foldid,
                        weights = weights[1:100])

    expect_is(wk, "wksvm")


})
