


context("weighted.ksvm")

test_that("weighted.ksvm fitting", {

    library(kernlab)

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


    foldid <- sample(rep(seq(5), length = 100))

    wk <- weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(1, 3),
                        foldid = foldid,
                        weights = weights[1:100])

    expect_is(wk, "wksvm")

    if (Sys.info()[[1]] != "windows")
    {

        expect_error(weighted.ksvm(x = x[1:100,], y = y[1:100], C = c(0.1),
                                   nfolds = 150,
                                   weights = weights[1:100]))

        wk <- weighted.ksvm(x = x[1:100,], y = as.factor(y[1:100]), C = c(1, 3),
                            foldid = foldid,
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        expect_error(weighted.ksvm(x = x[1:100,], y = c(1:5, y[5:100]), C = c(0.1),
                                   weights = weights[1:100]))


        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 3),
                            foldid = foldid,
                            weights = weights[1:100])

        expect_is(wk, "wksvm")


        wk <- weighted.ksvm(x = x[1:100,], y = as.factor(y[1:100]), C = c(1, 3),
                            foldid = foldid,
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        expect_warning(weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 3),
                                     nfolds = -5,
                                     weights = weights[1:100]))

        expect_error(weighted.ksvm(x = x[1:100,], y = y[1:100]/2 + 0.5, C = c(0.1),
                                   weights = weights[1:100]))




        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "polydot",
                            weights = weights[1:100])

        expect_is(wk, "wksvm")


        wk <- weighted.ksvm(x = x[1:100,], y = as.factor(y[1:100]), C = c(1, 3),
                            foldid = foldid,
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(10),
                            foldid = foldid,
                            kernel = "tanhdot",
                            weights = rep(1, 100),
                            margin = 0.5,
                            bound = 10,
                            maxiter = 200)

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "vanilladot",
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "laplacedot",
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "besseldot",
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:25,], y = as.character(y[1:25]), C = c(1, 10),
                            kernel = "tanhdot", maxiter = 500, bound = 10,
                            weights = weights[1:25])

        expect_is(wk, "wksvm")

        summary(wk)

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "anovadot",
                            weights = weights[1:100])

        expect_is(wk, "wksvm")

        wk <- weighted.ksvm(x = x[1:100,], y = as.character(y[1:100]), C = c(1, 10),
                            foldid = foldid,
                            kernel = "splinedot",
                            weights = weights[1:100],
                            margin = 0.1,
                            maxiter = 100)

        expect_is(wk, "wksvm")
    }



})
