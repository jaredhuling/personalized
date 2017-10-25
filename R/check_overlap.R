..density.. <- NULL # need to make variable defined

#' Check propensity score overlap
#'
#' @description Results in a plot to check whether the propensity score has adequate overlap between treatment groups
#'
#' @param x The design matrix (not including intercept term)
#' @param trt treatment vector with each element equal to a 0 or a 1, with 1 indicating
#'            treatment status is active.
#' @param propensity.func function that inputs the design matrix x and the treatment vector trt and outputs
#' the propensity score, ie Pr(trt = 1 | X = x). Function should take two arguments 1) x and 2) trt. See example below.
#' For a randomized controlled trial this can simply be a function that returns a constant equal to the proportion
#' of patients assigned to the treatment group, i.e.:
#' \code{propensity.func = function(x, trt) 0.5}.
#' @param type Type of plot to create. Options are either a histogram (\code{type = "histogram"}) for each treatment
#' group, a density (\code{type = "density"}) for each treatment group, or to plot both a density and histogram
#' (\code{type = "code"})
#' @param bins integer number of bins for histograms when \code{type = "histogram"}
#' @importFrom ggplot2 guides guide_legend
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 250
#' n.vars <- 15
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.25 + 0.5 * x[,11] - 0.5 * x[,12]
#' trt.prob <- exp(xbetat) / (1 + exp(xbetat))
#' trt01    <- rbinom(n.obs, 1, prob = trt.prob)
#'
#' # create function for fitting propensity score model
#' prop.func <- function(x, trt)
#' {
#'     # fit propensity score model
#'     propens.model <- cv.glmnet(y = trt,
#'                                x = x, family = "binomial")
#'     pi.x <- predict(propens.model, s = "lambda.min",
#'                     newx = x, type = "response")[,1]
#'     pi.x
#' }
#'
#' check.overlap(x = x,
#'               trt = trt01,
#'               propensity.func = prop.func)
#'
#' # now add density plot with histogram
#' check.overlap(x = x,
#'               trt = trt01,
#'               type = "both",
#'               propensity.func = prop.func)
#'
#'
#' # simulated non-randomized treatment with multiple levels
#' xbetat_1   <- 0.15 + 0.5 * x[,9] - 0.25 * x[,12]
#' xbetat_2   <- 0.15 - 0.5 * x[,11] + 0.25 * x[,15]
#' trt.1.prob <- exp(xbetat_1) / (1 + exp(xbetat_1) + exp(xbetat_2))
#' trt.2.prob <- exp(xbetat_2) / (1 + exp(xbetat_1) + exp(xbetat_2))
#' trt.3.prob <- 1 - (trt.1.prob + trt.2.prob)
#' prob.mat <- cbind(trt.1.prob, trt.2.prob, trt.3.prob)
#' trt    <- apply(prob.mat, 1, function(rr) rmultinom(1, 1, prob = rr))
#' trt    <- apply(trt, 2, function(rr) which(rr == 1))
#'
#' # use multinomial logistic regression model with lasso penalty for propensity
#' propensity.multinom.lasso <- function(x, trt)
#' {
#'     if (!is.factor(trt)) trt <- as.factor(trt)
#'     gfit <- cv.glmnet(y = trt, x = x, family = "multinomial")
#'
#'     # predict returns a matrix of probabilities:
#'     # one column for each treatment level
#'     propens <- drop(predict(gfit, newx = x, type = "response", s = "lambda.min",
#'                             nfolds = 5, alpha = 0))
#'
#'     # return the probability corresponding to the
#'     # treatment that was observed
#'     probs <- propens[,match(levels(trt), colnames(propens))]
#'
#'     probs
#' }
#'
#' check.overlap(x = x,
#'               trt = trt,
#'               type = "histogram",
#'               propensity.func = propensity.multinom.lasso)
#'
#'
#'
#' @export
check.overlap <- function(x,
                          trt,
                          propensity.func,
                          type = c("histogram", "density", "both"),
                          bins = 50L)
{
    type <- match.arg(type)
    bins <- as.integer(bins[1])

    # compute propensity scores
    pi.x <- drop(propensity.func(x = x, trt = trt))

    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity.func() should return values between 0 and 1")

    # should be FALSE for treatment/control scenario,
    # TRUE for multiple treatment scenario
    multiplot <- FALSE

    dim.pi.x <- dim(pi.x)
    if (!is.null(dim.pi.x))
    {
        if (length(dim.pi.x) == 1)
        {
            pi.x <- as.vector(pi.x)
            prop.scores <- data.frame(Treatment = as.factor(trt), prop.score = pi.x)

        } else if (length(dim.pi.x) > 2)
        {
            stop("propensity.func() returns a multidimensional array; it can only return a vector or matrix.")
        }


        trt.names <- colnames(pi.x)

        if (is.null(trt.names))
        {
            if (is.factor(trt))
            {
                # drop any unused levels of trt
                trt         <- droplevels(trt)
                trt.names   <- levels(trt)
            } else
            {
                trt.names   <- sort(unique(trt))
            }

        }

        prop.scores <- data.frame(Treatment_Received = as.factor(rep(trt, NCOL(pi.x))),
                                  Treatment = rep(trt.names, NROW(trt)),
                                  prop.score = as.vector(t(pi.x)))

        levels(prop.scores$Treatment_Received) <- paste(levels(prop.scores$Treatment_Received), "Group")

        multiplot <- TRUE

    } else
    {
        prop.scores <- data.frame(Treatment = as.factor(trt), prop.score = pi.x)
    }


    Treatment <- prop.score <- NULL

    if (type == "density")
    {
        pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
            geom_density(alpha = 0.5, colour = "grey50") +
            geom_rug(aes(colour = Treatment)) +
            theme(legend.position = "bottom") +
            ggtitle("Densities of propensity scores by treatment group") +
            xlab("Propensity Score")
    } else if (type == "histogram")
    {
        pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
            geom_histogram(bins = bins, alpha = 0.5, position = "identity") +
            geom_rug(aes(colour = Treatment)) +
            theme(legend.position = "bottom") +
            ggtitle("Histograms of propensity scores by treatment group") +
            xlab("Propensity Score")
    } else
    {
        pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
            geom_histogram(aes(y = ..density..), bins = bins, alpha = 0.25, position = "identity") +
            geom_rug(aes(colour = Treatment)) +
            geom_density(alpha = 0.25) +
            theme(legend.position = "bottom") +
            ggtitle("Densities and histograms of propensity scores by treatment group") +
            xlab("Propensity Score")
    }

    if (multiplot)
    {
        pl.obj <- pl.obj + facet_grid(Treatment_Received ~ .)
    }
    pl.obj
}


