

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
#' @examples
#' library(personalized)
#'
#' set.seed(123)
#' n.obs  <- 1000
#' n.vars <- 50
#' x <- matrix(rnorm(n.obs * n.vars, sd = 3), n.obs, n.vars)
#'
#'
#' # simulate non-randomized treatment
#' xbetat   <- 0.5 + 0.5 * x[,21] - 0.5 * x[,41]
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
    pi.x <- propensity.func(x = x, trt = trt)

    # make sure the resulting propensity scores are in the
    # acceptable range (ie 0-1)
    rng.pi <- range(pi.x)

    if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity.func() should return values between 0 and 1")

    prop.scores <- data.frame(Treatment = as.factor(trt), prop.score = pi.x)

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
    pl.obj
}
