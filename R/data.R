#' National Supported Work Study Data
#'
#'
#' @description The LaLonde dataset comes from the National Supported Work Study, which sought to
#' evaluate the effectiveness of an employment trainining program on wage increases.
#' @format A data frame with 722 observations and 12 variables:
#' \describe{
#'   \item{outcome}{whether earnings in 1978 are larger than in 1975; 1 for yes, 0 for no}
#'   \item{treat}{whether the individual received the treatment; "Yes" or "No"}
#'   \item{age}{age in years}
#'   \item{educ}{education in years}
#'   \item{black}{black or not; factor with levels "Yes" or "No"}
#'   \item{hisp}{hispanic or not; factor with levels "Yes" or "No"}
#'   \item{white}{white or not; factor with levels "Yes" or "No"}
#'   \item{marr}{married or not; factor with levels "Yes" or "No"}
#'   \item{nodegr}{No high school degree; factor with levels "Yes" (for no HS degree) or "No"}
#'   \item{log.re75}{log of earnings in 1975}
#'   \item{u75}{unemployed in 1975; factor with levels "Yes" or "No"}
#'   \item{wts.extrap}{extrapolation weights to the 1978 Panel Study for Income Dynamics dataset}
#' }
#' @source The National Supported Work Study.
#' @references LaLonde, R.J. 1986. "Evaluating the econometric evaulations of training programs with experimental data." American Economic Review, Vol.76, No.4, pp. 604-620.
#'
#' Egami  N,  Ratkovic  M,  Imai  K  (2017). "\pkg{FindIt}:  Finding  Heterogeneous  Treatment  Effects." \code{R} package version 1.1.2, \url{https://CRAN.R-project.org/package=FindIt}.
#' @examples
#' data(LaLonde)
#' y <- LaLonde$outcome
#'
# treatment assignment (employment training vs not)
#' trt <- LaLonde$treat
#'
#' x.varnames <- c("age", "educ", "black", "hisp", "white",
#'                 "marr", "nodegr", "log.re75", "u75")
#'
#' # covariates
#' data.x <- LaLonde[, x.varnames]
#'
#' # construct design matrix (with no intercept)
#' x <- model.matrix(~ -1 + ., data = data.x)
#'
#' const.propens <- function(x, trt)
#' {
#'     mean.trt <- mean(trt == "Trt")
#'     rep(mean.trt, length(trt))
#' }
#'
#' subgrp_fit_w <- fit.subgroup(x = x, y = y, trt = trt,
#'     loss = "logistic_loss_lasso",
#'     propensity.func = const.propens,
#'     cutpoint = 0,
#'     type.measure = "auc",
#'     nfolds = 10)
#'
#' summary(subgrp_fit_w)
"LaLonde"
