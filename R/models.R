#' Estimate truncated normal distribution
#'
#' Estimates the posterior modes for the mean (mu) and standard deviation
#' (sigma) of the underlying normal distribution, given truncated data with
#' known truncation point(s).
#'
#' @param x Vector of observations from truncated normal.
#' @param mu_start Initial value for mu.
#' @param sigma_start Initial value for sigma.
#' @param ci_level Confidence level of the interval – gives a 100*ci_level%
#'   symmetric HPD interval (defaults to 95%).
#' @param a Left truncation limit.
#' @param b Right truncation limit.
#' @param ... Parameters to pass to \code{rstan::sampling()}.
#'
#' @return A list with two elements:
#'  \describe{
#'    \item{stats}{A data frame with two rows and the columns \code{param}
#'                 (\code{mu}, \code{sd}), \code{mode} (posterior mode),
#'                 \code{mean} (posterior mean), \code{median}
#'                 (posterior median), \code{se} (standard error),
#'                 \code{ci_lower} (lower CI bound), \code{ci_upper}
#'                 (upper CI bound), \code{rhat}.}
#'    \item{fit}{A \code{stanfit} object (the result of fitting the model).}
#'  }
#'
#' @export
#'
#' @references
#' \insertRef{zhou2014}{truncnormbayes}
#'
#' \insertRef{stan2022}{truncnormbayes}
#'
#' @examples
#' set.seed(22)
#' x <- truncnorm::rtruncnorm(100, a = -1, b = 2, mean = 0.5, sd = 0.5)
#' trunc_est(x, a = -1, b = 2)
trunc_est <- function(x,
                      a,
                      b,
                      mu_start = 0,
                      sigma_start = 1,
                      ci_level = 0.95,
                      ...) {
  stopifnot(a < b)
  stopifnot(sigma_start > 0)
  stopifnot(all(x >= a))
  stopifnot(all(x <= b))

  # set start values for sampler
  init_fcn <- function() list(mean = mu_start, sd = sigma_start)

  stan_fit <- rstan::sampling(stanmodels$trunc_est,
                              init = init_fcn,
                              data = list(n = length(x), a = a, b = b, y = x),
                              ...)

  stan_extract <- rstan::extract(stan_fit)
  stan_summary <- as.data.frame(
    rstan::summary(stan_fit)$summary[c("mu", "sigma"), ]
  )
  means <- stan_summary$mean
  ses <- stan_summary$se_mean
  rhats <- stan_summary$Rhat

  medians <- c(median(stan_extract$mu), median(stan_extract$sigma))

  index_maxlp <- which.max(stan_extract$log_post)
  modes <- c(stan_extract$mu[index_maxlp], stan_extract$sigma[index_maxlp])

  stan_ci <- function(param, q) as.numeric(quantile(stan_extract[[param]], q))
  alpha <- 1 - ci_level
  cil <- c(stan_ci("mu", alpha / 2), stan_ci("sigma", alpha / 2))
  ciu <- c(stan_ci("mu", 1 - alpha / 2), stan_ci("sigma", 1 - alpha / 2))

  stan_stats <- data.frame(param = c("mu", "sigma"), mode = modes, mean = means,
                           median = medians, se = ses, ci_lower = cil,
                           ci_upper = ciu, rhat = rhats)

  return(list(stats = stan_stats, fit = stan_fit))
}
