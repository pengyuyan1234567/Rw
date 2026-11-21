#' Remove the seasonal cycle of a series
#'
#' The climatology of each label in `t_mmdd` is the mean of the observations
#' sharing that label. The anomalies are the deviations from their own
#' climatology, i.e. `y - mu` matched by `t_mmdd`.
#'
#' `NA` are dropped when averaging (`na.rm = TRUE`) and propagate to
#' `anomaly`.
#'
#' @param y Numeric vector of observations.
#' @param t_mmdd Integer vector, same length as `y`, giving the calendar
#'   position of each observation, typically `MMDD` (month-day) or the month
#'   number.
#'
#' @return A list of two numeric vectors:
#' - `mu`: the climatology, ordered by increasing `t_mmdd` and named after it;
#' - `anomaly`: the de-seasonalised series.
#'
#' @family piecewise regression helpers
#' @export
#' @examples
#' t_mmdd <- rep(1:12, each = 3)
#' y <- t_mmdd * 10 + rep(c(-1, 0, 1), times = 12)
#' RemoveSeasonal(y, t_mmdd)
RemoveSeasonal <- function(y, t_mmdd) {
  stopifnot(length(y) == length(t_mmdd))

  mu <- c(tapply(y, t_mmdd, mean, na.rm = TRUE))
  anomaly <- y - mu[as.character(t_mmdd)]
  list(mu = mu, anomaly = unname(anomaly))
}

#' Piecewise regression with one break point
#'
#' Fits `Y ~ X` allowing the intercept to shift at `Ic` while the slope is kept
#' constant on both segments, the classical step-change model of a
#' "same slope, different interception" series.
#'
#' @param Y Numeric response vector.
#' @param X Numeric predictor vector, usually time in the same length as `Y`.
#' @param Ic Integer position of the last observation before the break, so
#'   that `[1, Ic]` and `[Ic + 1, length(Y)]` are the two segments. Use `NA`
#'   for a model without any break point.
#' @param only.SSE Logical. If `TRUE`, only the residual sum of squares is
#'   returned; this is much cheaper when the fit is evaluated inside a search
#'   over candidates of `Ic`.
#'
#' @return A list with elements
#' - `coef`: the estimated coefficients,
#' - `fitted`: the fitted values,
#' - `resi`: the residuals,
#' - `SSE`: the residual sum of squares.
#'
#' Only `SSE` is returned when `only.SSE = TRUE`.
#'
#' @family piecewise regression helpers
#' @export
#' @examples
#' X <- seq_len(20)
#' Y <- c(seq_len(10), seq_len(10) + 15)
#' fit <- LS_OneBreak(Y, X, Ic = 10)
#' fit$coef
#' fit$SSE
LS_OneBreak <- function(Y, X, Ic, only.SSE = FALSE) {
  n <- length(Y)
  stopifnot(length(X) == n)

  design <- cbind(1, X)
  if (!is.na(Ic)) {
    design <- cbind(design, as.numeric(seq_len(n) > Ic))
  }
  lm_solve(design, Y, only.SSE = only.SSE)
}

#' Piecewise regression with multiple break points
#'
#' Extends [LS_OneBreak()] to any number of break points: each segment owns its
#' own intercept, while a single common slope is shared by all of them.
#'
#' @param Y Numeric response vector.
#' @param X Numeric predictor vector, same length as `Y`.
#' @param Ips Integer vector of the last position of each segment, in
#'   increasing order. Its last element must be `length(Y)`, so a two-segment
#'   series is specified as `c(Ic, length(Y))`; `length(Y)` alone means no break
#'   point at all.
#' @param only.SSE Logical, see [LS_OneBreak()].
#'
#' @return A list, see [LS_OneBreak()].
#'
#' @family piecewise regression helpers
#' @export
#' @examples
#' X <- seq_len(30)
#' Y <- c(seq_len(10), seq_len(10) + 5, seq_len(10) + 10)
#' fit <- LS_MultiBreaks(Y, X, Ips = c(10, 20, 30))
#' fit$coef
LS_MultiBreaks <- function(Y, X, Ips, only.SSE = FALSE) {
  n <- length(Y)
  stopifnot(
    length(X) == n,
    length(Ips) >= 1L,
    !is.unsorted(Ips),
    Ips[length(Ips)] == n
  )

  ## a new segment starts right after each break point
  segment <- 1L + cumsum(seq_len(n) %in% (Ips[-length(Ips)] + 1L))
  ## without a break point there is no dummy, factor(segment) would have one level
  dummies <- if (max(segment) > 1L) {
    stats::model.matrix(~ factor(segment))[, -1L, drop = FALSE]
  }
  design <- cbind(1, X - mean(X), dummies)
  lm_solve(design, Y, only.SSE = only.SSE)
}

#' Least squares on a given design matrix
#'
#' A thin wrapper around [stats::.lm.fit()] shared by the piecewise regression
#' helpers. It does not build a model frame, hence no `NA` handling and no
#' formula interface: the caller supplies a numeric design matrix.
#'
#' @param design Numeric design matrix, one row per observation.
#' @param response Numeric response vector.
#' @param only.SSE Logical, see [LS_OneBreak()].
#'
#' @return A list, see [LS_OneBreak()].
#'
#' @keywords internal
lm_solve <- function(design, response, only.SSE = FALSE) {
  fit <- stats::.lm.fit(design, as.matrix(response))
  SSE <- sum(fit$residuals^2)

  if (only.SSE) {
    return(list(SSE = SSE))
  }
  list(
    coef = fit$coefficients,
    fitted = as.vector(design %*% fit$coefficients),
    resi = drop(fit$residuals),
    SSE = SSE
  )
}
