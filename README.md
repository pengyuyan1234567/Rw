
# Rw

<!-- badges: start -->
[![R-CMD-check](https://github.com/pengyuyan1234567/Rw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pengyuyan1234567/Rw/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`Rw` provides fast piecewise ("segmented") regression with one or several
break points, where the slope is shared by all segments and each segment has
its own intercept, plus a small helper to remove the seasonal cycle of a
series.

The regression routines work directly on a design matrix built with
`stats::.lm.fit()`, so they are cheap enough to be called thousands of times
while searching for the best break points.

## Installation

You can install the development version of `Rw` like so:

``` r
remotes::install_github("pengyuyan1234567/Rw")
```

## Example

``` r
library(Rw)

set.seed(1)
X <- seq_len(20)
Y <- c(seq_len(10), seq_len(10) + 15) + rnorm(20, sd = 0.4)

## one break point, shared slope
fit <- LS_OneBreak(Y, X, Ic = 10)
round(fit$coef, 3)
#> [1] -0.061  1.021  4.840
fit$SSE
#> [1] 2.454179

## several break points
LS_MultiBreaks(Y, X, Ips = c(10, 20))$coef

## remove the seasonal cycle
RemoveSeasonal(y = c(11, 12, 13, 21, 22, 23), t_mmdd = c(101, 101, 101, 102, 102, 102))
#> $mu
#> 101 102
#>  12  22
#>
#> $anomaly
#> [1] -1  0  1 -1  0  1
```

## Documentation

The model equations, the design decisions and a walk-through of the four
functions are written up as a Typst note in `docs/note_01.typ`:

``` sh
typst compile --root . docs/note_01.typ
```