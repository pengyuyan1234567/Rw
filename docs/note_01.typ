#import "@local/modern-cug-report:0.1.3": *

#counter(heading).update(2)

#let size-config = define_size(
  12pt,
  (
    text: 12pt,
    math: (text: 12pt, block: 12pt),
    table: 13pt,
    figure: 10pt,
    raw: 11pt,
    heading: (H1: 13pt, H2: 13pt, H3: 13pt),
  ),
)
#show: doc => template(doc, size-config: size-config, footer: "CUG水文气象 2026", header: "Rw 包说明")

= 1 *引言*

气候与水文观测序列的均一性检验，往往要在成千上万个候选断点上反复拟合同一族模型：斜率不变、各段截距自由。拟合函数因此必须廉价，`lm()` 那一套 `model.frame` 与公式解析的开销在这里是纯浪费。

`Rw` 把这类回归直接写成设计矩阵上的最小二乘，交给 `stats::.lm.fit()` 求解；此外附带一个按日历标签去季节的辅助函数。包内只有四个函数，全部是纯函数，没有状态。

#box-blue[
  *核心思路*：分段回归 = 一次带指示变量的多元线性回归。断点一旦给定，模型就是线性的，剩下的只是把设计矩阵拼好。
]

= 2 *数学模型*

== 2.1 分段回归

设观测为 $(X_i, Y_i)$，$i = 1, dots, n$，候选断点为 $I_1 < dots < I_m$，其中 $I_m = n$。模型为

$
  Y_i = a + b X_i + sum_(j=1)^m c_j d_(j i) + e_i, quad d_(j i) = cases(1 & quad X_i > I_j, 0 & quad X_i <= I_j)
$

斜率 $b$ 全段共用，$c_j$ 是第 $j$ 段相对第一段的截距增量，$m = 1$（无断点）时模型退化为普通线性趋势。矩阵形式记作 $bold(Y) = bold(D) bold(beta) + bold(e)$，最小二乘解与残差平方和

$
  hat(bold(beta)) = (bold(D)' bold(D))^(-1) bold(D)' bold(Y), quad "SSE" = sum_i hat(e)_i^2
$

#box-green[
  `LS_MultiBreaks()` 估计时把自变量中心化，设计矩阵用 $tilde(X) = X - macron(X)$ 与 $bold(1)$ 拼成。这样斜率列与截距列正交，断点个数多时条件数不会爆炸，$c_j$ 的含义也更直观：#highlight()[每一段相对首段的截距差]。
]

== 2.2 去季节

给定日历标签（月日的 `MMDD`，或月份编号）$t_i$，记标签 $k$ 的气候态为

$
  mu_k = 1 / n_k sum_(i: t_i = k) y_i
$

其中 $n_k$ 为标签等于 $k$ 的观测个数。

去季节序列即距平 $a_i = y_i - mu_(t_i)$。

= 3 *函数一览*

#figure(
  caption: [导出的三个函数与一个内部求解器。],
  table(
    columns: (auto, 1fr, 1fr),
    align: horizon,
    table.header([*函数*], [*作用*], [*返回*]),
    [`lm_solve`], [设计矩阵上的最小二乘（内部）], [`coef`、`fitted`、`resi`、`SSE`],
    [`LS_OneBreak`], [单断点分段回归], [同上],
    [`LS_MultiBreaks`], [多断点分段回归], [同上],
    [`RemoveSeasonal`], [按日历标签去季节], [`mu`、`anomaly`],
  ),
)

=== 3.1 `LS_OneBreak(Y, X, Ic, only.SSE = FALSE)`

`Ic` 是断点前最后一个观测的位置，$[1, I_c]$ 与 $[I_c + 1, n]$ 为两段；取 `NA` 表示不含断点。系数向量依次为截距、斜率、截距增量。

=== 3.2 `LS_MultiBreaks(Y, X, Ips, only.SSE = FALSE)`

`Ips` 为各段末位，升序，末位须等于 `length(Y)`。内部以 `cumsum()` 一次性生成段编号，再用 `model.matrix(~ factor(segment))` 造出 $m - 1$ 个指示列。

```r
Ips <- c(10L, 20L, 30L)
segment <- 1L + cumsum(seq_len(30) %in% (Ips[-length(Ips)] + 1L))
```

=== 3.3 `lm_solve(design, response, only.SSE = FALSE)`

`stats::.lm.fit()` 的薄封装。它不建模型框，因此不处理 `NA`，也不接受公式；调用者负责给出数值设计矩阵。`only.SSE = TRUE` 时只返回残差平方和，供断点搜索内部使用——目标函数只需要一个标量。

=== 3.4 `RemoveSeasonal(y, t_mmdd)`

按标签求气候态，返回按标签升序、以标签命名的 `mu`，以及距平序列 `anomaly`。用 `tapply()` 一次算完，不做循环。

= 4 *使用示例*

```r
library(Rw)

X <- seq_len(20)
Y <- c(seq_len(10), seq_len(10) + 15)
n <- length(Y)

fit <- LS_OneBreak(Y, X, Ic = 10)
fit$coef          # 截距、斜率、截距增量
fit$SSE

LS_MultiBreaks(Y, X, Ips = c(10, 20))$coef

## 断点搜索：目标函数只要 SSE
sse <- \(Ic) LS_MultiBreaks(Y, X, c(Ic, n), only.SSE = TRUE)$SSE
which.min(vapply(1:19, sse, numeric(1)))   # 10

RemoveSeasonal(y = c(11, 12, 13, 21, 22, 23),
               t_mmdd = c(101, 101, 101, 102, 102, 102))
```

= 5 *实现要点*

+ *设计矩阵一次拼好*。断点给定后模型是线性的，没有理由再调用公式接口。
+ *开销与 $n$ 同阶*。`%in%`、`cumsum()` 均为线性，主要代价是 `.lm.fit()` 的 QR 分解。
+ *纯函数*。四个函数都不修改输入、不依赖全局状态，便于并行与重复调用。
+ *失败要早*。`stopifnot()` 检查长度、升序与末位位置，把非法断点在入口挡掉，而不是在 `:` 生成负长度序列时抛出难懂的报错。

= 6 *测试*

`tests/testthat/` 用 `testthat` 第三版，把每个函数与 `stats::lm()` 在等价公式下的结果对照：系数、拟合值、残差、`SSE` 全部 `expect_equal()`；此外覆盖无断点退化、`only.SSE` 一致性、观测次序无关性，以及非法 `Ips` 的报错。

```r
devtools::test()
```
