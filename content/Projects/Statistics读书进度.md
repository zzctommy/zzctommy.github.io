---
title: Statistics读书进度
---

## 进度

起到一个鞭策自己看书的作用，每天读一点记下来也许有成就感。


| Date        | Nonparametric Regression(190) | Casual Inferene(118) |
| ----------- | ----------- | --------------|
| ?-26.10.19    | 15 | 2 |
| 26.10.20    | 19 | 12 |
| 26.10.21    | 30 | 18 |
| 26.10.22    | 30 | 28 |



------

## 每日总结

我意识到应该每节学完之后写一个概括，否则我会迷失在技术细节里。

### Nonparametric Regression

#### 1.1

非参（参数甚至不可数）统计的典例，例如从函数空间recover probabilistic density。

#### 1.2

引入Kernel density estimator $\hat{p_n}=\frac{1}{nh}\sum_i K((X_i-x)/h)$，function class用于限制bias，variance比较好控制，主要结果有： MISE（MSE的推广） $\sup_\limits{x_0\in \mathbb{R}}\sup\limits_{p\in \mathcal{P}(\beta,L)}\mathbb{E}_p[(\hat{p_n}(x_0)-p(x_0))^2]\le C n^{-2\beta/(2\beta+1)}$ 。

#### 1.3

通过Fourier transform给出更精确的MISE估计，因为求导在Fourier transform下的简单表示主要考虑Sobolev class，主要结果是：$L^2$ 且对称的kernel，如果 $ess\sup_{t\not=0}\frac{1-\hat{K}(t)}{|t|^{\beta}}\le A$ ，那么 $\sup_{p\in \mathcal{P_s}(\beta,L)}\mathbb{E}_p[(\hat{p_n}(x_0)-p(x_0))^2]\le C n^{-2\beta/(2\beta+1)}$ 。

#### 1.4

cross-validation首次出现（看了因果推断才懂，不过这里只是构造estimator而不是减少过拟合问题）。目标是寻找恰当的h给出好的estimator。因为 MISE(h)依赖未知density p，这一节给出完全依赖data的unbiased estimator CV(h) of $MISE(h)-\int_\mathbb{R}p^2(x)dx$，所以问题变成最小化CV(h)。另外还有从Fourier transform角度出发构造的另一个 $MISE(h)-\int_\mathbb{R}p^2(x)dx$ 的无偏统计量。

#### 1.5

目标是估计回归函数 $f(x)=\mathbb{E}[Y|X=x]$，构造了Nadaraya-Watson estimator $f_n^{NW}(x)=\frac{\sum_i Y_iK((X_i-x)/h)}{\sum_i K((X_i-x)/h)}$(for random design)。fix design没太看懂，没细讲，遇到再说。

### Casual Inference

#### Lec1

ATE（average treatment effect） $\tau=\mathbb{E}[Y_i(1)-Y_i(0)]$。

$Y_i(W_i)$：outcome，是否康复，$W_i$ 是否接受治疗。

理想化情况 $\tau_{SATE}=\frac{1}{n}\sum_iY_i(1)-Y_i(0)$ （sample average treatment effect）就是一个估计，但是在 $\{Y_i(0),Y_i(1)\}$ 中只能选一个观测。

Difference-in-means estimator $\hat{\tau}_{DM}=\frac{1}{n_1}\sum_{W_i=1}Y_i-\frac{1}{n_0}\sum_{W_i=0}Y_i$。

后面拿linear model和nonlinear model的OLS方法（nonlinear用Huber-White）算了example，OLS给出的estimator都小于difference-in-means estimator，所以具体问题中DM往往不是最优的

#### Lec2

整体比较会丢失因果，因此通过加入协变量 $X_i$ 给被试者分组，探测协变量带来的因果。

group-wise ATE $\tau(x)=\mathbb{E}[Y_i(1)-Y_i(0)|X_i=x]$。

Aggregating difference-in-means estimator: $\hat{\tau}_{AGG}=\sum_{x\in \chi}\frac{n_x}{n}\hat{\tau}(x)$ where $\hat{\tau}(x)=\frac{1}{n_{x,1}}\sum\limits_{X_i=x,W_i=1}Y_i-\frac{1}{n_{x,0}}\sum\limits_{X_i=x,W_i=0}Y_i$，这个estimator渐进正态，要求协变量离散（要求每个协变量被重复观测的次数够多才精确）。

unconfoundedness：$\{Y_i(0),Y_i(1)\}\perp W_i\big|X_i$，也就是在协变量分组下，潜在结果与是否被治疗无关。

propensity score：$e(x)=\mathbb{P}[W_i=1\big|X_i=x]$。

在unconfoundedness成立的条件下，$\{Y_i(0),Y_i(1)\}\perp W_i\big|e(X_i)$。

Propensity stratification：按照propensity score的大小分为若干组，并把每一组的propensity score平均化处理成一样的，跑DM estimate，通过平衡组的个数（组越多对每一组propensity score的估计就越准；组越少每一组的数据就越多，DM estimate更准）来构造estimator。

IPW（inverse propensity weighting）$\hat{\tau}_{IPW}=\frac{1}{n}\sum_i \frac{W_iY_i}{\hat{e}(X_i)}-\frac{(1-W_i)Y_i}{1-\hat{e}(X_i)}$ （$\hat{e}$ 用非参回归计算）

oracle IPW：$\hat{\tau}_{IPW}^*=\frac{1}{n}\sum_i \frac{W_iY_i}{e(X_i)}-\frac{(1-W_i)Y_i}{1-e(X_i)}$

在unconfoundedness下，$\hat{\tau}_{IPW}$ 是ATE的unbiased estiamtor，$\frac{1}{\sqrt{n}}$ rate逼近 $\tau$。

协变量离散的时候 $\hat{\tau}_{AGG}$ 比 $\hat{\tau}_{IPW}^*$ 精确（方差更小）。但是由于IPW更“抽象”，只需要考虑propensity score的性质，所以更加可以研究。

#### Lec3

notations：$\mu_{(w)}(x)=\mathbb{E}[Y_i(w)|X_i=x],\sigma^2_w(x)=Var[Y_i(w)|X_i=x]$。

回归估计：对 $\mu_{(0)}(x),\mu_{(1)}(x)$  做非参回归得到 $\hat{\tau}_{REG}=\frac{1}{n}\sum_i (\hat{\mu}_{(1)}(X_i)-\hat{\mu}_{(0)}(X_i))$。

Augmented（混合） IPW：先做回归估计，再用IPW修正。
$$
\hat{\tau}_{AIPW}=\frac{1}{n}\sum_{i=1}^n\left(\hat{\mu}_{(1)}(X_i)-\hat{\mu}_{(0)}(X_i)+W_i\frac{Y_i-\hat{\mu}_{(1)}(X_i)}{\hat{e}(X_i)}-(1-W_i)\frac{Y_i-\hat{\mu}_{(0)}(X_i)}{1-\hat{e}(X_i)}\right)
$$
优美性质是，$\hat{\mu}$ 和 $\hat{e}$ 有一个是consistent的，那这就是 $\tau$ 的consistent estimator（weak robust），但是这个性质不给统计学家关心的收敛速度和confidence interval。

Oracle augmented IPW：
$$
\hat{\tau}_{AIPW}^*=\frac{1}{n}\sum_{i=1}^n\left(\mu_{(1)}(X_i)-\mu_{(0)}(X_i)+W_i\frac{Y_i-\mu_{(1)}(X_i)}{e(X_i)}-(1-W_i)\frac{Y_i-\mu_{(0)}(X_i)}{1-e(X_i)}\right)
$$
strong robust：
$$
|\hat{\tau}_{AIPW}-\hat{\tau}_{AIPW}^*|=\mathcal{O}_p\left(\max_{w}\mathbb{E}[(\hat{\mu}_{(w)}(X_i)-\mu_{(w)}(X_i))^2]^{1/2}\mathbb{E}[(\hat{e}(X_i)-e(X_i))^2]^{1/2} \right)
$$
也就是 $\mu,e$ 有一个估计的特别好，那么就可以弥补另一个估计。

Cross-fitting AIPW：防止过拟合的技术。将数据分组（一般比较少），$k:\{1,...,n\}\to\text{组别}$。
$$
\hat{\tau}_{AIPW}=\frac{1}{n}\sum_{i=1}^n\left(\hat{\mu}_{(1)}^{-k(i)}(X_i)-\hat{\mu}_{(0)}^{-k(i)}(X_i)+W_i\frac{Y_i-\hat{\mu}^{-k(i)}_{(1)}(X_i)}{\hat{e}^{-k(i)}(X_i)}-(1-W_i)\frac{Y_i-\hat{\mu}_{(0)}^{-k(i)}(X_i)}{1-\hat{e}^{-k(i)}(X_i)}\right)
$$
$-k(i)$：刨去第 $k(i)$ 组所有的数据学习。（然后把 $k_i$ 组的数据带进去用，有效防止过拟合）。

在三个基本假设下可得 $\sqrt{n}(\hat{\tau}_{AIPW}-\hat{\tau}_{AIPW}^*)\to_p 0$，并且我们可以给出渐进意义下 $\hat{\tau}_{AIPW}$ 的方差，使得我们可以给出asymptotic confidence interval。