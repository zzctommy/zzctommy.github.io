---
title: Statistics读书进度
---

# 进度

起到一个鞭策自己看书的作用，每天读一点记下来也许有成就感。


| Date        | Nonparametric Regression(190) | Casual Inferene(118) |
| ----------- | ----------- | --------------|
| ?-26.10.19  |15 | 2 |
| 26.10.20    |19 | 12 |
| 26.10.21    |30 | 18 |
| 26.10.22    |30 | 28 |
| 26.10.23    |30 | 33 |
| 26.10.24    |46 | 36 |
| 26.10.25-28    |在看论文、写作业|一点没看www|
| 26.10.29    |61 | 36 |
| 26.10.24    |76 | 36 |

------

# 每日总结

我意识到应该每节学完之后写一个概括，否则我会迷失在技术细节里。

## Nonparametric Regression

### 1.1

非参（无限维模型）统计的典例，例如从函数空间recover probabilistic density。

三个例子：

+ 有一些取样 $X_1,\cdots,X_n$, recover density $p$.
+ 有一些数据 $Y_i=f(X_i)+\xi_i$，$\xi_i$ 是噪声，recover $f$。
+ SDE $dY(t)=f(t)dt+\frac{1}{\sqrt{n}}dW(t)$，recover $f$.

### 1.2

引入Kernel density estimator $\hat{p_n}=\frac{1}{nh}\sum_i K((X_i-x)/h)$，function class用于限制bias，variance比较好控制，主要结果有：点态MSE $\sup\limits_{x_0\in \mathbb{R}}\sup\limits_{p\in \mathcal{P}(\beta,L)}\mathbb{E}_p[(\hat{p_n}(x_0)-p(x_0))^2]\le C n^{-2\beta/(2\beta+1)}$ 。

### 1.3

通过Fourier transform给出更精确的MISE估计，因为求导在Fourier transform下的简单表示主要考虑Sobolev class，主要结果是：$L^2$ 且对称的kernel，如果 $ess\sup_{t\not=0}\frac{|1-\hat{K}(t)|}{|t|^{\beta}}\le A$ ，那么MISE的上界 $\sup_{p\in \mathcal{P_s}(\beta,L)}\mathbb{E}_p[\int_{\mathbb{R}}(\hat{p_n}(x)-p(x))^2dx]\le C n^{-2\beta/(2\beta+1)}$ 。

### 1.4

Find optimal $h$ for kernel estimator. cross-validation首次出现（看了因果推断才懂，不过这里只是构造estimator而不是减少过拟合问题）。目标是寻找恰当的h给出好的kernel estimator。因为 MISE(h)依赖未知density p，这一节给出完全依赖data的unbiased estimator CV(h) of $MISE(h)-\int_\mathbb{R}p^2(x)dx$，所以问题变成最小化CV(h)：
$$
CV(h)=\int \hat{p_h}(x)^2dx-\sum_{i=1}^n\frac{2}{n}\hat{p}_{h,-i}(X_i)
$$

另外还有从Fourier transform角度出发构造的另一个 $MISE(h)-\int_\mathbb{R}p^2(x)dx$ 的无偏统计量：
$$
\tilde{J}(h):=\int_{\mathbb{R}} -2\hat{K}(hw)+\hat{K}(hw)^2(1-\frac{1}{n})|\phi_n(w)|^2dw+\frac{2}{n}\int_{\mathbb{R}}\hat{K}(hw)dw\\
\mathbb{E}_p[\hat{J}(h)]=2\pi(1-\frac{1}{n})[MISE(h)-\int_\mathbb{R}p(x)^2dx]
$$


### 1.5

两类design。random design：$X_i$ are random variables。fix design：$X_i$  预先给定，例如 $X_i=\frac{i}{n}$ 。

目标是估计回归函数 $f(x)=\mathbb{E}[Y|X=x]$，构造了Nadaraya-Watson estimator $f_n^{NW}(x)=\frac{\sum_i Y_iK((X_i-x)/h)}{\sum_i K((X_i-x)/h)}$(for random design)。fix design没太看懂，没细讲，遇到再说。

### 1.6

通过如下方式定义 local polynomial estimator of order $l$ （简记为 $LP(l)$ estimator）：kernel $K:\mathbb{R}\to\mathbb{R}$，$h>0$ bandwidth，$l\in \mathbb{Z}_{\ge 0}$，$U(u)=(1,u,u^2/2!,\cdots,u^l/l!)^T$。
$$
\hat{\theta}_n(x)=\arg\min_{\theta\in \mathbb{R}^{n+1}} \sum_{i=1}^n (Y_i-\theta^TU(\frac{X_i-x}{h}))^2K(\frac{X_i-x}{h})\in \mathbb{R}^{l+1}\\
\hat{f_n}(x)=U^T(0)\hat{\theta}_n(x)
$$
这么写的动机是，在局部想象 $f$ 是 $l$ 阶 Taolor 多项式，对系数做线性回归。乘以 $K$ 保证只用到局部的数据。但是为什么这真的可能是一个estimator，可以证明：如果f局部上是 $l$ 次多项式，并且线性回归矩阵  $B_{nx}$ 可逆，那么 $f(x)=\hat{f_n}(x)$。特别的，$LP(0)$ estimator 就是 Nadatraya-Watson estimator.

以下研究这个estimator的性质。

记 $B_{nx}=\frac{1}{nh}\sum_{i=1}^n U(\frac{X_i-x}{h})U^T(\frac{X_i-x}{h})K(\frac{X_i-x}{h}),a_{nx}=\frac{1}{nh}\sum_{i=1}^n Y_iU(\frac{X_i-x}{h})K(\frac{X_i-x}{h})$。

那么 $\hat{\theta}_n(x)=\arg\min_{\theta\in\mathbb{R}^{l+1}} -2\theta^T a_{nx}+\theta^T B_{nx}\theta$，可以得到正规方程 $B_{nx}\hat{\theta}_n(x)=a_{nx}$。

我们关心的情况：$K\ge 0$，保证 $B_{nx}$ 半正定；进一步假设 $B_{nx}$ 正定，可以得到唯一解 $\hat{\theta}_n=B_{nx}^{-1}a_{nx}$。

在此情况下 $\hat{f}_n$ 是 $f$ 的linear estimator：$\hat{f_n}(x)=\sum_{i=1}^n Y_i W_{ni}^*(x)$，其中
$$
W_{ni}^*(x)=\frac{1}{nh}U^T(0)B_{nx}^{-1}U(\frac{X_i-x}{h})K(\frac{X_i-x}{h})
$$
首先估计点态MSE的收敛速率
$$
MSE(x_0):=\mathbb{E}_f[(\hat{f_n}(x_0)-f(x_0))^2]=b(x_0)^2+\sigma^2(x_0)\\
\text{where}\ b(x_0)=\mathbb{E}_f[\hat{f_n}(x_0)]-f(x_0),\sigma^2(x_0)=Var_f[\hat{f_n}(x_0)]
$$
**定理（点态MSE一致趋于0的速率）**：$(f:[0,1]\to \mathbb{R})\in \Sigma(\beta,L)$（Holder class），$\hat{f_n}$ be the $LP(l)$ estimator for $f$ with $l=\lfloor \beta\rfloor$constructed from observations $(X_i,Y_i)$ s.t. $Y_i=f(X_i)+\xi_i$. 一些条件下

Choose $\psi_n=n^{-\frac{\beta}{2\beta+1}}$, we have
$$
\limsup_{n\to\infty} \sup_{f\in \Sigma(\beta,L)}\sup_{x_0\in [0,1]}\mathbb{E}_f[\psi_n^{-2}|\hat{f_n}(x_0)-f(x_0)|^2]\le C<\infty\\
\limsup_{n\to\infty} \sup_{f\in \Sigma(\beta,L)}\mathbb{E}_f[\psi_n^{-2}||\hat{f_n}(x_0)-f(x_0)||_{L^2}^2]\le C<\infty
$$
在 fixed design $X_i=\frac{i}{n}$ 时条件可以有简单刻画，只需要 $K$ 满足一些容易的条件。

我们还可以研究 $L^{\infty}$-risk of $\hat{f_n}$ ：$\mathbb{E}_f[||\hat{f_n}-f||_{L^\infty}^2]$ 的收敛速度。

**定理（$L^{\infty}$-risk 一致收敛速率）**：$(f:[0,1]\to \mathbb{R})\in \Sigma(\beta,L)$（Holder class），$\hat{f_n}$ be the $LP(l)$ estimator for $f$ with $l=\lfloor \beta\rfloor$constructed from observations $(X_i,Y_i)$ s.t. $Y_i=f(X_i)+\xi_i$. 一些条件下

Choose $\psi_n=(\frac{\log n}{n})^{\frac{\beta}{2\beta+1}}$ we have
$$
\limsup_{n\to\infty} \sup_{f\in \Sigma(\beta,L)}\mathbb{E}_f[\psi_n^{-2}||\hat{f_n}(x_0)-f(x_0)||_{L^\infty}^2]\le C<\infty
$$

### 1.7

For nonparametric regression model $Y_i=f(X_i)+\xi_i$, $X_i$ 在 $[0,1]$ 均布（否则estimate很差），取 $L^2([0,1])$ 的 ONB $\{\varphi_j\}$ 并计算傅里叶系数 $\theta_j=\int_0^1 f\varphi_j$，得到 $\hat{\theta}_j=\frac{1}{n}\sum_{i=1}^n Y_i\varphi_j(X_i)$，定义projection estimator $\hat{f}_{nN}=\sum_{j=1}^N \hat{\theta}_j \varphi_j(x)$ （估计前 $N$ 项Fourier系数）

For estimator of probability density：$X_1,\cdots,X_n$ iid of density $p$ supported on an interval $A\subseteq \mathbb{R}$, ONB $\{\varphi_j\}$ on $L^2(A)$. $c_j=\int p\varphi_j, \hat{c}_j=\frac{1}{n}\sum_{i=1}^n\varphi_j(X_i)$. Define $\hat{p}_{nN}(x)=\sum_{j=1}^N \hat{c}_j \varphi_j(x)$. 

常用两个ONB：

+ 三角级数 $\varphi_1(x)=1,\varphi_{2k}(x)=\sqrt{2}\cos(2\pi kx), \varphi_{2k+1}(x)=\sqrt{2}\sin(2\pi kx), k\ge 1$. ONB for $L^2([0,1])$.

+ wavelet base : $\psi$ 有较好性质，可以证明 $\psi_{jk}=2^{j/2}\psi(2^j x-k),\quad (j,k\in \mathbb{Z})$ 是 $L^2(\mathbb{R})$ 的 ONB。

常用函数类：

+ $W(\beta, L)=\{f:[0,1]\to\mathbb{R}\big|f^{(\beta-1)} A.C.,||f^{(\beta)}||_{L^2}\le L\}$.
+ $W^{per}(\beta,L)=\{f\in W(\beta,L):f^{(j)}(0)=f^{(j)}(1),j=0,1,\cdots,\beta-1 \}$.
+ Sobolev ellipsoid: $a_j:=\begin{cases}j^\beta, j\ even\\ (j-1)^\beta,j\ odd \end{cases}$, $\{\varphi_j\}$ trigonometric ONB. $\Theta(\beta,Q)=\{\theta\in l^2(\mathbb{N}):\sum_{j=1}^na_j^2\theta_j^2<\infty \}$. Then $W^{per}(\beta,L)\simeq\Theta(\beta,\frac{L^2}{\pi^{2\beta}})$ via $f=\sum \theta_j \varphi_j$.
+ $\tilde{W}(\beta,L)=\{f\in L^2([0,1]):\theta\in\Theta(\beta,Q=L^2/\pi^{2\beta}) \}$（$W^{per}$ 在 $\beta$ 不是整数时的推广）

在特殊情况 $X_i=\frac{i}{n}$，三角ONB以及其他一些条件下，可以拿到 MISE的收敛速度：

**定理** For  $N=\lfloor \alpha n^{\frac{1}{2\beta+1}}\rfloor, \beta\ge 1, L>0$, the projection estimator $\hat{f}_{nN}$ satisfies
$$
\limsup_{n\to \infty}\sup_{f\in \tilde{W}(\beta,L)} \mathbb{E}_f[n^{\frac{2\beta}{2\beta+1}}||\hat{f}_{nN}-f||_{L^2}^2]\le C=C(\beta,L,\alpha)
$$
三个Generalization：

+ least square estimator. $\hat{\theta}^{LS}=\arg\min_{\theta\in\mathbb{R}^N}\sum_{i=1}^N(Y_i-\theta^T\varphi(X_i))^2$ where $\varphi=(\varphi_1,\cdots,\varphi_N)$ be ONB. 在上面讨论的情况下这是 $\hat{f_{nN}}$ 的另一种表述。这个线性回归在比较好的情况下可以计算正规方程求解。
+ weighted projection estimator. $\hat{\theta_j}=\frac{1}{n}\sum_{i=1}^n Y_i\varphi_j(X_i)$ and $f_{n,\lambda}(x)=\sum_{j=1}^\infty \lambda_j\hat{\theta}_j\varphi_j(x)$. 之前的 $\hat{f}_{nN}$ 就是 $\lambda_j=I(j\le N)$ 的结果。取更好的 $\lambda$ 会给出更好的MISE收敛。
+ penalized least square estimator. A popular estimator is Tikhonov regularization estimator. $\hat{\theta}^{TR}=\arg\min_{\theta\in\mathbb{R}^N}\{\sum_{i=1}^N(Y_i-\theta^T\varphi(X_i))^2+\sum_{j=1}^N b_j\theta_j^2\}$. Lasso estimator $\hat{\theta}^{TR}=\arg\min_{\theta\in\mathbb{R}^N}\{\sum_{i=1}^N(Y_i-\theta^T\varphi(X_i))^2+\sum_{j=1}^N b_j|\theta_j|\}$

Lasso的好处多多，例如数值复杂度只有 $O(Nn)$，例如可以把某些权重学习成 $0$.

### 1.8

一般的oracle定义。暂时没用，不抄了。

### 1.9

Consider regression model $Y_i=f(X_i)+\xi_i, X_i\in[0,1],\mathbb{E}[\xi_i]=0$.

$\{f_{\tau}:\tau\in \mathcal{T} \}$ be a family of estimator based only on sample $\{(X_i,Y_i)\}$ and $\tau$.

discretized MISE $r_{n,\tau}^D(f)=\mathbb{E}_f[||f_{\tau}-f||_{2,n}^2]$ where $||f_\tau-f||_{2,n}=(\frac{1}{n}\sum_{i=1}^n (f_{\tau}(X_i)-f(X_i))^2)^{1/2}$.

minimizing $r_{n,\tau}^D(f)$ is equivalent to minimizing $J(\tau)=\mathbb{E}_{f}[||f_{\tau}||_{2,n}^2-\frac{2}{n}\sum_{i=1}^n f_{\tau}(X_i)f(X_i)]$.

We find an unbiased estimator $\hat{J}(\tau)=||f_{\tau}||_{2,n}^2-\frac{2}{n}\sum_{i=1}^n Y_i f_{\tau}(X_i)+\text{residual of }\xi$.

minimizing $\hat{J}(\tau)$ to give $\hat{\tau}$, equivalently to minimize $C_p(\tau)=\frac{1}{n}\sum_{i=1}^n (Y_i-f_{\tau}(X_i))^2+\text{residual}$

In the end consider the cross-validation technique. $CV_*(\tau)=\frac{1}{n}\sum_{i=1}^n (Y_i-f_{\tau,-i}(X_i))^2$.

$CV_*$ may not be unbiased, but still have good properties, can choose nearly optimal parameter.

### 1.10

非参回归、Gaussian white noise model, Gaussian sequence model 之间的相互转化。没有严格的推导，感觉和因果推断关系也不大，略看了。





## Casual Inference

### Lec1

ATE（average treatment effect） $\tau=\mathbb{E}[Y_i(1)-Y_i(0)]$。

$Y_i(W_i)$：outcome，是否康复，$W_i$ 是否接受治疗。

理想化情况 $\tau_{SATE}=\frac{1}{n}\sum_iY_i(1)-Y_i(0)$ （sample average treatment effect）就是一个估计，但是在 $\{Y_i(0),Y_i(1)\}$ 中只能选一个观测。

Difference-in-means estimator $\hat{\tau}_{DM}=\frac{1}{n_1}\sum_{W_i=1}Y_i-\frac{1}{n_0}\sum_{W_i=0}Y_i$。

后面拿linear model和nonlinear model的OLS方法（nonlinear用Huber-White）算了example，OLS给出的estimator方差都小于difference-in-means estimator，所以具体问题中DM往往不是最优的

### Lec2

整体比较会丢失因果，因此通过加入协变量 $X_i$ 给被试者分组，探测协变量带来的因果。

group-wise ATE $\tau(x)=\mathbb{E}[Y_i(1)-Y_i(0)|X_i=x]$。

Aggregating difference-in-means estimator: $\hat{\tau}_{AGG}=\sum_{x\in \chi}\frac{n_x}{n}\hat{\tau}(x)$ where $\hat{\tau}(x)=\frac{1}{n_{x,1}}\sum\limits_{X_i=x,W_i=1}Y_i-\frac{1}{n_{x,0}}\sum\limits_{X_i=x,W_i=0}Y_i$，这个estimator渐进正态，要求协变量离散（要求每个协变量被重复观测的次数够多才精确）。

unconfoundedness：$\{Y_i(0),Y_i(1)\}\perp W_i\big|X_i$，也就是在协变量分组下，潜在结果与是否被治疗无关。

propensity score：$e(x)=\mathbb{P}[W_i=1\big|X_i=x]$。

在unconfoundedness成立的条件下，$\{Y_i(0),Y_i(1)\}\perp W_i\big|e(X_i)$。

Propensity stratification：按照propensity score的大小分为若干组，每一组跑DM estimate，通过平衡组的个数（组越多对每一组propensity score差别越小，DM越准；组太多每一组数据太少，样本不足带来大误差）来构造estimator。

IPW（inverse propensity weighting）$\hat{\tau}_{IPW}=\frac{1}{n}\sum_i \frac{W_iY_i}{\hat{e}(X_i)}-\frac{(1-W_i)Y_i}{1-\hat{e}(X_i)}$ （$\hat{e}$ 用非参回归计算）

oracle IPW：$\hat{\tau}_{IPW}^*=\frac{1}{n}\sum_i \frac{W_iY_i}{e(X_i)}-\frac{(1-W_i)Y_i}{1-e(X_i)}$

在unconfoundedness下，$\hat{\tau}^*_{IPW}$ 是ATE的unbiased estiamtor，$\frac{1}{\sqrt{n}}$ rate逼近 $\tau$。

协变量离散的时候 $\hat{\tau}_{AGG}$ 比 $\hat{\tau}_{IPW}^*$ 精确（方差更小，因为AGG使用的是精确的每组人数，IPW使用的是估计每组人数）。但是由于IPW更“抽象”，只需要考虑propensity score的性质而不需要准确的每组人数，所以更加可以研究。

### Lec3

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
|\hat{\tau}_{AIPW}-\hat{\tau}_{AIPW}^*|=\mathcal{O}_p\left(\max_{w}\mathbb{E}[(\hat{\mu}_{(w)}(X_i)-\mu_{(w)}(X_i))^2]^{1/2}\mathbb{E}[(\hat{e}(X_i)-e(X_i))^2]^{1/2} \right)+\color{red}{o_p(n^{-1/2})}
$$
（红色：AI指出讲义应该加上余项）也就是 $\mu,e$ 有一个估计的特别好，那么就可以弥补另一个估计。

Cross-fitting AIPW：防止过拟合的技术。将数据分组（一般比较少），$k:\{1,...,n\}\to\text{组别}$。
$$
\hat{\tau}_{AIPW}=\frac{1}{n}\sum_{i=1}^n\left(\hat{\mu}_{(1)}^{-k(i)}(X_i)-\hat{\mu}_{(0)}^{-k(i)}(X_i)+W_i\frac{Y_i-\hat{\mu}^{-k(i)}_{(1)}(X_i)}{\hat{e}^{-k(i)}(X_i)}-(1-W_i)\frac{Y_i-\hat{\mu}_{(0)}^{-k(i)}(X_i)}{1-\hat{e}^{-k(i)}(X_i)}\right)
$$
$-k(i)$：刨去第 $k(i)$ 组所有的数据学习。（然后把 $k(i)$ 组的数据带进去用，有效防止过拟合）。

在三个基本假设下可得 $\sqrt{n}(\hat{\tau}_{AIPW}-\hat{\tau}_{AIPW}^*)\to_p 0$，并且我们可以给出渐进意义下 $\hat{\tau}_{AIPW}$ 的方差，使得我们可以给出asymptotic confidence interval。

### Lec4

Statistics setting：observed iid $(X_i,Y_i,W_i)$，$X_i\in \chi$ personal features.

CATE: $\tau(x)=\mathbb{E}[Y_i(1)-Y_i(0)|X_i=x]$，i.e. ATE over a target group.

现在人们想要估计 $\tau(x)$，这有明确的动机，例如我们想知道某个药对不同人群的副作用程度。

一个直接方法是非参回归估计 $\mu_{(w)}(x)=\mathbb{E}[Y_i|X_i=x,W_i=w]$ 得到 $\hat{\mu}_{(w)}(x)$，得到回归估计 $\hat{\tau}_T(x)=\hat{\mu}_{(1)}(x)-\hat{\mu}_{(0)}(x)$。

这个估计有两个问题，一是（在观察性研究中）如果对照组和实验组人数差距过大，会有一个estimator相比另一个estimator非常不准，从而会导致误判因果（例如真实值是0.8，0.2，学得好的直接学出了0.2，学得差的学出了一个0.8附近波动的函数，于是在波峰波谷等位置人们会认为有因果，但实际没有）。二是，如果 $W_i=0,1$ 的数据的组别差异很大，也就是集中在不同的 $X_i$，那么在每一个组的对照组实验组人数差距都很大都会出现上面引起的因果假象。

Semiparametric modeling：$Y_i(w)=f(X_i)+w\tau(X_i)+\epsilon_i(w)$，where $\tau(x)=\psi(x)\cdot \beta, \beta\in \mathbb{R}^k$。也就是人们已经有了一些 $\tau(x)$ 的候选模型，认为这部分是线性的。

Robinson研究过这类问题。在unconfoundedness下可以重写模型为
$$
Y_i-m(X_i)=(W_i-e(X_i))\psi(X_i)\cdot \beta+\epsilon_i\\
\text{where}\  m(x)=\mathbb{E}[Y_i\big|X_i=x]=f(X_i)+e(X_i)\tau(X_i)
$$
这给出估计 $\beta$ 的oracle algorithm：$\tilde{Y}_i^*=Y_i-m(X_i),\tilde{Z}_i^*=(W-e(X_i))\psi(X_i)$，estimator
$$
\hat{\zeta}_R^*=arg\min_{b\in \mathbb{R}^k}(\sum_{i=1}^n\tilde{Y}_i^*-(\tilde{Z}_i^*)^Tb)^2
$$
可以证明 $\sqrt{n}(\hat{\zeta}_R^*-\beta)\Rightarrow \mathcal{N}(0,V_R)$ where $V_R=Var[\hat{Z}_i^*]^{-1}Var[\tilde{Z}_i^*\tilde{Y}_i^*]Var[\hat{Z}_i^*]^{-1}$。

实操的时候我们用cross-fitting：

1. 非参回归 $Y\sim X,W\sim X$ 得到 $\hat{m}(x), \hat{e}(x)$。

2. 定义 $\tilde{Y}_i=Y_i-\hat{m}^{(-k(i))}(X_i),\tilde{Z}_i=(W_i-\hat{e}^{-(k(i))}(X_i))\psi(X_i)$。

3. OLS回归 $\tilde{Y}_i\sim \tilde{Z}_i$ 得到 $\hat{\zeta}_R$。

可以证明如果m和e的非参回归的MSE都是 $o_p(n^{-1/2})$ 的，那么 $\sqrt{n}(\hat{\zeta}_R^*-\hat{\zeta}_R)\to_p 0$，此时 $\hat{\zeta}_R$ 也有上面的渐进正态。

但是在复杂的covariate环境下，我们更希望机器自己学习CATE而不是我们预设 $\psi$。

模型写成 $Y_i(w)=\mu_{(0)}(X_i)+w\tau(X_i)+\epsilon_i(w)$，在unconfoundedness下 $\mathbb{E}[\epsilon_i\big|X_i,W_i]=0$。

重写为 $Y_i-m(X_i)=(W_i-e(X_i))\tau(X_i)+\epsilon_i$。

注意到 $\tau(X_i)=arg\min_{\tau'}\mathbb{E}[((Y_i-m(X_i))-(W_i-e(X_i))\tau'(X_i))^2]$。

所以定义empirical loss function
$$
\hat{\tau}_R^*(\cdot)=arg\min_{\tau'}\frac{1}{n}\sum_{i=1}^n ((Y_i-m(X_i))-(W_i-e(X_i))\tau'(X_i))^2+\Lambda_n(\tau'(\cdot))
$$
其中 $\Lambda_n$ 是给候选函数复杂度增加的惩罚（Lasso的思想）。

实操的时候依旧cross-fitting
$$
\hat{\tau}_R(\cdot)=arg\min_{\tau} \hat{L}_n(\tau(\cdot))+\Lambda(\tau(\cdot))\\
\text{where}\ \hat{L}_n(\tau(\cdot))=\frac{1}{n}\sum_{i=1}^n[(Y_i-\hat{m}^{(-k(i))}(X_i))-(W_i-\hat{e}^{(-k(i))}(X_i))\tau(X_i)]^2
$$
$\hat{L}_n$ 被称为 R-loss。

Now we can say something on "treatment heterogeneity"

1. 比较候选模型准确性：在训练数据上得到若干候选模型（例如不同惩罚强度的lasso），在验证集（未参与训练）上对每个候选模型算 R-loss，认为loss越小越准确。
2. 检查模型是否具有区分heterogeneity的能力：先估计 $\hat{\tau}(x)$，按 $\hat{\tau}(x)$ 的大小分组（stratification），每一组用上一讲的robust方法计算 $\hat{\tau}_{AIPW}(x)$。例如，for simplicity只有H（high）和L（low）两组，其中 $G_H=\{x:\tau(x)>c\},G_L=\{x:\tau(x)\le c\}$，Hypothesis：$\tau_H-\tau_L>0$。在每一组内分别跑AIPW，可以算出 $\hat{\tau}_H-\hat{\tau}_L$ 的渐进正态得到CI，观察CI是否包含0得到结论。
3. 衡量模型准确性：训练出 $\hat{\tau}(x)$ 之后，考虑半参数模型 $\tau(x)\approx \alpha+\beta\hat{\tau}(x)$ ，run Robinsin method。
   + $\hat{\beta}\approx 1$ 说明 $\hat{\tau}$ nice
   + $\hat{\beta}$ 显著大于 $0$ ：$\hat{\tau}$ 含真实信号，不是“pure noise”
   + $0<\hat{\beta}<1$：$\hat{\tau}$ 需要缩小，训练得到的模型高估了差异
   + $\hat{\beta}>1$：$\hat{\tau}$ 需要放大，训练得到的模型低估了差异

