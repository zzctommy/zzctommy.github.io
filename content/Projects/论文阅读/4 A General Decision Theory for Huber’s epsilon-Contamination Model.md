---
title: 4 A General Decision Theory for Huber’s epsilon-Contamination Model
date: 2026-11-20
---

[论文地址](https://arxiv.org/abs/1511.04144)

---

## 1 Introduction

Huber $\epsilon$-contamination model
$$
\mathbb{P}_{(\epsilon,\theta,Q)}=P=(1-\epsilon)P_\theta+\epsilon Q
$$
Want to investigate the quantity
$$
\inf_{\hat{\theta}}\sup_{\theta\in\Theta, Q}\mathbb{E}_{\epsilon,\theta, Q}L(\hat{\theta},\theta)
$$

## 2 Robust Testing

Given iid observations $X_1,\cdots,X_n\sim P$, test
$$
H_0:\quad P\in \{(1-\epsilon)P_0+\epsilon Q:Q\}
$$

$$
H_1:\quad P\in\{(1-\epsilon)P_1+\epsilon Q:Q\}
$$

If $\epsilon$ is known, Huber proved the optimal estimator is the clipped LRT(c, C depends on $\epsilon$).
$$
\phi_{Huber}=I\left(\prod_{i=1}^n (\frac{dP_1}{dP_0}(X_i)\lor c)\land C>t\right)
$$
The paper provides an adaptive estimator.
$$
\phi=I(|P_n(A)-P_0(A)|>|P_n(A)-P_1(A)|)
$$
where $P_n(A)=\frac{1}{n}\sum_{i=1}^n I(X_i\in A)$ be the empirical distribution and $A$ is the set that achieves TV distance: $A=\arg\max_{A}|P_0(A)-P_1(A)|$.

**Theorem2.1** Assume $TV(P_0,P_1)>2\epsilon$, then we have
$$
\sup_{P\in\{(1-\epsilon)P_0+\epsilon Q:Q\}}P\phi+\sup_{P\in\{(1-\epsilon)P_1+\epsilon Q:Q\}}P(1-\phi)\le 4\exp(-\frac{1}{2}n (TV(P_0,P_1)-2\epsilon)^2)
$$
大概证明思路：如果错判，比如真实情形是 $H_0$，错判为 $H_1$，则 $|P_n(A)-P_0(A)|\ge|P_n(A)-P_1(A)|$，结合 $|P(A)-P_0(A)|\le \epsilon$ 得到 $|P_n(A)-P(A)|\ge \frac{1}{2}TV(P_0,P_1)-\epsilon$. 然后直接放 Hoeffding $\mathbb{P}[|P_n(A)-P(A)|\ge t]\le 2e^{-2nt^2}$，然后合并两类错误。

## 3 Construction of Upper Bounds

设数据为 iid $X_1,\cdots,X_n\sim (1-\epsilon)P_\theta+\epsilon Q$. 假设parameter space $\Theta$ totally bounded：任意精度 $\eta>0$, 存在有限个球覆盖参数空间 $\Theta\subseteq \cup_{j=1}^N B_d(s_j,\eta)$.

定义 $m=\mathcal{M}(\delta,\Theta,TV)$ 为最小的正整数使得存在 $\theta_1,\cdots,\theta_m\subseteq \Theta$ 使得对任意 $\theta\in\Theta$，存在 $j\in\{1,\cdots,m\}$ 使得 $TV(P_{\theta},P_{\theta_j})\le\delta$. Say $\{\theta_j\}_{1\le j\le m}$ a $\delta$-covering set and $\mathcal{M}(\delta,\Theta,TV)$. 

For $j\not=k$ define the testing function
$$
\phi_{jk}=I(|P_n(A_{jk})-P_{\theta_j}(A_{jk})|>|P_n(A_{jk})-P_{\theta_k}(A_{jk})|)
$$
where $A_{jk}=\{p_{\theta_j}>p_{\theta_k} \}$ is the set achieves $TV(P_{\theta_j},P_{\theta_k})$ and $P_n=\frac{1}{n}\sum_{i=1}^n I(X_i\in A_{jk})$ be the empirical probability.

Remark: the author's notation is
$$
\phi_{jk}=I\left(|\frac{1}{n}\sum_{i=1}^n I(p_{\theta_j}(X_i)>p_{\theta_k}(X_i))-P_{\theta_j}(p_{\theta_j}(X)>p_{\theta_k}(X))|\\ \quad\quad\quad\quad>|\frac{1}{n}\sum_{i=1}^n I(p_{\theta_j}(X_i)>p_{\theta_k}(X_i))-P_{\theta_k}(p_{\theta_j}(X)>p_{\theta_k}(X))| \right)
$$
where $p_{\theta}=dP_{\theta}/d\mu$ is the density function w.r.t. a common measure.

Estimator $\hat{\theta}$ is defined to be $\theta_{\hat{j}}$ where 

$$
\hat{j}=\arg\min_{1\le j\le m}\sum_{k\not=j}\phi_{jk}
$$

**Theorem 3.1** Assume $\eta>8(\epsilon+\delta)$. For the estimator $\hat{\theta}$ define above, we have
$$
\sup_{\theta\in\Theta, Q}\mathbb{P}_{(\epsilon, \theta, Q)} \{TV(P_{\hat{\theta}},P_{\theta})>\eta+\delta \}\le4\mathcal{M}^2(\delta,\Theta,TV)\exp(-\frac{1}{2}n(\eta/4-2(\epsilon+\delta))^2)
$$
大概证明思路：如果估计的不好，$TV(P_{\hat{\theta}},P_{\theta})>\eta+\delta$。由covering的性质，取 $\theta_{j^*}$ s.t. $TV(P_{\theta_{j^*}},P_{\theta})\le\delta$，就有 $TV(P_{\theta_{j^*}},P_{\hat{\theta}})>\eta$ ，这会要求某个 $P_n(A_{jk})$ 偏离期望很多。theorem2.1给出这种情形的概率，呈指数下降，全部控制就得到这个定理中的界。注意条件 $\eta>8(\epsilon+\delta)$ 目的是让 $\exp$ 里的项为正，希望 $\eta$ 是主项。 

**Corollary** Take（统计误差和污染误差取max。统计误差的取法是想让 $n\eta^2\asymp\log\mathcal{M}$）
$$
\eta\asymp \left\{\frac{\log\mathcal{M}(\delta,\Theta,TV)}{n}+\delta^2 \right\}\lor\epsilon^2
$$
Then we have the convergence rate of our risk with loss function be $TV^2$. The rate is(need some calculation)
$$
\min_{\delta>0}\left\{\frac{\log\mathcal{M}(\delta,\Theta,TV)}{n}+\delta^2 \right\}\lor\epsilon^2
$$
**Theorem3.2(Chen,Gao&Ren 2015)** $L(\cdot,\cdot)$ is a loss function defined on parameter space $\Theta$.  Define
$$
\omega(\epsilon,\Theta) =\sup\{L(\theta_1,\theta_2):TV(P_{\theta_1},P_{\theta_2})\le\frac{\epsilon}{1-\epsilon};\theta_1,\theta_2\in\Theta \}
$$
Suppose $\mathcal{R}(0)$ s.t. 
$$
\inf_{\hat{\theta}}\sup_{\theta\in\Theta, Q}\mathbb{P}_{(\epsilon,\theta,Q)}\{L(\hat{\theta},\theta)\ge\mathcal{R}(\epsilon) \}\ge c
$$
holds for $\epsilon=0$, then it holds for $\mathcal{R}(\epsilon)=\mathcal{R}(0)\lor\omega(\epsilon,\Theta)$.

For loss $L=TV$, $\omega(\epsilon,\Theta)\asymp \epsilon$. 在本文的loss $L=TV^2$ 就给出污染下界 $\epsilon^2$，因此污染下界不能改进。

另一个引用的结果给出 $\mathcal{R}(0)\asymp \min_{\delta>0}\left\{\frac{\log\mathcal{M}(\delta,\Theta,TV)}{n}+\delta^2 \right\}$ ，因此统计误差下界不能改进，作者的estimator已经达到最优收敛速率。

作者指出在高维问题中theorem3.1未必很精细，对他的定理3.1给出一个优化：

Let $\Theta'=\{\theta_1,\cdots,\theta_m\}$ be a $\delta$-covering set of $\Theta$. For any integer $l$ define 
$$
D_l(\delta)=\max_{\theta_0\in\Theta}|\{\theta\in\Theta':l\delta< TV(P_{\theta},P_{\theta_0})\le(l+1)\delta\}|
$$
也就是每一层最多有几个 $\Theta'$ 的元素。

**Theorem3.3** Let $L$ s.t. $\frac{L}{4}\delta-2(\epsilon+\delta)>0$ and $\frac{L}{4}$ is an integer. For the estimator $\hat{\theta}$ we have
$$
\sup_{\theta\in\Theta,Q}\mathbb{P}_{(\epsilon,\theta,Q)}\{TV(P_{\hat{\theta}},P_{\theta})>(l+1)\delta \}\\
\le2\sum_{l\ge L/4} D_l(\delta)\exp(-\frac{1}{2}n(l\delta-2(\epsilon+\delta))^2\\
+2\left(\sum_{l=0}^{L/4-1}D_l(\delta) \right)\sum_{l\ge L}D_l(\delta)\exp(-\frac{1}{2}n((l-3L/4)\delta-2(\epsilon+\delta))^2)
$$

## 4 Applications

Apply the theorems to three problems: density estimation with Holder smoothness, sparse linear regression and low-rank trace regression

### 4.1 Density Estimation

Consider iid observation $X_1,\cdots,X_n\sim (1-\epsilon)P_f+\epsilon Q$ where $f=\frac{dP}{d\lambda}$ be the density w.r.t. Lebesgue measure supported on $[0,1]$.

Construct the standard wavelet basis. Define the density class
$$
\mathcal{H}_{den}(\beta,M)=\left\{f=\sum_{l\ge 0,0\le k\le 2^l-1} f_{lk}\psi_{lk}:f\ge 0, \int_{0}^1 f=1,\sup_{l\ge 0, 0\le k\le 2^l-1}2^{l(1/2+\beta)}|f_{lk}|\le M \right\}
$$
引用他人结果有
$$
\log \mathcal{M}(\delta,\mathcal{H}_{den}(\beta,M),TV)\asymp \delta^{-1/\beta}
$$
因此统计误差部分
$$
\frac{\log \mathcal{M}(\delta,\mathcal{H}_{den}(\beta,M),TV)}{n}+\delta^2\asymp \frac{\delta^{-1/\beta}}{n}+\delta^2
$$
为了平衡，取 $\delta=n^{-\beta/(2\beta+1)}$。

**Corollary 4.1** For the Holder class $\mathcal{H}_{den}(\beta,M)$, there exists constants $C,C'$ s.t. $||\hat{f}-f||_1^2\le C(n^{-\frac{2\beta}{2\beta+1}}\lor \epsilon^2)$  with $\mathbb{P}_{(\epsilon,f,Q)}$-probability  at least $1-\exp(-C(n^{\frac{1}{2\beta+1}}+n\epsilon^2))$ uniformly over $Q$ and $f\in\mathcal{H}_{den}(\beta,M)$.

这与经典的非参回归下界吻合，因此最优。

### 4.2 Sparse Linear Regression

Linear regression model 
$$
y_i=X_i^T\theta+z_i
$$
without contamination: $X_i\sim N(0,\Sigma)$ and $z_i\sim N(0,\sigma^2)$ are independent.

Now we add contamination, both the design and response can be contaminated, that is, $(X_1,y_1),\cdots,(X_n,y_n)\sim \mathbb{P}_{(\epsilon,\theta,Q)}=(1-\epsilon)P_\theta+\epsilon Q$.

我都没学过无污染的，看这个是不是没有意义啊。。。

## 5 Discussion with an Example under Supreme Norm

Consider the model: (denoted by $P_f$)
$$
dY_t=f(t)dt+\frac{1}{\sqrt{n}}dW_t\quad t\in[0,1]
$$
We observe iid processes $dY_{t,i}=f(t)dt+dW_{t,i} \sim (1-\epsilon)P_f+\epsilon Q$.

We use a slightly modified Holder class
$$
\mathcal{H}(\beta,M)=\left\{f=\sum_{l\ge 0, 0\le k\le 2^l-1} f_{lk}\psi_{lk}: \sup_{l\ge 0, 0\le k\le 2^l-1}2^{l(1/2+\beta)}|f_{lk}|\le M \right\}
$$
Target: use supereme loss $||\hat{f}-f||_{\infty}$, want to achieve the optimal rate.

 Let $L$ be the largest integer s.t. $2^L\le (\frac{\log n}{n}\lor\epsilon^2)^{-\frac{1}{2\beta+1}}$. 

Estimator is
$$
\hat{f}=\sum_{0\le l\le L,0\le k\le 2^l-1}\hat{f}_{lk}\psi_{lk}
$$
for $\hat{f}_{lk}=Median(\{y_{lk,i}\}_{i=1}^n)$, $y_{lk,i}=\int_0^1\psi_{lk}(t)dY_{t,i}$ are empitical wavelet coefficients.

**Theorem 5.1** Assume $\epsilon<\frac{1}{4}$, for Holder class $\mathcal{H}(\beta,M)$, there are constants $C,C'$ s.t.
$$
||\hat{f}-f||_{\infty}^2\le C\left[\left(\frac{n}{\log n}\right)^{-\frac{2\beta}{2\beta+1}}\lor \epsilon^{\frac{4\beta}{2\beta+1}} \right]
$$
with $\mathbb{P}_{\epsilon,f,Q}$-probability at least $1-\exp(-C'(\log n +n\epsilon^2))$ uniformly over $f\in\mathcal{H}(\beta,M)$ and $Q$.

The rate is optimal by theorem 3.2.









