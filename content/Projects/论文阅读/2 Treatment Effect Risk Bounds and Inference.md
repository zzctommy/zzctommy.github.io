---
title: Treatment Effect Risk Bounds and Inference
date: 2026-11-04
---

$\newcommand{\CVaR}{\operatorname{CVaR}}$

$\newcommand{\Var}{\operatorname{Var}}$

$\newcommand{\E}{\mathbb{E}}$

-----

[论文地址](https://arxiv.org/abs/2201.05893)

---

### 1 Introduction

一个值得考虑的问题是，ATE是正的，并不代表对所有人都有利，我们也希望知道对效果最差的10%或者20%的人的治疗效果。

### 2 PROBLEM SET UP AND DEFINITIONS

potential outcome $Y^*(0),Y^*(1)\in\mathbb{R}$, covariates $X\in\chi$.

ITE(individual treatment effect) $\delta=Y^*(1)-Y^*(0)$.

ATE $\bar{\tau}=\mathbb{E}[Y^*(1)-Y^*(0)]$.

CATE $\tau(X)=\mathbb{E}[\delta|X]=\mu(X,1)-\mu(X,0)$ where $\mu(X,a)=\mathbb{E}[Y^*(a)|X]$.

Average effect among the most worst $(100\times\alpha) \%$：$\CVaR_{\alpha}(Z)=\sup_{\beta}(\beta+\frac{1}{\alpha}\mathbb{E}[(Z-\beta)_-])$ where $x_-=\min(x,0)$. The $\sup$ is attained by $\beta=F_Z^{-1}(\alpha)=\inf\{\beta:F_Z(\beta)\ge\alpha\}$ i.e. $\alpha$-quantile.

Provided by $F_Z$ continuous, then $\CVaR_{\alpha}(Z)=\mathbb{E}[Z|Z\le F_Z^{-1}(\alpha)]$.

Each individual is assigned with an assignment $A\in\{0,1\}$.

Observed data: $(X_i,A_i,Y_i)\sim(X,A,Y)$.

Assume unconfoundedness throughout $Y^*(a)\perp A|X$. Under unconf. ATE and CATE are identifiable: $\mu(X,a)=\mathbb{E}[Y|X,A=a]$. 

Propensity score: $e(X)=\mathbb{P}[A=1|X]$.

Marginal outcome $\bar{\mu}(X)=\mathbb{E}[Y|X]=e(X)\mu(X,1)+(1-e(X))\mu(X,0)$.

### 3 BOUNDS

#### 3.1 Upper Bound: The CATE-CVaR

**Theorem3.1** $\CVaR_\alpha(\delta)\le\CVaR_\alpha(\tau(X))$. 并且在某些情况下可以取等（只用这些信息不能改进）。

证明核心是，$\mathbb{E}[(\delta-\beta)_-|X]\le\min\{0,\mathbb{E}[\delta-\beta|X]=\tau(X)-\beta \}=(\tau(X)-\beta)_-$.

直观是 $\tau(X)$ 把治疗效果好的人也平均进来了，所以对最差的一批人的估计会偏乐观。

remark3：$\CVaR_{\alpha}(\tau(X))$ 本身也具备研究价值，作为治疗异质性的反映。

remark4：如果能很好的估计 $\CVaR_\alpha(\tau(X))$，那么我们也能够估计
$$
\mathbb{E}[\delta|F_{\tau(X)}^{-1}(\alpha)<\tau(X)\le F_{\tau(X)}^{-1}(\alpha')]=\frac{\alpha'\CVaR_{\alpha'}(\tau(X))-\alpha\CVaR_{\alpha}(\tau(X))}{\alpha'-\alpha}
$$
remark5：发现负面效果之后可以定位到是哪些特征：$\{X:\tau(X)<F_{\tau(X)}^{-1}(\alpha)\}$.

#### 3.2 Lower Bounds under Limited Residual Heterogeneity Range

**Theorem3.2** Suppose $|\tau(X)-\delta|\le b$, then（在某些情况可以取等）
$$
\CVaR_{\alpha}(\delta)\ge\sup_{\beta}\left(\beta+\frac{1}{2\alpha}\mathbb{E}[(\tau(X)-b-\beta)_-]+\frac{1}{2\alpha}\mathbb{E}[(\tau(X)+b-\beta)_-]\right)
$$
证明核心：$(x-\beta)_-$ 是上凸函数，$\delta=\frac{\tau(X)+b-\delta}{2b}(\tau(X)-b)+\frac{\delta-\tau(X)+b}{2b}(\tau(X)+b)$ 所以
$$
\mathbb{E}[(\delta-\beta)_-|X]\ge\frac{1}{2}\mathbb{E}[(\tau(X)-b-\beta)_-|X]+\frac{1}{2}\mathbb{E}[(\tau(X)+b-\beta)_-|X]
$$
**Theorem3.3** Suppose $\tau(X)-\delta\le b$, then（在某些情况可以取等）
$$
\CVaR_{\alpha}(\delta)\ge \CVaR(\tau(X))-b
$$

也是对 $x_-$ 进行简单放缩。

#### 3.3 Lower Bounds under Limited Residual Heterogeneity Variance

上面的绝对偏差限制会被个案影响到。

**Theorem3.4** Suppose  $\Var[\delta|X]\le\bar{\sigma}^2(X)$ for some integrable $\bar{\sigma}^2:\mathcal{X}\to\mathbb{R}_+$. Then（某些情况可以逼近）
$$
\CVaR_{\alpha}(\delta)\ge\sup_{\beta}\left(\beta+\frac{1}{2\alpha}\E[\tau(X)-\beta-\sqrt{(\tau(X)-\beta)^2+\bar{\sigma}^2(X)}]\right)
$$
作者的办法很复杂，用到一些Sion minimax theorem和凸优化，没读。

如果只要得到下界，一步Cauchy-Schwarz就可以做到，但是没法推出这个界是精确的。

Corollary3.5 
$$
0\le \CVaR_\alpha(\tau(X))-\CVaR(\delta)\le\frac{1}{2\alpha}\E[\sqrt{\Var[\delta|X]}]
$$

所以作者后面都在研究 $\CVaR(\tau(X))$.

### 4 INFERENCE

Target: $\Psi=\CVaR_{\alpha}(\tau(X))=\beta^*+\frac{1}{\alpha}\E[(\tau(X)-\beta^*)_-]$ where $\beta^*=F_{\tau(X)}^{-1}(\alpha)=\inf\{\beta:\mathbb{P}[\tau(X)\le\beta]\ge\alpha\}$.

最容易想到的实践方法是先估计 $\tau$ 得到 $\hat{\tau}$ ，然后得到 $\hat{\Psi}_{plug-in}=\sup_\beta\{\beta+\frac{1}{n\alpha}\sum_{i=1}^n (\hat{\tau}(X_i)-\beta)_-\}$. 但是如果 $\tau$ 估计的不好/收敛慢，会导致 $\hat{\Psi}_{plug-in}$ 估计的不好。

作者给出的统计方法是分别估计 $e,\tau,\mu,\beta^*$ 再整合（注意是四个参数全部直接从数据中学习），这里的动机应该是教材上double-robust的方法。
$$
\phi(X,A,Y;\hat{e},\hat{\mu},\hat{\tau},\hat{\beta})=\widehat\beta+
\frac{\mathbf1\{\widehat\tau(X)\le\widehat\beta\}}{\alpha}
\left[
\widehat\mu(X,1)-\widehat\mu(X,0)
+
\frac{A-\widehat e(X)}
{\widehat e(X)(1-\widehat e(X))}
\bigl(Y-\widehat\mu(X,A)\bigr)
-\widehat\beta
\right].
$$

$$
\hat{\Psi}=\frac{1}{n}\sum_{i=1}^n\phi_i\\
\hat{se}=\sqrt{\frac{1}{n(n-1)}\sum_{i=1}^n(\phi_i-\hat{\Psi})^2}
$$

并在估计 estimators 的时候使用cross-fitting。

作者画饼：可以让这个方法的结果对四类参数全部正交（e.g. $\frac{d}{dt}\big|_{t=0}\E[\phi(X,T,Y;\eta+th)]=0$），从而误差的一阶项全部消失，有可能达到希望的 $o_p(n^{-1/2})$ 收敛（因为CLT给的置信区间也只有这个收敛速度，达到 $n^{-1/2}$ 一般认为最优）。

remark7：附录A给了估计lower bound的统计方法，二者组合可以给出覆盖 ITE-CVaR 的置信区间。

remark8：如果对多个 $\alpha$ 计算，$\CVaR_{\alpha}$ 总单调递增，但是我们统计出来的 $\alpha$ 未必递增，给他们重排成递增的，就可以减少整体的（平方损失）误差。应该是个标准方法。

#### 4.1 Local Robustness and Confidence Intervals

Assumption1：$\bar{e}\le e\le 1-\bar{e}, \bar{e}>0$, $|Y|\le B$, $F_{\tau(X)}$ is $C^1$ at $F_{\tau(X)}^{-1}(\alpha)$.

Lemma4.1：For each fold $k=1,\cdots,K$, the cross-fitting estimator $\hat{\beta}^{(k)}$ satisfies
$$
|\hat{\beta}^{(k)}-\beta^*|=O_p(n^{-1/2}\or ||\hat{\tau}^{(k)}-\tau||_r^{\frac{r}{r+1}}),\forall r\in[1,\infty]
$$
$r<\infty$ 是 nontrivial 的。

ai觉得附录的证明伪了。一个小问题是assumption不够强，delta method未必能用（可能出现除以0）。

大问题是，作者写了一步
$$
\left(
\frac1{|I_{-k}|}
\sum_{i\in I_{-k}}
|\widehat\tau^{(k)}(X_i)-\tau(X_i)|^r
\right)^{1/r}
=
O_p\!\left(\|\widehat\tau^{(k)}-\tau\|_r\right).
$$
$\hat{\tau}^{(k)}$ 是 $I_{-k}$ 的数据训练的，右边是直接对真实分布积分。可以构造非常病态的训练数据让这个不成立（例如刚好在训练数据上全部偏差了1，但是样本至多可数，所以零测，不影响整体，那么左边会远大于右边）。

但是ai告诉我论文主要用的是 $r=\infty$，所以我继续看。注意 $r=\infty$ 的时候仍然要用delta method，所以需要 $F'_{\tau(X)}(\beta^*)\not=0$，不过这个感觉是小问题。

**Theorem4.2** 假设 Assumption1，for $k=1,\cdots,K$, $||\hat{e}^{(k)}-e||_2=o_p(1),||\hat{\mu}^{(k)}-\mu||_2=o_p(1)$, $||\hat{e}^{(k)}-e||_2||\hat{\mu}^{(k)}-\mu||_2=o_p(n^{-1/2})$, $||\hat{\tau}^{(k)}-\tau||_{\infty}=o_p(n^{-1/4})$, $\mathbb{P}[||\hat{\mu}^{(k)}||_{\infty}\le B]\to 1$, $\mathbb{P}[\bar{e}\le \hat{e}^{(k)}\le 1-\bar{e}]\to1$. Then
$$
\hat\Psi=\frac{1}{n}\sum_{i=1}^n\phi(X,A,Y;e,\mu,\tau,\beta^*)+o_p(n^{-1/2})=\Psi+o_p(n^{-1/2})
$$

$$
\mathbb{P}(\Psi\in[\hat\Psi\pm\Phi^{-1}(\frac{1+\gamma}{2})\hat{se}])\to\gamma, \forall \gamma
$$

#### 4.2 Double Robustness and Double Validity

接下来建立鲁棒性，如果 $\hat{\mu},\hat{e}$ 有一个收敛到错误的极限也可以得出正确的结果。

**Theorem4.3** Fix any $\widetilde e,\widetilde\mu$ with $\bar e \le \widetilde e \le 1-\bar e$ and $\|\widetilde\mu\|_\infty \le B$. Let $r_n \to 0$ be a deterministic sequence. Suppose Assumption 1 holds and for $k=1,\ldots,K$, $||\hat{e}^{(k)}-\tilde{e}||_2=o_p(1)$, $||\hat{\mu}^{(k)}-\tilde{\mu}||_2=o_p(1)$, $||\hat{\tau}^{(k)}-\tau||_\infty=O_p(r_n^{1/2})$, $\mathbb{P}[||\hat{\mu}^{(k)}||_{\infty}\le B]\to 1$, $\mathbb{P}[\bar{e}\le \hat{e}^{(k)}\le 1-\bar{e}]\to1$, and

either $||\hat{e}^{(k)}-e||_2=O_p(r_n)$ or $||\hat{\mu}^{(k)}-\mu||_2=O_p(r_n)$.

Then $\hat{\Psi}=\Psi+O_p(r_n\or n^{-1/2})$.

如果 CATE-CVaR 估计的也不好，我们也可以有控制（如果 $e,\mu$ 之一估计的好）

 **Theorem4.5** Fix any $\widetilde e,\widetilde\mu,\tilde{\tau}$ with $\bar e \le \widetilde e \le 1-\bar e$ and $\|\widetilde\mu\|_\infty \le B,||\tilde{\tau}||_{\infty}\le2B$. Let $r_n \to 0$ be a deterministic sequence. Suppose Assumption 1(replace $\tau$ by $\tilde{\tau}$) holds and for $k=1,\ldots,K$, $||\hat{e}^{(k)}-\tilde{e}||_2=o_p(1)$, $||\hat{\mu}^{(k)}-\tilde{\mu}||_2=o_p(1)$, $||\hat{\tau}^{(k)}-\tilde\tau||_\infty=O_p(r_n)$, $\mathbb{P}[||\hat{\mu}^{(k)}||_{\infty}\le B]\to 1$, $\mathbb{P}[\bar{e}\le \hat{e}^{(k)}\le 1-\bar{e}]\to1$, and

either $||\hat{e}^{(k)}-e||_2=O_p(r_n)$ or $||\hat{\mu}^{(k)}-\mu||_2=O_p(r_n)$.

Then $\hat{\Psi}\ge\Psi-O_p(r_n\or n^{-1/2})$.

#### 4.3 CATE-Estimation and Rates

主要讨论了什么学习方法可以满足上面的收敛条件，引了很多论文结果。看了晕晕。

### 5 CASE STUDY

#### 5.1 Background and Setup

 如果把私人机构提供的就业辅导改成公共机构提供，会怎样影响失业者的再就业？ 

$A=0:$  分配到私人机构强化辅导项目

$A=1$：分配到公共机构强化辅导项目

$Y=1/0$：六个月内是/否重新就业

$\widehat{ATE}$ 为 1.22%（90% CI $[-0.35,2.8]$）。数据没有明确排除负效应或者零效应。如果公共就业辅导可以节约成本，那么政策制定者也可以接受，但是，这个结论是否掩盖了某些人群的不利影响。

作者使用4.3.1（对已知 $e$ 的回归方法）的方法处理了数据，作者自己分了很多协变量。

#### 5.2 Upper bounds

作者枚举 $\alpha$ 计算 90% CI得到了Fig1，识别出56%的群体处理效果为负。而且计算的是CATE-CVaR，也就是真实的CVaR还会更小。

Fig2是最不受益的这部分人与平均的差距。

Fig3是只用少量协变量来run，就没有检测到Fig1的结论，所以这种方法在协变量分的不好的时候也没啥用。（这不是显然吗）

#### 5.3 Lower bounds

做了“敏感性分析”，用thm3.4和cor3.5（3.3比较难控制），看看允许多少的方差偏差下，对最差的 $\alpha$ 部分的人的判断会有多少改变。

使用cor3.5，需要 $\Var[\delta|X]$，但是算这个需要 $Y(1),Y(0)$ 的correlation，不可识别，所以作者枚举这个值画了Fig5.

### My questions

+ 作者的研究重心非常快的从 $\CVaR_{\alpha}(\delta)$ 转到 $\CVaR_{\alpha}(\tau(X))$ ，二者的联系只有一个简单的不等式（thm3.1），虽然只给定CATE的时候界是精确的，但是给定整个 $(X,A,Y)$ 分布之后未必是精确的，光凭这个显然的“取等”并不让我感觉研究CATE-CVaR可以很好地代替ITE-CVaR。或者说，整篇文章更像是研究了CATE-CVaR然后推出了ITE-CVaR的性质，而不是在研究CATE-CVaR。
+ lemma4.1怎么还伪了，不过ai检查了，没用到 $r<\infty$ 的情况，就当是小问题。
+ 假设unconfoundedness，合理吗。














