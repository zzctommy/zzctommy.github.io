---
date: 2026-11-07
title: 3 Doubly-Valid Doubly-Sharp Sensitivity Analysis for Causal Inference with Unmeasured Confounding
---

[论文地址](https://doi.org/10.1080/01621459.2024.2335588) published on JASA 2025

作者信息：

Kevin Guo：2023PhD in Princeton. Now OpenAI

Jacob Dorn: Assistant professor  at Cornell University

Nathan Kallus:  Associate Professor at Cornell Tech, Cornell University, and Director of Machine Learning & Inference Research at Netflix.

---

### 1. Introduction

They will establish a new estimator with double sharpness and double validity(under weaker condition).

他们认为double sharpness是半参回归的标准结论，double validity是他们的主要贡献。

### 2. Sensitivity Assumption and Problem Statement

#### 2.1. Setting

$n$ units $(X_i,Y_i(0),Y_i(1),Z_i)$ are sampled independently from a distribution $P_{full}$.

$X_i\in\mathbb{R}^d$ is a vector of measured confounders, $(Y_i(0),Y_i(1))\in\mathbb{R}^2$ are potential outcomes, $Z_i\in\{0,1\}$ treatment indicator. The data analyst only observe $(X_i,Y_i,Z_i)$ (distribution $P$) where $Y_i=Y_i(Z_i)$.

Goal: ATE $\psi_{ATE}(P_{full})=\mathbb{E}_{P_{full}}[Y(1)-Y(0)]$; counterfactual means $\psi_z(P_{full})=\mathbb{E}_{P_{full}}[Y(z)],z\in\{0,1\}$.

#### 2.2. The Marginal Sensitivity Model

$MSM(P,\Lambda):$ Those $P_{full}$ incudes $P$ and s.t.
$$
\Lambda^{-1}\le \frac{P_{full}(Z=1|X,Y(0),Y(1))/P_{full}(Z=0|X,Y(0),Y(1))}{P_{full}(Z=1|X)/P_{full}(Z=0|X)}\le \Lambda
$$
Compared with Rosenbaum's sensitivity model $R(\Gamma):$ Assume $Z\perp (Y(0),Y(1))|X,U$.
$$
\Gamma^{-1}\le \frac{P(Z=1|X,U=u)/P(Z=0|X,U=u)}{P(Z=1|X,U=u')/P(Z=0|X,U=u')}\le\Gamma
$$
Relation:
$$
MSM(\sqrt{\Gamma})\subseteq\R(\Gamma)\subseteq MSM(\Gamma)
$$
proof: The first inclusion. take $U=(Y(0),Y(1))$.

The second inclusion. 


$$
P_{full}(Z=1|X,Y(0),Y(1))
=\mathbb{E}[Pr(Z=1|X,U,Y(0),Y(1))|X,Y(0),Y(1)]\\
=\mathbb{E}[Pr(Z=1|X,U)|X,Y(0),Y(1)]
$$

$$
e(X)=Pr(Z=1|X)=\mathbb{E}[Pr(Z=1|X,U)|X]
$$

Since for each $X,U$ we have the Rosenbaum's control, take average we get the MSM.    $\square$

A quoted result says: the assumption $P_{full}\in MSM(P,\Lambda)$ restricts $\psi_0,\psi_1,\psi_{ATE}$ to closed intervals known as partially identified set
$$
[\psi_{\diamond}^-,\psi_{\diamond}^+]=[\inf_{P_{full}\in MSM(P,\Lambda)}\psi_{\diamond}(P_{full}),\sup_{P_{full}\in MSM(P,\Lambda)}\psi_{\diamond}(P_{full})],\diamond\in\{0,1,ATE\}
$$

#### 2.3. Notation

固定 $\Lambda\ge 1,\tau=\frac{\Lambda}{1+\Lambda}$.

propensity score $e(x)=P(Z=1|X=x)$.

outcome regression $\mu(x,z)=\mathbb{E}[Y|X=x,Z=z]$.

cdf of $\mu$: $F(y|x,z)=P(Y\le y|X=x,Z=z)$.

conditional outcome quantiles: $Q_+(x,z)=\inf \{q:F(q|x,z)\ge \tau \}$; $Q_-(x,z)=\inf\{q:F(q|x,z)\ge 1-\tau \}$.

Assume also $Y$ is $L^1$, $P(0<e(X)<1)=1$.

### 3. Partial Identification and Singly-Valid Estimation

#### 3.1. Formulation Using Adversarial Propensities

**Proposition1** the worst case of confounding structure. There exist a function $e_+:\mathbb{R}^d\times \mathbb{R}\to(0,1)$ s.t.
$$
\frac{e_+(x,y)}{1-e_+(x,y)}/\frac{e(x)}{1-e(x)}=\begin{cases}\Lambda^{-1}, y>Q_+(x,1)\\
\gamma(x,y), y=Q_+(x,1)\\ \Lambda, y<Q_+(x,1) \end{cases}
$$
where $\gamma$ can be chosen s.t.  $\mathbb{E}[\frac{Z}{e_+(x,y)}|X]=1$. Then $\psi_1^+=\mathbb{E}[YZ/e_+(X,Y)]$.

证明：作者引用了别人的结果。

也就是依照治疗效果重新分配权重，治疗效果好的人被治疗的概率小，权重 $1/e_+$ 就大，会把大结果解释为，很难观察到，但是观察到就可以代表很多人，从而推高 $\psi_1$.

Add $0=\mathbb{E}[Q_+(X,1)Z/e(X)-Q_+(X,1)Z/e_+(X,Y)]$ to get （注意 $Y=Q_+(X,1)$ 的时候第二项为0，所以等式成立。这个办法消去了未知的 $\gamma$，只需要 $e,Q_+$ 就可以计算）
$$
\psi_1^+=\mathbb{E}[\frac{Q_+(X,1)Z}{e(X)}+\frac{(Y-Q_+(X,1))Z}{e_+(X,Y)}]\\
=\mathbb{E}[\frac{Q_+(X,1)Z}{e(X)}+(Y-Q_+(X,1))Z\cdot(1+\frac{1-e(X)}{e(X)}\Lambda^{sgn(Y-Q_+(X,1))})]
$$
**Proposition2** Plug any integrable $\hat{Q}_+$ we get upper bound for $\psi_1^+$:
$$
\psi_1^+\le\mathbb{E}[\frac{\hat{Q}_+(X,1)Z}{e(X)}+(Y-\hat{Q}_+(X,1))Z\cdot(1+\frac{1-e(X)}{e(X)}\Lambda^{sgn(Y-\hat{Q}_+(X,1))})]
$$
证明：依旧引用同一个结果，探讨 $e_+$ 的性质进行放缩。

#### 3.2. Formulation Using Adversarial Regressions

$$
\mathbb{E}_{P_{full}}[Y(1)]=\mathbb{E}[ZY+(1-Z)\mathbb{E}_{P_{full}}[Y(1)|X,Z=0]]
$$

按照是否被治疗分组，被治疗的可观测，控制组需要估计。
$$
\mathbb{E}_{P_{full}}[Y(1)|X,Z=0]=\int y dF_{full}(y|x,0)
$$
由 MSM以及一些简单恒等变形 $\Lambda^{-1}\le\frac{dF_{full}(y|x,0)}{dF(y|x,1)}\le\Lambda$（本质是把 $Y(0)$ 平均掉了，仍然有MSM的界）

现在问题变成，给已知分布 $dF(y|x,1)$，乘以一个取值在 $[\Lambda^{-1},\Lambda]$ 之内的倍率函数，使得总概率仍为 $1$。如何分配使得均值最大。

设 $F_{full}(y|x,0)=\Lambda^{-1}F(y|x,1)+(1-\Lambda^{-1})G(y)$, $dG(y)/dF(y|x,1)\le\frac{1}{1-\tau}=1+\Lambda$ 

因此控制组最大可能的 $\mathbb{E}_{P_{full}}[Y(1)|X,Z=0]$ 是
$$
\Lambda^{-1}\int y dF(y|x,1) +(1-\Lambda^{-1})\sup_{\frac{dG(y)}{dF(y|x,1)}\le \frac{1}{1-\tau}}\int ydG(y)
$$
note（答案当然是找 dF最大的 $1-\tau$ 部分集体放大 $\frac{1}{1-\tau}$ 倍）
$$
CVaR_{\tau}(x,z)=\sup_{\frac{dG(y)}{dF(y|x,1)}\le \frac{1}{1-\tau}}\int ydG(y)\\
=\inf_{q\in\mathbb{R}}\int\{q+\frac{1}{1-\tau}(y-q)_+ \}dF(y|x,z)
$$
the  $\inf$ is attained at $q=Q_+(x,z)$.

**Proposition3** Define $\rho^+(x,1)=\Lambda^{-1}\mu(x,1)+(1-\Lambda^{-1})CVaR_{\tau}(x,1)$, then $\psi_1^+=\mathbb{E}[ZY+(1-Z)\rho^+(x,1)]$. 之前的讨论保证了这个是上界，追溯等号可知有取等。

Specifically
$$
\rho^+(x,z)=\int\left(\Lambda^{-1}y+(1-\Lambda^{-1})(Q_{+}(x,z)+\frac{1}{1-\tau}(y-Q_+(x,z))_+)\right)dF(y|x,z)
$$
于是实操上，先估计quantile $\hat{Q}_+$ 再跑回归得到  $\rho^+$ 的估计。
$$
\Lambda^{-1}Y_i+(1-\Lambda^{-1})(\hat{Q}_{+}(X_i,Z_i)+\frac{1}{1-\tau}(Y_i-\hat{Q}_+(X_i,Z_i))_+)\sim(X_i,Z_i)
$$
**Proposition4** For any integrable $\hat{Q}_+(x,z)$, define
$$
\varrho_+(x,z;\hat{Q}_+)=\int\left(\Lambda^{-1}y+(-\Lambda^{-1})(\hat{Q}_{+}(x,z)+\frac{1}{1-\tau}(y-\hat{Q}_+(x,z))_+)\right)dF(y|x,z)
$$
then $\psi_1^+\le\mathbb{E}[ZY+(1-Z)\varrho_+(x,1;\hat{Q}_+)]$.

我感觉这个和prop3说的几乎是同一件事，$Q_+$ 都是minimizer了这两个命题不是显然等价吗，。

### 4. The Doubly-Valid/Doubly-Sharp Estimators

#### 4.1. First-Step: Nuisance Estimation

DVDS estimator需要估计三类辅助函数：propensity score（prop2），transformed outcome-regressions（prop4），conditional quantiles（两个prop都有 ）。

设 $\rho_{\pm}(x,z)=\varrho(x,z;Q_{\pm})$，$\eta=(e,Q_+,Q_-,\rho_+,\rho_-)$ full vector of true DVDS nuisance functions.

Remark：quantile regression：estimate $Q_{\alpha}(x,z)=\inf\{q:F(q|x,z)\ge \alpha\}$. 

回归方法是，$Q_{\alpha}(x,z)\in\arg\min_q\mathbb{E}[l_{\alpha}(Y-q)|X=x,Z=z]$ where loss function $l_{\alpha}(u)=u(\alpha-1_{u<0})$.

+ 对continuous outcome：估计 $\hat{e},\hat{Q}_{\pm}$，然后用 $\hat{\rho}_{\pm}$(prop3之后建议的回归方法)估计 $\rho_{\pm}$.
+ binary outcome: $\rho_\pm$ 的积分变成两项求和，有好看的形式。

注意所有回归用cross-fitting防止overfitting。

#### 4.2. Second-Step: Influence Function Aggregation

divide data into $K$ folds: $\mathcal{F}_1,\cdots,\mathcal{F}_K$.
$$
\hat{\psi}_{\diamond}^{\pm}=\frac{1}{n}\sum_{k=1}^K\sum_{i\in\mathcal{F}_k}\hat{\phi}_{\diamond}^{\pm}(X_i,Y_i,Z_i;\hat{\eta}^{(k)})\quad \diamond\in\{0,1,ATE\}
$$
where 

$$
\hat{\phi}_{1}^{\pm}(x,y,z;\eta)=zy+(1-z)\hat{\rho}_{\pm}(x,1)\\
+\frac{(1-\hat{e}(x))z}{\hat{e}(x)}(\hat{Q}_{\pm}(x,1)+\Lambda^{\pm sgn(y-\hat{Q}_{\pm}(x,1))}(y-\hat{Q}_{\pm}(x,1))-\hat{\rho}_{\pm}(x,1))
$$



$$
\hat{\phi}_{0}^{\pm}(x,y,z;\eta)=(1-z)y+z\hat{\rho}_{\pm}(x,0)\\
+\frac{\hat{e}(x)(1-z)}{1-\hat{e}(x)}(\hat{Q}_{\pm}(x,0)+\Lambda^{\pm sgn(y-\hat{Q}_{\pm}(x,0))}(y-\hat{Q}_{\pm}(x,0))-\hat{\rho}_{\pm}(x,0))
$$

$$
\hat{\rho}_{ATE}^{\pm}=\hat{\phi}_1^{\pm}-\hat{\phi}_0^{\mp}
$$

#### 4.3. Intuition

先看 $Z=1$ 的情况，设
$$
\tilde{Y}=\Lambda^{-1}Y+(1-\Lambda^{-1})(\hat{Q}_{+}(X,1)+\frac{1}{1-\tau}(Y-Q_+(X,1))_+)
$$
recall常规AIPW:  $\mathbb{E}_X[\mathbb{E}[Y|X,Z=1]]\leftrightarrow \hat{\mu}(X,1)+\frac{Z}{\hat{e}(X)}(Y-\hat{\mu}(X,1))$, $\mu(X,1)=\mathbb{E}[Y|X,Z=1]$，关于 $\hat{e},\hat{\mu}$ 双鲁棒。

Note $\mathbb{E}[\tilde{Y}|X,Z=1]=\varrho_+(X,1;\hat{Q}_+)$.

So $\mathbb{E}_X[\varrho_+(X,1;\hat{Q}_+)]\leftrightarrow \hat{\rho}_+(X,1)+\frac{Z}{\hat{e}(X)}(\tilde{Y}-\hat{\rho}_+(X,1))$.

现在尝试估计 $\mathbb{E}[ZY+(1-Z)\varrho_{+}(X,1;\hat{Q}_+)]$.

注意到 $\mathbb{E}[Z\varrho_+(X,1;\hat{Q}_+)]=\mathbb{R}[Z\tilde{Y}]$.（两边condition on $Z,X$）

所以等价于估计 $\mathbb{E}[Z(Y-\tilde{Y})+\varrho_+(X,1;\hat{Q}_+)]$，第一项直接取样本平均，第二项用AIPW。

整理（省略下标）：
$$
ZY-Z\tilde{Y}+\hat{\rho}_++\frac{Z}{\hat{e}}(\tilde Y-\hat{\rho}_+)\\
=ZY+(1-Z)\hat{\rho}_++\frac{(1-\hat{e})Z}{\hat{e}}(\tilde{Y}-\hat{\rho}_+)
$$
然后把 $\tilde{Y}$ 展开就得到了作者给出的estimator。

此时直观上有这个estimator关于 $\hat{e}, \varrho_+$ double robust。

设
$$
 \Psi_1^+(\widehat Q_+) = E\!\left[ ZY+(1-Z)\varrho_+(X,1;\widehat Q_+) \right]
$$

那么
$$
\hat{\psi}_1-\psi_1=(\hat{\psi}_1-\Psi_1^+(\hat{Q}_+))+(\Psi_1^+(\hat{Q}_+)-\Psi_1^+(Q_+))
$$
第一项由之前的AIPW控制；第二项，由于 $Q_+$ 是minimizer，所以期望 $\Psi_1^+$ 附近泰勒展开一阶项消失，误差小。

### 5. Theory for Estimation of Bounds

For brevity, we only state results for the ATE bounds $(\hat{\psi}_{ATE}^-,\hat{\psi}_{ATE}^+)$.

#### 5.1. Sharpness and Validity

**Condition1**: outcome $Y$ has $q>2$ finite moments. For some $\epsilon>0$, $\mathbb{P}[\epsilon<e(X)<1-\epsilon]=1$.

**Condition2**: $||\hat{Q}_{\pm}||_q=O_p(1),||\hat{\rho}_{\pm}||_q=O_p(1)$, $\mathbb{P}[\epsilon<\inf_x\hat{e}(x),\sup_x\hat{e}(x)<1-\epsilon]\to1$.

**Theorem1** (double sharpness and double validity) Assume condition1,2.

+ Suppose $||\hat{Q}_+-Q_+||_2\lor||\hat{Q}_--Q_-||_2=o_P(1)$, so that quantiles are estimated consistently. Then the DVDS estimators $(\hat{\psi}_{ATE}^-,\hat{\psi}_{ATE}^+)$ converge to the sharp ATE bounds in probability if either $||\hat{e}-e||_2=o_P(1)$ or $||\hat{\rho}_+-\rho_+||_2\lor||\hat{\rho}_--\rho_-||_2=o_P(1)$.
+ Supppose $||\hat{Q}_+-Q_+||_2\lor||\hat{Q}_--Q_-||_2\not=o_P(1)$, so  that quantile models are mis-specified. Then $\hat{\psi}_{ATE}^+\ge\psi_{ATE}^+-o_P(1)$ and $\hat{\psi}_{ATE}^-\le\psi_{ATE}^-+o_P(1)$ if either $||\hat{e}-e||_2=o_P(1)$ or $||\hat{\rho}_+-\rho_+||_2\lor||\hat{\rho}_--\rho_-||_2=o_P(1)$.

#### 5.2. Asymptotic Normality and Efficiency

**Condition3**: One of the following holds.

+ continuous case: For each $(x,z)$, $F(y|x,z)$ is continuous with a uniformly bounded density $f(y|x,z)$ that is positive near $\tau$ and $1-\tau$ quantiles.
+ binary case: $Y$ is binary and the regression function $\mu(X,Z)$ has bounded density.

**Theorem2** (Asymptotic normality and semiparametric efficiency) Assume Condition1-3 and at

+ continuous case: $||\hat{Q}_+-Q_+||_2=o_P(n^{-1/4})$ and $||\hat{e}-e||_2\times||\hat{\rho}_{\pm}-\rho_{\pm}||_2=o_P(n^{-1/2})$.
+ binary case: $||\hat{\mu}-\mu||_{\infty}=o_P(n^{-1/4})$ and $||\hat{e}-e||_2\times||\hat{\mu}-\mu||=o_P(n^{-1/2})$.

Then the DVDS estimators are $\sqrt{n}$ - consistent, asymptotically normal and semiparametrically efficient.
$$
\sqrt{n}\begin{pmatrix}\hat{\psi}_{ATE}^--\psi_{ATE}^-\\ \hat{\psi}_{ATE}^+-\psi_{ATE}^+ \end{pmatrix}=\frac{1}{\sqrt{n}}\sum_{i=1}^n \begin{pmatrix}\phi_{ATE}^-(X_i,Y_i,Z_i;\eta)-\psi_{ATE}^-\\ \phi_{ATE}^+(X_i,Y_i,Z_i;\eta)-\psi_{ATE}^+ \end{pmatrix}+o_P(1)\overset{d}{\to} N(0,\Sigma)
$$

### 6. Theory for Inference on Bounds

standard error
$$
\hat\sigma_\pm^2
=
\frac1{n(n-1)}
\sum_{k=1}^K\sum_{i\in F_k}
\left[
\phi^\pm_{\mathrm{ATE}}(X_i,Y_i,Z_i;\hat\eta^{(k)})
-\hat\psi^\pm_{\mathrm{ATE}}
\right]^2
$$
**Theorem3** Assume the conditions of theorem2. Then $\mathbb{P}[\psi_{ATE}^+\le\hat{\psi}_{ATE}^++z_{1-\alpha}\hat{\sigma}_+]\to1-\alpha$. Moreover, for testing the null hypothesis $\hat{\psi}_+\le 0$ at level $\alpha$, the test is asymptotically most powerful.





**Theorem4** Assume condition 1-3 and the following:

+ continuous case: $||\hat{e}-e||_2\times ||\rho_{\pm}-\varrho(\cdot,\cdot;\hat{Q}_{\pm})||_2=o_P(n^{-1/2})$, and $\hat{Q}_{\pm}$ tends to a possibly mis-specified limit $\bar{Q}_{\pm}$.
+ binary case: $||\hat{e}-e||_2\times||\hat{\mu}-\mu||_2=o_P(n^{-1/2})$.

Then $\mathbb{P}[\psi_{ATE}^+\le\hat{\psi}_{ATE}^++z_{1-\alpha}\hat{\sigma}_+]\ge1-\alpha-o(1)$ for $\forall \alpha\in(0,1)$.

*"We consider this to be the most surprising finding of this article."*

尽管estimator未必渐进正态，但是正态分位数仍然适用，这是因为与渐进正态的oracle相比，这个估计只会更保守。

### My questions

+ 如何陈述双鲁棒性，该加什么条件得到什么结论。
+ 我没有看任何有关双鲁棒性的证明。感觉这个才是文章的核心，但是我未必用得到。


