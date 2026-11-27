---
date: 2026-11-24
title: Adaptive Robust Confidence Intervals
---

[论文地址](https://arxiv.org/abs/2410.22647)

2026.6发表在  The Annals of Statistics；2024.10 arxiv

作者信息： 

**Yuetian Luo** 现为罗格斯大学统计系助理教授。2017 年本科毕业于中国人民大学统计学专业，2022 年在威斯康星大学麦迪逊分校获得统计学博士学位，导师为 Anru Zhang；随后在芝加哥大学与 Rina Foygel Barber 从事博士后研究。他的研究方向包括高维统计推断、无分布假设的推断和稳健统计。 

 **Chao Gao** 是芝加哥大学统计系教授。他毕业于耶鲁大学，师从 Harrison H. Zhou，研究方向包括稳健统计、高维与非参数统计、网络分析和贝叶斯理论。 [个人主页](https://www.stat.uchicago.edu/~chaogao/).

---

这篇文章相比于我上周搞出来的方法，巨大的改进是，如果我们只知道 $\epsilon\le 0.05$，但是实际的 $\epsilon=0$，我的方法会导致置信区间长度趋近于一个常数（由污染上界确定），但是作者的方法会导致置信区间长度以 $1/\sqrt{\log n}$ 的速度趋近于 $0$（在单参数正态分布的case）。这才是真正意义上的 adaptive。

---

## 1 Introduction

Classical result：$X_1,\cdots,X_n\overset{iid}{\sim}N(\theta,1)$. Then $[\bar{X}-\frac{1.96}{\sqrt{n}},\bar{X}+\frac{1.96}{\sqrt{n}}]$ is $95\%$ CI and have the shortest length.

Will consider Huber's contamination model
$$
X_1,\cdots,X_n\overset{iid}{\sim}P_{\epsilon,\theta,Q}=(1-\epsilon)N(\theta,1)+\epsilon Q
$$
A result says (use sample median to estimate)
$$
\inf_{\theta,Q}P_{\epsilon,\theta,Q}[|\hat{\theta}_{median}-\theta|\le\sqrt{2.1\pi}(\frac{1.36}{\sqrt{n}}+\epsilon)]\ge 0.95
$$
And $\frac{1}{\sqrt{n}}+\epsilon$ is known as the optimal rate.

一个被诟病的点是，区间长度依赖于 $\epsilon$，现实中没用。并且，由于 $Q$ 是任意的，$\epsilon$ 是不可识别的。（例如上面的例子，$Q=N(\theta,1)$ 不改变观测）

### 1.3 Notation

$[n]:=\{1,2,\cdots,n\}$.

For two sequences $\{a_n\},\{b_n\}$, write $a_n\asymp b_n$ if exists constants $c,C>0$ s.t. $ca_n\le b_n\le Ca_n$.

$a_n\lesssim b_n$ if exists $C>0$ s.t. $a_n\le Cb_n$.

Given cdf $F(t)$, $q$-quantile is defined as $F^{-1}(q)=\inf\{t\in \mathbb{R};F(t)\ge q\}$.

Given data $X_1,\cdots,X_n$, empirical cdf $F_n(t)=\frac{1}{n}\sum_{i\in[n]}I(X_i\le t)$, $F_n^{-1}(q)=X_{\lceil nq\rceil}$.

## 2 Problem Setting

Given iid observations
$$
X_1,\cdots,X_n\overset{iid}{\sim}P_{\epsilon,\theta,Q}=(1-\epsilon)P_\theta+\epsilon Q
$$
with some parametric model $P_{\theta}$ and a pre-specified error probability $\alpha\in(0,1)$, want to construct $\widehat{CI}=\widehat{CI}(\{X_i\}_{1\le i\le n})$ that covers $\theta$ with probability at least $1-\alpha$, while optimizing its length quantified by Lebesgue measure.

Since $\epsilon$ is unknown, we require the CI be uniform over $\epsilon\in\mathcal{E}\subseteq [0,1]$.

Call $\widehat{CI}$ be $(1-\alpha)$-level adaptive robust confidence interval(ARCI) over $\mathcal{E}$, if
$$
\inf_{\epsilon\in\mathcal{E}}\inf_{\theta,Q}(\theta\in\widehat{CI})\ge 1-\alpha
$$
The set of $(1-\alpha)$-level ARCIs is defined by
$$
\mathcal{I}_\alpha(\mathcal{E})=\{\widehat{CI}=[l(\{X_i\}),u(\{X_i\})]:\inf_{\epsilon\in\mathcal{E}}\inf_{\theta,Q}(\theta\in\widehat{CI})\ge 1-\alpha \}
$$


The optimal length at $\epsilon\in\mathcal{E}$ over all $(1-\alpha)$-level ARCIs is（高概率长度不超过 $r$）
$$
r_{\alpha}(\epsilon, \mathcal{E})=\inf\left\{r\ge 0: \inf_{\widehat{CI}\in\mathcal{I}_{\alpha}(\mathcal{E})}\sup_{\theta,Q}P_{\epsilon,\theta,Q}(|\widehat{CI}|\ge r)\le \alpha \right\}
$$
We will mainly focus on $\mathcal{E}=[0,\epsilon_{max}]$. $P_{\epsilon,\theta,Q}\subseteq P_{\epsilon_{max},\theta,Q}$ implies $r_{\alpha}(\epsilon,[0,\epsilon_{max}])=r_{\alpha}(\epsilon,\{\epsilon_{max}\})$

Moreover
$$
r_{\alpha}(\epsilon,\{\epsilon\})\subseteq r_{\alpha}(\epsilon,\{\epsilon_{max}\})\subseteq r_{\alpha}(\epsilon_{max},\{\epsilon_{max}\})
$$

## 3 Adaptation Cost in Gaussian Location Model

In this section consider $P_{\theta}=N(\theta,1)$. Goal is to construct optimal ARCI. Will choose $\epsilon_{max}=0.05$ for convenience.

### 3.1 A Lower Bound of Optimal Length

**Theorem1** Suppose $n \ge 100$. For any $\alpha\in(0,1/4)$ exists constant $c=c(\alpha)$ s.t.
$$
r_{\alpha}(\epsilon, [0,0.05])\ge c\left(\frac{1}{\sqrt{\log n}}+\frac{1}{\log (1/\epsilon)} \right)
$$
for any $\epsilon\in[0,0.05]$.(recall if $\epsilon$ is known, the optimal rate is $\frac{1}{\sqrt{n}}+\epsilon$)

proof: 

> Lemma4：For any $n\ge 100, \alpha\in(0,1),\epsilon\in[0,0.05]$, there exists a consatnt $c>0$ only depends on $\alpha$ such that as long as
> $$
> r\le c\left(\frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}} \right)
> $$
> we have
> $$
> \inf_{Q_0,Q_1} TV(P_{0.05,\theta,Q_0}^{\otimes n},P_{\epsilon,\theta\pm r,Q_1}^{\otimes n})\le \alpha
> $$
> lemma4的直观：区间太短会导致无法区分 $\theta,\theta\pm r$. 证明是构造性的.

Now by lemma4 take such constant $c$.

Suppose “$r_\alpha(\epsilon, 0.05)>c(\frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}})$ for all $\epsilon\in[0,0.05]$ ” does not hold. That is, exists $\epsilon\in[0,0.05]$ and $\widehat{CI}$ s.t. （1. coverage for all $\epsilon$; 2. short for $\epsilon$）
$$
\inf_{\epsilon\in[0,0.05],\theta,Q}P_{\epsilon,\theta,Q}(\theta\in\widehat{CI})\ge 1-\alpha \quad and\quad \sup_{\theta,Q}P_{\epsilon,\theta,Q}(|\widehat{CI}|\ge r)\le \alpha
$$
with $r=c(\frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}})$.

> Proposition2：Suppose $\widehat{CI}$ is an interval satisfies(coverage)
> $$
> \inf_{\epsilon\in[0,\epsilon_{max}],\theta,Q}P_{\epsilon,\theta,Q}(\theta\in\widehat{CI})\ge 1-\alpha
> $$
>
> Then the testing functions $\phi_{\theta,\theta\pm r_{\epsilon},\epsilon}:=I(\theta\not\in\widehat{CI})$ satisfies（对所有 $\epsilon$ 用了同样的test）
> $$
> \sup_{Q}P_{\epsilon_{max},\theta,Q}(\sup_{\epsilon\in[0,\epsilon_{max}]}\phi_{\theta,\theta\pm r_\epsilon,\epsilon}=1)\le \alpha
> $$
> for all $\theta\in\mathbb{R}$. Suppose in addition we have(length)
> $$
> \inf_{\theta,Q}P_{\epsilon,\theta,Q}(|\widehat{CI}|\le r_{\epsilon})\ge 1-\alpha
> $$
> for some $\epsilon\in [0,\epsilon_{max}]$. Then for the same $\epsilon$ 
> $$
> \sup_{Q}P_{\epsilon,\theta\pm r_{\epsilon},Q}(\phi_{\theta,\theta\pm r_{\epsilon},\epsilon}=0)\le 2\alpha
> $$
> Where the label denotes the test 
> $$
> \mathcal{H}(\theta,\theta\pm r,\epsilon): H_0:P\in\{P_{\epsilon_{max},\theta,Q}:Q \}\quad v.s. \quad H_1:P\in\{P_{\epsilon,\theta\pm r,Q}:Q\}
> $$

Then given $\theta\in\mathbb{R}$ we can define a test s.t.
$$
\sup_{Q}P_{\epsilon_{max},\theta,Q}(\phi_{\theta,\theta-r,\epsilon}=1)+\sup_{Q} P_{\epsilon,\theta-r,Q}(\phi_{\theta,\theta-r,\epsilon}=0)\\
\le \sup_{Q}P_{\epsilon_{max},\theta,Q}(\sup_{\epsilon\in[0,\epsilon_{max}]}\phi_{\theta,\theta-r,\epsilon}=1)+\sup_{Q} P_{\epsilon,\theta-r,Q}(\phi_{\theta,\theta-r,\epsilon}=0)\le 3\alpha
$$
By the construction of $c$, for the test $\mathcal{H}(\theta,\theta-r,\epsilon)$
$$
\inf_{\phi\in\{0,1\}}\left\{\sup_{Q}P_{\epsilon_{max},\theta,Q}(\phi=1)+\sup_QP_{\epsilon,\theta-r,Q}(\phi=0) \right\}\\
\ge 1-\inf_{Q_0,Q_1} TV(P_{\epsilon_{max},\theta,Q_0}^{\otimes n},P_{\epsilon,\theta-r,Q_1}^{\otimes n})\ge 1-\alpha
$$
This gives contradiction when $\alpha<1/4$.

### 3.2 An Optimal Adaptive Confidence Interval

Set empirical cdf be $F_n(\cdot)$. Standard Gaussian cdf $\Phi(t)=\mathbb{P}[N(0,1)\le t]$. Define
$$
\widehat{CI}=\left[\max_{1.6\le t\le \Phi^{-1}\left(1-\sqrt{\log (2/\alpha)/(2n)}\right)}\left(F_n^{-1}(2(1-\Phi(t)))+t-\frac{2}{t}\right),\\
\quad\quad\quad \min_{1.6\le t\le \Phi^{-1}\left(1-\sqrt{\log (2/\alpha)/(2n)}\right)}\left(F_n^{-1}(1-2(1-\Phi(t)))-t+\frac{2}{t}\right)\right]
$$

**Theorem2** There exists some universal constant $C>0$ s.t. for any $\alpha\in(0,1)$ and $n\ge C\log(1/\alpha)$, the $\widehat{CI}$ defined above satisfies 
$$
\text{(coverage)}\quad \inf_{\epsilon\in[0,0.05],\theta,Q}P_{\epsilon,\theta,Q}(\theta\in\widehat{CI})\ge 1-\alpha
$$

$$
\text{(length)}\quad \inf_{\epsilon\in[0,0.05], \theta,Q}P_{\epsilon,\theta,Q}\left(|\widehat{CI}|\le C'(\frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}}) \right)\ge 1-\alpha
$$

where $C'=C'(\alpha)$.

Thus
$$
r_{\alpha}(\epsilon,[0,0.05])\lesssim \frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}}
$$

### 3.3 Uncertainty Quantification of Quantiles

We should intepret the $\widehat{CI}$ above as
$$
\widehat{CI}=\bigcap_{1.6\le t\le \Phi^{-1}(1-\sqrt{\log(2/\alpha)/(2n)})}[\hat{\theta}_L(t)-\frac{2}{t},\hat{\theta}_R(t)+\frac{2}{t}]
$$
where $\hat\theta_L,\hat\theta_R$ can be viewed as estimators for $\theta$.
$$
\hat\theta_L(t)=F_n^{-1}(2(1-\Phi(t)))+t
$$

$$
\hat{\theta}_R(t)=F_n^{-1}(1-2(1-\Phi(t)))-t
$$

**Lemma1** Simultaneously for $t\in[1.6,\Phi^{-1}(1-\sqrt{\log(2/\alpha)/(2n)})]$, we have
$$
\hat\theta_L-\theta\le \frac{2}{t}\quad and\quad \hat\theta_R-\theta\ge -\frac{2}{t}
$$
with $P_{\epsilon,\theta,Q}$-probability at least $1-\alpha$, uniformly over $\epsilon\in[0,0.05]$, all $\theta\in\mathbb{R}$ and all $Q$.

**Lemma2** For any $\epsilon\in[0,0.05]$, define 
$$
t_{\epsilon}=\Phi^{-1}(1-\epsilon-\sqrt{\log(2/\alpha)/(2n)})
$$
Then we have
$$
\hat\theta_L(t_\epsilon)-\theta\ge 0, \hat{\theta}_R-\theta\le 0
$$
with $P_{\epsilon,\theta,Q}$-probability at least $1-\alpha$ uniformly over $\theta\in\mathbb{R}$ and all $Q$.

Lemma2 leads to the length in theorem2（找到了一个足够短的区间）
$$
\widehat{CI}\le (\hat{\theta}_R(t_\epsilon)+\frac{2}{t_\epsilon}-(\hat{\theta}_L(t_\epsilon)-\frac{2}{t_\epsilon}))\le \frac{4}{t_\epsilon}\asymp\frac{1}{\sqrt{\log n}}+\frac{1}{\sqrt{\log(1/\epsilon)}}
$$

## 4 A Robust Testing Problem

### 4.1 Equivalence between Robust Testing and Robust Confidence Interval

如果一个置信区间必须对 $\epsilon\in[0,\epsilon_{max}]$ 都满足覆盖率要求，那么数据实际污染比例只有 $\epsilon$ 的时候，这个置信区间可以多短？

Start with general contamination $P_{\epsilon,\theta,Q}=(1-\epsilon)P_\theta+\epsilon Q$ with some parametric model $P_\theta$. With $n$ iid observations $X_1,\cdots,X_n$ drawn from  some distribution $P$, consider the following two pairs of robust hypothesis testing
$$
\mathcal{H}(\theta,\theta-r,\epsilon): H_0:P\in\{P_{\epsilon_{max},\theta,Q}:Q\}\quad v.s.\quad H_1:P\in\{P_{\epsilon,\theta-r,Q}:Q\}
$$

$$
\mathcal{H}(\theta,\theta+r,\epsilon): H_0:P\in\{P_{\epsilon_{max},\theta,Q}:Q\}\quad v.s.\quad H_1:P\in\{P_{\epsilon,\theta+r,Q}:Q\}
$$

如果可以高概率相信 $H_0$，就说明对所有的 $\epsilon\in[0,\epsilon_{max}]$，$\theta$ 是更可能的parameter而不是 $\theta-r$，因此置信区间就可以把 $\theta-r$ 排除了。

#### 4.1.1 ARCI via Robust Testing

For each $\epsilon\in[0,\epsilon_{max}]$, let $\phi_{\theta,\theta-r_\epsilon,\epsilon}$ be a binary testing function for $\mathcal{H}(\theta,\theta-r_\epsilon,\epsilon)$.

A CI inverting test procedures $\{\phi_{\theta,\theta\pm r_\epsilon,\epsilon}:\epsilon\in[0,\epsilon_{max}]\}$ is given by (not necessarily be interval)
$$
\widehat{CI}=\{\theta:\phi_{\theta,\theta-r_\epsilon,\epsilon}=\phi_{\theta,\theta+r_\epsilon,\epsilon}=0\ \text{for all}\ \epsilon\in[0,\epsilon_{max}]\}
$$
**Proposition1** Suppose for some $\alpha\in(0,1/4)$, the testing functions $\{\phi_{\theta,\theta\pm r_\epsilon,\epsilon}:\epsilon\in[0,\epsilon_{max}],\theta\in\mathbb{R}\}$ satisfy
$$
\text{(Simultaneous Type-1 error)}\quad \sup_Q P_{\epsilon_{max},\theta,Q}\left(\sup_{\epsilon\in[0,\epsilon_{max}]}\phi_{\theta,\theta\pm r_{\epsilon},\epsilon}=1 \right)\le \alpha
$$

$$
\text{(Type-2 error)}\quad \sup_Q P_{\epsilon,\theta\pm r_\epsilon,Q}\left(\phi_{\theta,\theta\pm r_\epsilon,\epsilon}=0 \right)\le \alpha
$$

for all $\theta\in\mathbb{R}$ and $\epsilon\in[0,\epsilon_{max}]$. Then the CI $\widehat{CI}$ satisfies the coverage property
$$
\inf_{\epsilon\in[0,\epsilon_{max}],\theta,Q}P_{\epsilon,\theta,Q}(\theta\in\widehat{CI})\ge 1-2\alpha
$$
If in addition, $\widehat{CI}$ is an interval, then its length can also be controlled
$$
\inf_{\epsilon\in[0,\epsilon_{max}],\theta,Q} P_{\epsilon,\theta,Q}(|\widehat{CI}|\le 2r_{\epsilon})\ge 1-4\alpha
$$

#### 4.1.2 Robust Testing via ARCI

Given $\widehat{CI}$, construct test functions
$$
\phi_{\theta,\theta-r_\epsilon,\epsilon}=\phi_{\theta,\theta+r_\epsilon,\epsilon}=I(\theta\not\in\widehat{CI})
$$
**Proposition2** 见 theorem1的证明。

If we can show that the robust testing problem $\mathcal{H}(\theta,\theta-r_\epsilon,\epsilon)$ or $\mathcal{H}(\theta,\theta+r_\epsilon,\epsilon)$ can not be solved with small Type1 and Type2 errors, this will exclude the possibility of ARCI with length smaller than $r_\epsilon$.

### 4.2 Gaussian Location Model: Upper Bound

In the special case $P_{\theta}=N(\theta,1), \epsilon_{max}=0.05$, we can show formula
$$
\widehat{CI}=\{\theta:\phi_{\theta,\theta-r_\epsilon,\epsilon}=\phi_{\theta,\theta+r_\epsilon,\epsilon}=0\ \text{for all}\ \epsilon\in[0,\epsilon_{max}]\}
$$
recover the optimal ARCI in section 3.2 with approriate testing functions.

### 4.3 Gaussian Location Model: Lower Bound

**Lemma4** see proof of theorem1.

The result in Lemma4
$$
\inf_{Q_0,Q_1} TV(P_{0.05,\theta,Q_0}^{\otimes n},P_{\epsilon,\theta\pm r,Q_1}^{\otimes n})\le \alpha
$$
implies the sum of Type-1 and Type-2 errors must be greater than $1-\alpha$ for any test.

（For any test $\varphi$, set RR $A=\{\varphi=1\}$, then sum of two types of error 
$$
P_0(A)+P_1(A^c)=1+P_0(A)-P_1(A)\ge1-TV(P_0,P_1)
$$
）

## 5 General Location Families

之前的结果都依赖Gaussian distribution的quantile。有没有方法做更一般的情形。

### 5.1 A Laplace Example

Consider $P_\theta=Laplace(\theta,1)$ given by density $p_\theta(x)=\frac{1}{2}\exp(-|x-\theta|)$.

Recall the definition of the optimal ARCI length
$$
r_{\alpha}(\epsilon, \mathcal{E})=\inf\left\{r\ge 0: \inf_{\widehat{CI}\in\mathcal{I}_{\alpha}(\mathcal{E})}\sup_{\theta,Q}P_{\epsilon,\theta,Q}(|\widehat{CI}|\ge r)\le \alpha \right\}
$$
**Proposition3** There exists some constant $c>0$ s.t. $$
$$
r_{\alpha}(\epsilon,[0,0.05])\ge c
$$
for any $\epsilon\in[0,0.05]$ and any $\alpha\in(0,1/4)$.

也就是未知 $\epsilon$ 的时候不存在长度趋于 $0$ 的CI。

不想看了，那不就是否定了一般的adaptive estimator存在的可能性吗。。。

然后theorem5对各种 $P_\theta$ 算了对应 $r_{\alpha}$。


















