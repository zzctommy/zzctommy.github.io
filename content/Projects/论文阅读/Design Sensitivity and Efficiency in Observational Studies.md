---
title: Design Sensitivity and Efficiency in Observational Studies
date: 2026-10-25
---

## 阅读随记

[论文地址](https://doi.org/10.1198/jasa.2010.tm09570)

----

### 1. INTRODUCTION

直观上：

+ 研究喝咖啡与心肌梗死的联系，在未探测的协变量下，treatment分布有一点不均匀就可能导致错误因果
+ 研究大量吸烟与肺癌的联系，需要在未探测的协变量下，treatment分布非常不均匀才可能导致错误因果

作者拿出一个例子，研究苯对基因的危害，协变量只有酒精摄入和吸烟。想要问这个例子更偏向哪一种。

78个被试者被分为39组，每组有相同的协变量，但是每组有一个人接触苯的时间更长。

### 2. REVIEW

#### 2.1 Notations

$I$ pairs of two subjects, one treated $Z_{ij}=1$, one controlled $Z_{ij}=0$. So $Z_{i1}+Z_{i2}=1$. Assume matched for measured covariates: $x_{i1}=x_{i2}$ but may exist unmeasured covariates $u_{i1}\not=u_{i2}$. Potential outcome $r_{Tij}$ is assigned to treated group and $r_{Cij}$ is assigned to controlled group, then the observed response is $R_{ij}=Z_{ij}r_{Tij}+(1-Z_{ij})r_{Cij}$, treatment effect is $r_{Tij}-r_{Cij}$. Fisher's sharp hypothesis $H_0:r_{Tij}=r_{Cij}, \forall i,j$. Write $\mathcal{F}=\{(r_{Tij},r_{Cij},x_{ij},u_{ij}):i=1,2,\cdots,I,j=1,2\}$. $\Omega:=\{\text{valid}\ (Z_{11},Z_{12},\cdots,Z_{I2}):Z_{i1}+Z_{i2}=1\}, |\Omega|=2^I$. Write $\mathcal{Z}$ be the event $\mathbf{Z}\in\Omega$.

注意这里的潜在结果不再是二元分布而是一个实数，$\mathcal{F}$ 理解为包含潜在结果和所有协变量的信息集合。

If $r_{Tij}-r_{Cij}\equiv\tau$, define the difference of observed response in i-th pair $Y_i=(Z_{i1}-Z_{i2})(R_{i1}-R_{i2})=\tau+\epsilon_i$ where $\epsilon_i=(Z_{i1}-Z_{i2})(r_{Ci1}-r_{Ci2})$.（保证始终是处理者减对照者）

$Y_i>0$ 的可能来源：(i) $\tau>0$ (ii) $r_C$ 更高的人更容易接受处理，$\epsilon>0$（选择偏差）

敏感性分析关心 $Y_i>0$ 能否只用第二种机制解释。

#### 2.2 Randomization Inference in  Randomized Experiments

RCT下，$\mathbb{P}[Z_{ij}=1|\mathcal{F},\mathcal{Z}]=1/2,\mathbb{P}[\mathbf{Z}=\mathbf{z}|\mathcal{F},\mathcal{Z}]=2^{-I}, \forall \mathbf{z}\in \Omega$.

Under Fisher’s sharp null hypothesis $H_0$, $Y_i=\epsilon_i=(Z_{i1}-Z_{i2})(r_{Ci1}-r_{Ci2})$: $\tau\equiv 0$.

Let $(I+1)q_i$ be the rank of $|Y_i|$(rank:从小到大排序第几个，编号从1到 $I$)， $sgn(y)=1_{y>0}$.

signed rank statistic $T=\sum_{i=1}^I sgn(Y_i)\varphi(q_i)$, where $\varphi:[0,1]\to \mathbb{R}_{\ge 0}$ is a score function.

$T$ 衡量处理者结果比对照者高的证据有多强，sgn统计处理组更高的部分，$\varphi$ 加权。

Wilcoxon’s signed rank statistic: $\varphi(q_i)=q_i$. Sign test: $\varphi(q_i)=1$.

Under RCT and $H_0$, we can calculate the null distribution of $T|\mathcal{F},\mathcal{Z}$ : $Y_i=|r_{Ci1}-r_{Ci2}|$ is fixed under $\mathcal{F}$, $sgn(Y_i)$ has exactly $1/2$ probability to be $1$ and $1/2$ to be $0$ if $Y_i\not=0$. Set $s_i=1_{Y_i\not=0}$ we get $\mathbb{E}[T|\mathcal{F},\mathcal{Z}]=\frac{1}{2}\sum_i s_i\varphi(q_i), Var[T|\mathcal{F},\mathcal{Z}]=\frac{1}{4}\sum_i s_i\varphi(q_i)^2$.

#### 2.3 Sensitivity Analysis in Observational Studies

敏感性分析要求知道偏离RCT $\mathbb{P}[Z_{ij}=1|\mathcal{F},\mathcal{Z}]=1/2$ 的程度(for all $i,j$)

simple model:在配对之前人们是否被治疗是被分别独立以某个概率assign的：$\pi_{ij}=\mathbb{P}[Z_{ij}=1|\mathcal{F}]$. 满足，如果两个人的已观测协变量相同 $x_{ij}=x_{i'j'}$ 那么odds 的比值被 $\Gamma$ 控制：$\Gamma^{-1}\le \frac{\pi_{ij}(1-\pi_{i'j')}}{\pi_{i'j'}(1-\pi_{ij})}\le \Gamma$. 

这个假设等价于($\gamma=\log\Gamma\ge 0,\mathbf{u}\in [0,1]^{2I}$,大概技巧是给odds取个min，每个odds写成$\min\cdot e^{\gamma u}$算出来的）
$$
\mathbb{P}[\mathbf{Z}=\mathbf{z}|\mathcal{F},\mathcal{Z}]=\frac{\exp(\gamma\mathbf{z}^T\mathbf{u})}{\sum_{\mathbf{b}\in\Omega}\exp(\gamma\mathbf{b}^T\mathbf{u})}=\prod_{i=1}^I\frac{\exp(\gamma z_{i1}u_{i1}+\gamma z_{i2}u_{i2})}{\exp(\gamma u_{i1})+\exp(\gamma u_{i2})}
$$
这里的 $\mathbf{u}$ 刻画了被治疗概率倾向的差异。

$$
\mathbb{P}[Z_{ij}=1|\mathcal{F},\mathcal{Z}]=\exp(\gamma u_{ij})/(\exp(\gamma u_{i1})+\exp(\gamma u_{i2}))\in [\frac{1}{1+\Gamma},\frac{\Gamma}{1+\Gamma}]
$$

每次给一个 $\Gamma$ 我们可以去算结果的 $p$-value，如果 $p$-value很小就说明在偏差 $\Gamma$ 下的unmeasured covariates并不能影响结论（结论对隐藏协变量的敏感程度）

定义 $\bar{\bar{T}}_{\Gamma}=\sum_{i}sgn(Y_i)\cdot\begin{cases}\varphi(q_i), s_i=1, probability \frac{\Gamma}{1+\Gamma}\\ 0, s_i=1,probability \frac{1}{1+\Gamma}\\ 0, s_0=0\end{cases}$

$\bar{T}_{\Gamma}=\sum_i sgn(Y_i)\cdot\begin{cases}\varphi(q_i), s_i=1, probability \frac{1}{1+\Gamma}\\ 0, s_i=1,probability \frac{\Gamma}{1+\Gamma}\\ 0, s_0=0\end{cases}$

Then $\mathbb{P}[\bar{T}_{\Gamma}\ge k|\mathcal{F},\mathcal{Z}]\le \mathbb{P}[T_{\Gamma}\ge k|\mathcal{F},\mathcal{Z}]\le \mathbb{P}[\bar{\bar{T}}_{\Gamma}\ge k|\mathcal{F},\mathcal{Z}]$, $\forall \mathbf{u}\in [0,1]^{2I}$.

Add mild conditions on score function $\varphi$ we can use normal approximation to estimate $\mathbb{P}[\bar{\bar{T}}_{\Gamma}\ge k|\mathcal{F},\mathcal{Z}]$。

$T$ 的null distribution是未知的，用已知的 $\bar{T},\bar{\bar{T}}$ 去比较尾概率（算p-value）

#### 2.4 Sensitivity Analysis in the Benzene Example

检验 $H_0:$ 苯对人的基因没有损伤。计算单侧p-value的上界尝试否定 $H_0$：$\mathbb{P}[\bar{\bar{T}}_{\Gamma}\ge T|\mathcal{F},\mathcal{Z}]$, $T$ 由数据以及选择的权重直接算出来。

然后作者算了若干不同统计方法在不同 $\Gamma$ 下给出的 p-value，发现不同的方法 p=0.05的 $\Gamma$ 各不相同，每种统计量都会报告一个敏感界 $\Gamma$，作者想研究如何选择统计方法。

**作者指出，有一些统计量给较小的 $|Y_i|$ 赋予了  $0$ 权，例如 Brown 统计量  $\varphi([0,\frac{1}{3}])=0,\varphi([(\frac{1}{3},\frac{2}{3}])=1,\varphi((\frac{2}{3},1])=2$，可以抵御更大的隐藏偏差。 **

#### 2.5 The Power of a Sensitivity Analysis and  Design Sensitivity

定义power需要先给定数据（$\mathcal{F},\mathbf{Z}$）的生成方式。在观察性研究中，哪怕covariates已经做到unconfoundedness了，我们也不可能断定这一点，因此最好的结论是unconfoundedness下可以得出可以抵御一定程度的隐藏偏差。

有利情形（一种数据生成方式）：$Z_{i1}-Z_{i2}=\pm 1$ with probability $1/2$, $\mathcal{F}$ 按照某个模型生成，且additive：$\tau\equiv r_{Tij}-r_{Cij}$. 再加上其他一些条件。

定义power $\mathbb{P}[\bar{p}_{\Gamma}\le \alpha|\text{有利情形},H_0]$（$\bar{p}_{\Gamma}$：给定 $\Gamma$ 下最坏的 p-value，也就是 $R$ 更高的人总是以最高概率 $\frac{\Gamma}{1+\Gamma}$ 被处理时的p-value）

在确定的统计量和有利情形下，定义 $\tilde{\Gamma}$ 使得在 $I\to\infty$ 时，如果 $\Gamma<\tilde{\Gamma}$ 就有 power $\to 1$ ；如果 $\Gamma>\tilde{\Gamma}$ 就有 power $\to 0$.（存在性暂时未知）

### 3 Intuition

尝试给出“更不关注小的 $|Y_i|$ 导致 less sensitive”的直观。

如果 $H_0$ 是对的，给定 $\Gamma$ 能否用未观测偏差 来解释 得到的结果？强度为 $\Gamma$ 的未观测偏差在 $H_0$ 下至多把出现正 $Y_i$ 的概率推高到 $\frac{\Gamma}{1+\Gamma}$。设RCT下 $Y_i$ 的density为 $g(\cdot)$，定义 $\operatorname{abz}(y)=\frac{g(y)}{g(y)+g(-y)}$，也就是给定 $|y|$ 的时候取正号的概率。如果这个值比 $\frac{\Gamma}{1+\Gamma}$ 大，那么就不是强度 $\Gamma$ 的未观测偏差能解释的了，因而这部分数据是非常重要的，需要调节权重使他们突出出来。

作者计算了若干种分布的 abz（$Y_i=\frac{3}{4}+\epsilon_i$, $\epsilon_i$ 的分布人为给定），有一些会先增后降（但是无一例外在 $Y_i$ 非常小的时候都是小的），所以并非越大的 $|Y_i|$ 配越大的权重是一个好策略，定量分析在下一节。

### 4 DESIGN SENSITIVITY AND EFFICIENCY

#### 4.1 A Formula for the Design Sensitivity With General Scores

我们给已经定义的变量加下标 $I$ 表示根据已有的 $I$ 对数据计算得到。

定义 score function $\varphi(\cdot)$ 是 simple 的，如果 (i) $\varphi:[0,1]\to[0,\infty)$; (ii) $\varphi$ piecewise continuous and right-continuous, i.e. exists $0=\xi_0<\xi_1<\cdots<\xi_A=1$ s.t. $\varphi$ continuous at $[\xi_{i-1},\xi_{i}),i\le A-1$ and $[\xi_{A-1},1]$； (iii) 非严格递增; (iv) 不恒为零。

Fix $\alpha=0.05$, write $\Phi$ for standard normal cdf. Define
$$
c_{\Gamma,I}=\frac{\Gamma}{1+\Gamma}\sum_{i=1}^I s_i\varphi(q_{i,I})+\Phi^{-1}(1-\alpha)\sqrt{\frac{\Gamma}{(1+\Gamma)^2}\sum_{i=1}^I s_i\varphi(q_i)^2}
$$
**Proposition**: If $Y_i,i=1,\cdots,I$ are iid with continuous density $G(\cdot)$. $\varphi(\cdot)$ is a simple score function. Then for a.s. $\mathcal{F}_I$ as $I\to \infty$ we have
$$
\mathbb{P}[\bar{\bar{T}}_{\Gamma,I}\ge c_{\Gamma,I}|\mathcal{F}_I,\mathcal{Z}_I]\to \alpha\\
\mathbb{P}[T_{I}\ge c_{\Gamma,I}|\mathcal{F}_I,\mathcal{Z}_I]\to 1\qquad \text{if}\ \Gamma<\tilde{\Gamma}\\
\mathbb{P}[T_{I}\ge c_{\Gamma,I}|\mathcal{F}_I,\mathcal{Z}_I]\to 0\qquad \text{if}\ \Gamma>\tilde{\Gamma}\\
$$
where 

$$
\tilde{\Gamma}=\frac{\int_0^\infty \varphi(G(y)-G(-y)) dG(y)}{\int_0^1 \varphi(y)dy-\int_0^\infty \varphi(G(y)-G(-y))dG(y)}
$$

pf: 连续的density导致我们可以只考虑 $|Y_i|$ 互不相同的情况（不影响a.s.性）

从而有 $\sum_{i=1}^I s_i\varphi(q_{i,I})=\sum_{i=1}^I\varphi(i/(I+1))$ 和 $\sum_{i=1}^I s_i\varphi(q_{i,I})^2=\sum_{i=1}^I\varphi(i/(I+1))^2$.

由$\varphi$ simple的假设每一个求和都是 $O(I)$ 的阶（无界）。

$\bar{\bar{T}}$ 是 $I$ 个独立随机变量的和，每一个有 $\frac{\Gamma}{1+\Gamma}$ 的概率是 $\varphi(i/(I+1))$，$\frac{1}{1+\Gamma}$ 概率是 $0$。

计算mean和var就可以得到 $\mathbb{P}[\bar{\bar{T}}_{\Gamma,I}\ge c_{\Gamma,I}|\mathcal{F}_I,\mathcal{Z}_I]\to \alpha$.

Now consider 
$$
D_I=\frac{\frac{1}{I}T_I-\frac{1}{I}\frac{\Gamma}{1+\Gamma}\sum_{i=1}^I s_i\varphi(q_{i,I})}{\sqrt{\frac{\Gamma}{I^2(1+\Gamma)^2}\sum_{i=1}^I s_i\varphi(q_{i,I})^2}}
$$
note $T_I\ge c_{\Gamma,I}$ is equivalent to $D_I\ge\Phi(1-\alpha)$.

denominator of $D_I$ is $O(1/\sqrt{I})$ tends to $0$.

the second term in numerator tends to $\frac{\Gamma}{1+\Gamma}\int_0^1 \varphi(x)dx$. 

the first term in numerator
$$
\frac{T_I}{I}=\frac{1}{I}\sum_{i=1}^Isgn(Y_i)\varphi(q_{i,I})\to\int_0^\infty \varphi(G(y)-G(-y))dG(y)
$$
直观上，$G(y)-G(-y)=\mathbb{P}[|Y_i|\le y]$. 在 $[0,\infty)$ 上积 $dG$ 意为 $Y_i>0$. $\varphi$ 的间断点由于有限会被忽略。

所以numerator of $D_I$ 的正负性决定了 $D_I$ 的极限 $\pm \infty$，反解临界情形得到 $\tilde{\Gamma}$.

#### 4.2 Design Sensitivities of Several Statistics

作者计算了Table3，其中Pitman ARE指的是两个检验 $A,B$ 要达到相同power需要的样本量比值：$H_0:\lambda=0$. 当真实的 $\lambda>0$ 时要达到power $\alpha$ 需要的样本量记为 $n_A(\lambda),n_B(\lambda)$. $ARE(A,B)=\lim_{\lambda\to 0}\frac{n_A(\lambda)}{n_B(\lambda)}$.一般认为ARE小于1意味着相同样本量下 $A$ 的power更大。

表格里 $B$ 恒为Wilcoxon。

其他的量我都见过。

结论有：

+ Wilcoxon的 $\tilde{\Gamma}$ 最小，敏感度最高。
+ 损失很少的 Pitman ARE可以大幅提高 $\tilde{\Gamma}$.
+ 降低 $\operatorname{abz}(y)$ 小的部分的权重可以提高 $\tilde{\Gamma}$.
+ 如果样本无限多，那么应该只考虑 $\operatorname{abz}$ 最大的部分，但是实际样本有限，不能只考虑很窄的一段区域，引出下一节的有限样本分析。

#### 4.3 Simulated Power of a Sensitivity Analysis

Table4,5算的是有限样本下，用 $Y_i=\tau+\epsilon_i$ 生成数据（$\tau\not=0$, 也就是 $H_1$） 敏感性分析拒绝 $H_0$ 的概率（power的定义）。

结论有：

+ Table3,4,5一起看会发现，每一种情形中存在Pitman ARE相近的统计量，$\tilde{\Gamma}$ 不相同，而且 $\tilde{\Gamma}$ 的大小顺序的确给出了power的顺序，所以 $\tilde{\Gamma}$ 某种意义上比 Pitman ARE更能反映power。因此可以恰当选择统计量来同时拥有高Pitman ARE和高power。$\color{red}\text{并没有告诉人们怎么选统计量。}$
+ 应当避免给 $\operatorname{abz}$ 小的 $Y_i$ 太多的权重，会降低敏感性分析的power。
+ 大的 $Y_i$ 是否有用要根据abz选，但是实际操作时并不知道分布，abz不能预先知道，作者在本文没有给出解决方法，引用了他的一篇文章，方向是从所有数据中分出一小部分来指导其余部分的选择。

### 5. DISCUSSION

1. 讨论的起点是认真考虑非随机化分配带来的偏差，形式化其对因果的影响。
2. 研究者可以自由选择权重来达到目标
3. 不必局限于 rank statistics $T$。说是有个什么 Huber-M statistics，我不知道是什么。
4. 不必局限于 $\tau$ 恒为常数。另一个有趣的情形是只有少部分被处理的人有强烈反应而其他人都没有变化。

### 6. SUMMARY

假设unconfoundedness在观察性研究中是不安全的，除去假设unconfoundedness的分析应该加入敏感性分析。$\tilde{\Gamma}$ 在很多统计量上可以比Pitman ARE给出更好的与敏感性分析power的关系，尽管二者都是与power有关的。















 