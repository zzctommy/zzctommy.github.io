---
date: 2026-11-22
title: Improving randomized controlled trial analysis via data-adaptive borrowing
---

[论文地址](https://arxiv.org/abs/2306.16642)

## 1 Introduction

RCT在临床研究中被认为是“黄金标准”，但是往往很贵，而且未必道德（比如给重病患者分配了不治疗）。

想用Observational study进行一些改进。

## 2 Methodology

### 2.1 Notation, assumptions, and objectives

$\mathcal{R}$: RCT；$\mathcal{E}$: external data.

$N_\mathcal{R},N_\mathcal{E}$ be the size resp., $N=N_\mathcal{R}+N_\mathcal{E}$.

$N_t,N_c$ be the size of treat/control in RCT. Assume all $\mathcal{E}$ receives control.

data: $(X,A,Y,R)$, $R=1$ means in RCT, otherwise means in external data.

$\mu_a=\mathbb{E}[Y(a)|R=1]$, ATE $\tau=\mu_1-\mu_0$.

**Assumption 1** 

(i) $Y=AY(1)+(1-A)Y(0)$.

(ii) $Y(0),Y(1)\perp A|(X,R=1)$ (RCT内unconfoundedness)

(iii) propensity score is nontrivial: $0<\pi_A(x)=\mathbb{P}[A=1|X=x,R=1]<1$.

**Assumption 2**

(i) $\mathbb{E}[Y(0)|X=x,R=0]=\mathbb{E}[Y(0)|X=x,R=1]$.

(ii) $\mathbb{P}[R=1,X=x]>0$ if $\mathbb{P}[R=0,X=x]>0$.

### 2.2 Semiparametric efficient estimation under the ideal situation

事已至此先去学半参理论吧。读不了一点。











