---
title: 'Self-Independence by Ancillarity and Completeness'
date: 2022-07-14
permalink: /posts/2022-self-independence-ancillarity
tags:
  - probability
excerpt: 'Back in 2020, I taught [STA261](https://rob-zimmerman.github.io/teaching/STA261) for the first time. The first part of that course deals with statistics (i.e., functions of random samples, not the subject as a whole!) and I chose to provide a light introduction to completeness because of how elegant the [Lehmann-Scheffé theorem](https://en.wikipedia.org/wiki/Lehmann%E2%80%93Scheff%C3%A9_theorem) and related results in point estimation are down the road...'
---
$$\newcommand{\N}{\mathbb{N}}$$
$$\newcommand{\E}{\mathbb{E}}$$
$$\newcommand{\R}{\mathbb{R}}$$
$$\newcommand{\bX}{\mathbf{X}}$$
$$\renewcommand{\P}{\mathbb{P}}$$
Back in 2020, I taught [STA261](https://rob-zimmerman.github.io/teaching/STA261) for the first time. The first part of that course deals with statistics (i.e., functions of random samples, not the subject as a whole!) and I chose to provide a light introduction to completeness because of how elegant the [Lehmann-Scheffé theorem](https://en.wikipedia.org/wiki/Lehmann%E2%80%93Scheff%C3%A9_theorem) and related results in point estimation are down the road, despite the unintuitive definition of completeness itself (I've kept up this choice in my subsequent offerings of the course). By the same token, I also introduced ancillarity in order to apply the extremely slick [Basu's theorem](https://en.wikipedia.org/wiki/Basu%27s_theorem) to several nice problems.

The first assessment in that course (in the form of an online quiz --- this was during lockdown, after all) was intended to test the students' understanding of these concept; with lots of difficulty, I managed to come up with several original questions. One of them was meant to be on the tricker side:

> Let $X_1,\ldots,X_n \stackrel{iid}{\sim} F_\theta$, where the $X_i$'s have finite first moments. Suppose that there exist some $j,k \in \{1,\ldots,n\}$ such that the statistic $S(\bX) = X_j$ is ancillary for $\theta$ and the statistic $T(\bX) = X_k$ is complete for the family $\\{F_\theta: \theta \in \Theta\\}$. Prove that all of the $X_i$'s must be constant with probability 1.

<i>Solution:</i> Let $h(T(\bX)) = T(\bX) - \E[X_j]$, which is free of $\theta$ because $X_j$ is ancillary for $\theta$. Then for any $\theta \in \Theta$, we have $\E[h(T(\bX))] = \E[X_k] - \E[X_j] = 0$, and by completeness it follows that $1 = \P_\theta( h(T(\bX)) = 0) = \P_\theta(X_k = \E[X_j])$ for all $\theta \in \Theta$. That is, $X_k$ is the constant $\E[X_j]$ with probability 1. Since the $X_i$'s are iid, they're all equal to that same constant with probability 1. $\square$

The question originally had a number of red herrings thrown in, but fortunately I removed those before going live. In a subsequent discussion with pals [Yanbo](https://yanbotang.github.io/) and [Michaël](https://mic-lalancette.github.io/), we noticed a connection to a basic result in probablity theory:

> <b>Theorem:</b> A random variable $Y$ is independent of itself if and only if $Y$ is almost surely constant.

One direction of this is quite trivial. The other direction is substantially more difficult and is typically proved by appealing to [Kolmogorov's zero-one law](https://en.wikipedia.org/wiki/Kolmogorov%27s_zero%E2%80%93one_law) (of course one easily shows that $$\P(Y \in A) \in \{0,1\}$$ for any event $$A$$ --- perhaps specializing to $$\P(Y = y) \in \{0,1\}$$ for any $$y \in \R$$ --- but clearly this is not enough to immediately conclude that $$\P(Y = y) = 1 $$ for some $$y \in \R$$, and some topological argument is required). However, if we're willing to further assume that $Y$ has a finite first moment, then we can produce an elementary proof using the result from the STA261 quiz:

$(\Rightarrow)$: Suppose that $Y$ is independent of itself. Define the "parameter space" $\Theta = \\{ \theta \\}$, where $\theta \in \R$ is arbitrary. We can vacuously associate $\theta$ to the distribution of $Y$ and write $Y \sim F_\theta$ without any ambiguity, even though the distribution of $Y$ is free of $\theta$. By the last remark, the "statistic" $S(Y) = Y$ is ancillary for $\theta$. Now, let $h: \R \to \R$ be some function such that $\E[h(Y)] = 0$. Since $Y$ is independent of itself, so too is $h(Y)$ and we have $0 = \E[h(Y)] \cdot \E[h(Y)] = \E[h(Y)^2]$, and it follows[^1] that $\P(h(Y) = 0) = 1$, and this holds vacuously for all $\theta \in \Theta$. Thus, the "statistic" $T(Y) = Y$ is a complete (and obviously sufficient) for the family $\\{ F_\theta: \theta \in \Theta \\}$. By the result from the quiz with $Y$ in place of each $X_i$, we see that $Y$ must be constant. $\square$

Unfortunately, this will not win the "most elegant proof" award for that direction of the theorem under the assumption of a finite first moment: since $Y$ is independent of itself, $\E[Y^2] = \E[Y] \cdot \E[Y] = \E[Y]^2$ so that $\text{Var}(Y) = \E[Y^2] - \E[Y]^2 = 0$ and the result follows.

[^1]: One could reasonably protest that the argument here is not completely elementary because the statement $\P(X \geq 0) = 1 \implies (\E[X] = 0 \iff \P(X) = 0) = 1$ requires some basic measure theory to prove. But it's still much simpler than Kolmogorov's zero-one law!