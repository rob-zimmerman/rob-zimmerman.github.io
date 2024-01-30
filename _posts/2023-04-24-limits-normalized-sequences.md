---
title: 'Limits of Self-Normalized Random Variables'
date: 2023-04-24
permalink: /posts/2023-self-normalized-rvs
excerpt: 'I recently tweeted something very silly (redundant information, I know). The [tweet in question](https://twitter.com/mr_roberts_z/status/1650471367299440641) asked the following:

<blockquote>Let $X_0$ be supported on some nonempty $A \subseteq \mathbb{N}_{>0}$ with $\mathbb{P}(X_0 = k) = p_{0,k}$ and $\mathbb{E}[X_0] < \infty$. For each $n \geq 1$, recursively define $X_n$ on $\mathbb{N}_{>0}$ by $\mathbb{P}(X_n = k) = c_n \cdot k \cdot p_{n-1,k}$, where $c_n$ is a normalizing constant. Then, as $n \to \infty$</blockquote>

...then what?'
tags:
  - probability
---
$\newcommand{\N}{\mathbb{N}}$
$\newcommand{\E}{\mathbb{E}}$
$\renewcommand{\P}{\mathbb{P}}$
I recently tweeted something very silly (redundant information, I know). The [tweet in question](https://twitter.com/mr_roberts_z/status/1650471367299440641) asked the following:

> Let $X_0$ be supported on some nonempty $A \subseteq \N_{>0}$ with $\P(X_0 = k) = p_{0,k}$ and $\E[X_0] < \infty$. For each $n \geq 1$, recursively define $X_n$ on $\N_{>0}$ by $\P(X_n = k) = c_n \cdot k \cdot p_{n-1,k}$, where $c_n$ is a normalizing constant. Then, as $n \to \infty$,

...then what? What happens to these random variables in the pointwise limit? I remember asking myself this years ago when I first learned about discrete random variables, and, given the simplicity of the formulation, I was surprised when a search for an answer yielded nothing.

Anyways, the answer I tweeted, "<i>you will waste a huge amount of time</i>," was perhaps somewhat unsatisfying. This post is an attempt to do a little better.

Of course, the normalizing constant here is simply $\left(\sum_{k \in A} k \cdot p_{n-1,k} \right)^{-1} = \E[X_{n-1}]^{-1}$. So one way to look at this construction is to view the mass that $X_n$ places on $k$ as the corresponding (normalized) summand in the expectation of $X_{n-1}$. Observing that $p_{1,k} = \E[X_0]^{-1} \cdot k \cdot p_{0,k}$, we might suspect that this identity holds more generally, and in fact an easy induction argument shows this is true:

> <b>Proposition:</b> For all $n \geq 1$, we have $p_{n,k} = \E[X_0^n]^{-1} \cdot k^n \cdot p_{0,k}$.

<i>Proof:</i> The base case $n=1$ is obvious. Now, assume that $p_{n-1,k} = \E[X_0^{n-1}]^{-1} \cdot k^{n-1} \cdot p_{0,k}$. Then, simply plugging this in, cancelling out normalizing constants and manipulating gives us

$$\begin{align*} 
p_{n,k} &= \left( \sum_{j \in A} j \cdot p_{n-1, j} \right)^{-1} \cdot k \cdot p_{n-1, k} \\
&= \left( \sum_{j \in A} j \cdot \left[\E[X_0^{n-1}]^{-1} \cdot j^{n-1} \cdot p_{0,j} \right] \right)^{-1} \cdot k \cdot \E[X_0^{n-1}]^{-1} \cdot k^{n-1} \cdot p_{0,k} \\
&= \left( \sum_{j \in A} j \cdot \left[ j^{n-1} \cdot p_{0,j} \right] \right)^{-1} \cdot k  \cdot k^{n-1} \cdot p_{0,k} \\
&= \left( \sum_{j \in A} j^{n} \cdot p_{0,j} \right)^{-1} \cdot  k^{n} \cdot p_{0,k} \\
&= \E[X_0^n]^{-1} \cdot k^n \cdot p_{0,k},
\end{align*}$$

as desired. $\square$

We now return to the original question: what can we say about $X_n$ as $n \to \infty$? The answer to this, as it turns out, depends solely on whether $A$ is infinite or not. If $A$ is infinite, then for any fixed $k \in A$ we have

$$\begin{align}
p_{n,k} &= \frac{k^n \cdot p_{0,k}}{ \sum_{j \in A} j^n \cdot p_{0,j}} \nonumber \\
&= p_{0,k} \cdot \left(\sum_{\substack{i \in A \\ i < k}} \left(\frac{i}{k} \right)^n \cdot p_{0,i} + 1 +  \sum_{\substack{j \in A \\ j > k}} \left(\frac{j}{k} \right)^n \cdot p_{0,j}\right)^{-1} \nonumber\\
&< p_{0,k} \cdot \left(\sum_{\substack{i \in A \\ i < k}} \left(\frac{i}{k} \right)^n \cdot p_{0,i} + 1 +  \underbrace{\left(\frac{k+1}{k} \right)^n \cdot \sum_{\substack{j \in A \\ j > k}} p_{0,j}}_{\to \, \infty}\right)^{-1} \label{eq:ineq}\\
& \xrightarrow{n \to \infty} 0 \nonumber.
\end{align}$$

That is, $\P(X_n = k) \xrightarrow{n \to \infty} 0$. Observe that that this works because 

$$ \begin{equation*} \sum_{\substack{j \in A, \\ j > k}} p_{0,j}> 0, \end{equation*}$$

which follows from our assumption that $A$ is infinite. Thus, for any fixed $k \in A$ we have $\P(X_n \leq k) \xrightarrow{n \to \infty} 0$, or equivalently $\P(X_n > k) \xrightarrow{n \to \infty} 1$, and we see that $X_n$ does not converge to a random variable as $n \to \infty$. 

If $A$ is finite, then the situation is quite different. In this case, let $k^* = \max A$. Then we find that 

$$\begin{align*}
p_{n,k^*} &=  \left( \sum_{\substack{i \in A \\ i \neq k^*}} \left(\frac{i}{k^*} \right)^n \cdot \frac{p_{0,i}}{p_{0,k^*}} + 1 \right)^{-1}\\
&= \left( \underbrace{\sum_{\substack{i < k^*}} \left(\frac{i}{k^*} \right)^n \cdot \frac{p_{0,i}}{p_{0,k^*}}}_{\to \, 0} + 1 \right)^{-1}\\
&\xrightarrow{n \to \infty} 1.
\end{align*}$$

That is, $\P(X_n = k^* ) \xrightarrow{n \to \infty} 1$. We can also see this with an approach from the other direction: for any $k \in A$ with $k < k^* $, then we return to the situation of $\eqref{eq:ineq}$: 

$$ \begin{equation*} p_{n,k} \leq \left(\sum_{\substack{i \in A \\ i < k}} \left(\frac{i}{k} \right)^n \cdot \frac{p_{0,i}}{p_{0,k}} + 1 +  \left(\frac{k+1}{k} \right)^n \cdot \sum_{\substack{j \in A \\ j > k}} \frac{p_{0,j}}{p_{0,k}}\right)^{-1} \xrightarrow{n \to \infty} 0. \end{equation*} $$

While the second sum in the parentheses is finite this time, we can be sure it's non-empty because it must include a term for $k^*$ (the inequality is not strict this time!).

We can summarize our answer to the silly tweet question in a neat statement:

<blockquote><b>Theorem:</b> Let $X_0$ be supported on some non-empty $A \subseteq \N^{>0}$ with $\P(X_0 = k) = p_{0,k}$ and $\E[X_0] < \infty$. For each $n \geq 1$, recursively define $X_n$ on $\N^{>0}$ by $\P(X_n = k) = c_n \cdot k \cdot p_{n-1,k}$, where $c_n$ is a normalizing constant. Then, as $n \to \infty$, the following holds:<br><br>

1) If $|A| = \infty$, then $X_n$ diverges.<br> 
2) If $|A| < \infty$, then $X_n$ converges to a point mass at $\max A$. </blockquote>