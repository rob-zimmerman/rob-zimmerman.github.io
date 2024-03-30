---
title: 'Non-Integer Moments and Distributions'
date: 2024-03-29
permalink: /posts/2024-non-integer-moments
excerpt: 'Does there exist a distribution which is determined only by its non-integer moments? To put it another way, for $p \geq 0$, do there exist random variables $X$ and $Y$ supported on $(0, \infty)$ such that $\mathbb{E}[X^p] = \mathbb{E}[X^p]$ if and only if $p \not \in \mathbb{N}$?'
tags:
  - probability
---
$\newcommand{\N}{\mathbb{N}}$
$\newcommand{\R}{\mathbb{R}}$
$\newcommand{\E}{\mathbb{E}}$
$\newcommand{\C}{\mathbb{C}}$
$\renewcommand{\P}{\mathbb{P}}$
$\newcommand{\one}[1]{\boldsymbol{1}_{#1}}$

## The Moment Problem

In probability theory, a rigorous definition of expectation is followed almost immediately by a proof of the law of the unconscious statistician; following that, one is shown the definition of moments. One encounters integer moments throughout probability and statistics, while non-integer moments are substantially rarer. A natural question that arises is whether a distribution is characterized by its integer moments. That is, given a sequence of real numbers $\\{\mu_j\\} _{j=1}^\infty$, is there (at most) one distribution $F$ with $\mu_p = \int x^p \, \mathrm{d}F$? This question is called a <i>moment problem</i>.

One learns by example that the answer is no, in general. Work on the moment problem reaches back to Stieltjes in 1894 (who himself invented the term <i>moment</i>), with precursors in Chebyshev and Markov.[^1] Using the theory of continued fractions, Stieltjes himself solved the moment problem for distributions supported on $(0, \infty)$, ultimately showing that it relied on the positivity of the determinants of what we now call [Hankel matrices](https://en.wikipedia.org/wiki/Hankel_matrix) built up from the prescribed moments (in the days before measure theory, Stieltjes' moments were defined by what we know today as Riemann-Stieltjes integrals). The standard counterexample came in 1963, when Chris Heyde asked about the moment problem for "commonly used distributions in statistics" and presented the famous log-normal family:[^2]

> <b>Example 1:</b> If $X$ is a standard log-normal random variable --- that is, if $X$ has density $f(x) = \frac{1}{\sqrt{2\pi}x} \exp\left(-\frac{\log(x)^2}{2}\right)$ for $x \in \R$ --- then for all $\varepsilon \in [-1,1]$, the function $f \cdot (1 + \epsilon \cdot \sin(2\pi \log{x}))$ is also a density with the same moments as $X$.

According to Stoyanov[^3], Stieltjes himself actually introduced that case in the same 1894 paper (in a non-statistical context, of course), so Heyde's example was apparently a rediscovery. Durrett[^4] gives us another ones:

> <b>Example 2:</b> Let $\lambda \in (0,1)$. If $X$ has density $f(x) = \left(\int \exp(-\|x\|^\lambda) \, \mathrm{d}x \right)^{-1} \exp\left(-\|x\|^{\lambda}\right)$ for $x \in \R$, then for all $\varepsilon \in [-1,1]$, the function $f \cdot (1 + \varepsilon \cdot \sin( \tan(\lambda \pi /2) \cdot \|x\|^\lambda \cdot \mathrm{sgn}(x)))$ is also a density with the same moments as $X$.

It turns out that this idea generalizes to <i>Stieltjes classes</i>:[^3] if $f$ is a density and $h$ is a non-zero continuous function taking values in $[-1,1]$ such that $\E_{X \sim f}[h(X) \cdot X^n] = 0$ for all $n \in \N$, then for all $\epsilon \in [-1,1]$, the function $f_\varepsilon = f \cdot (1 + \varepsilon \cdot h)$ is a density with the same moments as $f$. The conditions required of $h$ essentially forces it to be periodic, as it is in the two examples above.


## Non-Integer Moments

Going through the explicit proofs of Examples 1 and 2 in Durrett's book, I was struck that the equivalance of the moments ultimately relies on the fact that we're really talking about <i>integer</i> moments here (in Example 1, the integers play into the sine function; in Example 2, they essentially make $x^n \cdot \varepsilon \cdot h \cdot f$ into an odd function that integrates to $0$). Clearly, non-integer moments are a different story. This observation led me to a kind of inverse question: does there exist a distribution which is determined only by its <i>non-integer</i> moments? To put it another way, for $p \in \R$, do there exist positive random variables $X$ and $Y$ such that $\E[X^p] = \E[X^p]$ if and only if $p \not \in \N$? The answer, as it turns out, is no. In fact, we only need the moments to coincide on some dense subset of $\R$ such as the rationals:

> <b>Proposition:</b> Let $X$ and $Y$ be positive random variables. If $A \subseteq \R$ is dense in $\R$ and $\E[X^p] = \E[Y^p]$ for all $p \in A$, then $\E[X^p] = \E[Y^p]$ for all $p \in \R$.

<i>Proof:</i> Fix $p \in \R$ and let $\\{p_j\\}_{j=1}^\infty \subseteq A$ be a sequence with $p_j \nearrow p$ as $j \to \infty$, which must exist by the density of $A$ in $\R$. We have two cases to consider: either $\E[X^p] = \infty$ or $\E[X^p] < \infty$. The first case is fairly easy to handle, while the second is only a bit trickier.


<u> Case 1: $\E[X^p] = \infty$</u>.

Fix a sequence $\\{q_j\\}_{j=1}^\infty$ such that $q_j \to 0$ and $p + q_j \in A$ for all $j$ (we can do this, again, because $A$ is dense in $\R$). We have that

$$\begin{align*}
\infty &= \E[X^p] &&\\
&= \E[\liminf_{j \to \infty} \, X^{p + q_j}] &\\ 
&\leq \liminf_{j \to \infty} \, \E[ X^{p + q_j}] &&\mbox{by Fatou's lemma}\\
&= \liminf_{j \to \infty} \, \E[ Y^{p + q_j}] &&\mbox{because $p + q_j \in A$}\\
&\leq \liminf_{j \to \infty} \, \E[ Y^p]^{(p + q_j)/p}  &&\mbox{by Hölder's inequality}\\
&= \E[Y^p].
\end{align*}$$

Thus $\E[Y^p] = \infty$, which gives $\E[X^p] = \E[Y^p]$.


<u> Case 2: $\E[X^p] < \infty$</u>.

To begin with, let's show that $\E[Y^p] < \infty$ using essentially same technique as we used above, but in the contrapositive direction. Fix another sequence $\\{r_j\\}_{j=1}^\infty$ such that $r_j \to 0$ and $p - r_j \in A$ for all $j$. Then 

$$\begin{align*}
\E[Y^p] &= \E[ \liminf_{j \to \infty} \, Y^{p - r_j}] &&\\
&\leq \liminf_{j \to \infty} \, \E[Y^{p - r_j}] &&\mbox{by Fatou's lemma}\\
&= \liminf_{j \to \infty} \, \E[X^{p - r_j}]  &&\mbox{because $p - r_j \in A$}\\
&\leq \liminf_{j \to \infty} \, \E[X^p]^{(p - r_j)/p}  &&\mbox{by Hölder's inequality}\\
&= \E[X^p] &&\\
&< \infty. &&
\end{align*}$$

Now, observe that for any $j \in \N$, we have

$$\begin{align*}
X^{p_j} &= X^{p_j} \cdot \one{X \in (0,1]} + X^{p_j} \cdot \one{X > 1} &&\\
&\leq 1 + X^{p_j} \cdot \one{X > 1} &&\mbox{since $x^{p_j} \leq 1$ when $x \in (0,1]$}\\
&\leq 1 + X^p \cdot \one{X > 1} &&\mbox{since $p_j \leq p$ implies $x^{p_j} \leq x^p$ when $x > 1$}\\
&\leq 1 + X^p. &&
\end{align*}$$

Since $\E[1 + X^p] < \infty$, we see that the $X^{p_j}$ are dominated by an integrable random variable. By the dominated convergence theorem, we get $\E[X^{p_j}] \to \E[X^p]$ as $j \to \infty$. The exact same argument applied to the $Y^{p_j}$ shows us that $\E[Y^{p_j}] \to \E[Y^p]$ as $j \to \infty$. Since $\E[X^{p_j}] = \E[Y^{p_j}]$ for all $j$, the uniqueness of limits (in metric spaces!) gives $\E[X^p] = \E[Y^p]$. $\square$

Taking $A = \R^{\geq 0} \setminus \N$ resolves the original question.

The condition $\E[X^p] = \E[Y^p]$ for all $p \in \R$ is very strong. When all moments exist, this is equivalent to saying that $\E[e^{tX}] = \E[e^{tY}]$ for all $t \in \R$ (use the Lambert-$W$ function), which is the same as saying that the moment generating functions (mgfs) of $X$ and $Y$ are equal. Since mgfs characterize distributions --- when the mgfs exist --- this gives us a nice corollary:

> <b> Corollary:</b> Let $X$ and $Y$ be positive random variables. If $A \subseteq \R$ is dense in $\R$ and $\E[X^p] = \E[Y^p] < \infty$ for all $p \in A$, then $X \stackrel{d}{=} Y$. 

This result isn't quite groundbreaking, as Gwo Dong Lin proved something quite a bit stronger in 1992:[^5]

> <b>Theorem:</b> Let $X$ and $Y$ be positive random variables, and suppose there exists some $\alpha > 0$ such that $\E[X^\alpha]$ and $\E[Y^\alpha] < \infty$. Let $\\{s_j\\}_{j=1}^\infty \subseteq (0, \infty)$ be a sequence of distinct numbers such that $s_j \to s \in (0, \alpha)$. If $\E[X^{s_j}] = \E[Y^{s_j}]$ for all $j$, then $X \stackrel{d}{=} Y$.

Lin's proof is quite slick. Here's a sketch of it:

<i>Proof (sketch):</i> The assumptions imply that mgfs of $\log{X}$ and $\log{Y}$ viewed as functions over $\C$ --- that is, the functions $z \mapsto \E[X^z]$ and $z \mapsto \E[Y^z]$ --- are analytic in the strip $S = \\{z \in \C: 0 < \Re(z) < \alpha\\}$ (because $\E[X^\alpha], \E[Y^\alpha] < \infty$), and moreover these functions agree on $S$ as well (because the equality of moments assumption activates the [identity theorem](https://en.wikipedia.org/wiki/Identity_theorem)). A (right)-continuity argument then shows that the functions agree on $\\{z \in \C: \Re(z) = 0\\}$. Replacing $z$ with $it$ for $t \in \R$, we see that the <i>characteristic functions</i> of $\log{X}$ and $\log{Y}$ agree everywhere. Thus $\log{X} \stackrel{d}{=} \log{Y}$, and from that comes $X \stackrel{d}{=} Y$. $\square$


[^1]: Shohat, James Alexander, and Jacob David Tamarkin. The problem of moments. Vol. 1. American Mathematical Society (RI), 1950.
[^2]: Heyde, Chris C. "On a property of the lognormal distribution." Journal of the Royal Statistical Society Series B: Statistical Methodology 25.2 (1963): 392-393.
[^3]: Stoyanov, Jordan. "Stieltjes classes for moment-indeterminate probability distributions." Journal of Applied Probability 41.A (2004): 281-294.
[^4]: Durrett, Rick. Probability: theory and examples. Vol. 49. Cambridge university press, 2019.
[^5]: Lin, Gwo Dong. "Characterizations of distributions via moments." Sankhyā: The Indian Journal of Statistics, Series A (1992): 128-132.

