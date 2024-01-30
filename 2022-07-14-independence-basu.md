---
title: 'Self-Independence by Basu'
date: 2023-04-24
permalink: /posts/2023-self-independence-via-Basu
tags:
  - probability
---
$$\newcommand{\N}{\mathbb{N}}$$
$$\newcommand{\E}{\mathbb{E}}$$
$$\newcommand{\R}{\mathbb{R}}$$
$$\renewcommand{\P}{\mathbb{P}}$$
Back in 2020, I taught STA261 for the first time. In the midst of the pandemic, some exam questions were a bit overkill. The first iteration of one was given by...

The question was toned down substantially before going live. However, in subsequent discussions with my pals Yanbo and Michaël, we noticed a connection to a famous result in probaiblity theory:

> <b>Theorem:</b> A random variable $Y$ is independent of itself if and only if $Y$ is almost surely constant.

One direction is quite trivial. The other direction is substantially more difficult and is usually proved by appealing to Kolmogorov's zero-one law (of course one easily shows that $$\P(Y \in A) \in \{0,1\}$$ for any event $$A$$ --- perhaps specializing to $$\P(Y = y) \in \{0,1\}$$ for any $$y \in \R$$ --- but clearly this is not enough to immediately conclude that $$\P(Y = y) =1 $$ for some $$y \in \R$$). It turns out that Basu's theorem from statistical inference makes this easy!

> <b>Theorem (Basu):</b> A complete sufficient statistic is independent of every ancillary statistic.

The proof of Basu's theorem in the discrete case is quite beautiful and is easily understood by students armed with a first course in elementary probability theory. The proof in the general case is substantially more complicated due to a number of measure-theoretic hurdles (see Shao, 1999). However, it does <i>not</i> rely on Kolmogorov's zero-one law and as far as I can tell, the argument is not circular.