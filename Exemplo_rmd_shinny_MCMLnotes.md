---
title: "Monte Carlo likelihood approach for finite population theory"
author: "Markus Stein"
date: "11 de outubro de 2018"
output: html_document
runtime: shiny
---



The basic principle in statistical inference is to learn about the probability 
distribution of a random variable from observing a random sample from it, and 
possibly also incomplete populational data. Traditional methods like 
maximum likelihood methods can show a bad performance even for simple problems. 
For example, to estimate the probability of success from a sample of binary 
variables when it is assumed to be close to zero in the population, we might 
expect to end up with a monotone (flat) likelihood... Only for this basic 
problem several adjustments/approximations of the likelihood function were suggested 
in the literature. Our aim is to discuss the likelihood construction for a sample of 
binary variables when the population total is known, and in particular to show how 
we can approximate the true likelihood via Monte Carlo simulations (reference MClikelihood)... 
We are motivated to allow for any complex sampling scheme... this problem is originally ... we also discuss survey-based approaches, simplest case of random 
variables, when it assume a Bernoulli

## Finite population sampling example

Let $Y$ be a binary (usually response) variable and denote $\pi = P(Y=1)$ the 
probability of $y$ being a success. Also assume that a finite population (or cohort) 
$\mathcal{U}=\left\{ 1, \ldots, N\right\}$ is sampled from a hypothetical 
superpoluation. Let $y_u$ be the value of variable $Y$ for the $u$th population 
unit. Very often a inference target in design-based approach is a population 
(non-random) total or a mean, $\bar{Y} = N^{-1} \sum_{u \: \in \: \mathcal{U}} y_u$.

In practice we do not observe all units in the complete finite population, 
instead we only observe a sample $\mathcal{S}$ of size $n$ from that. Assume 
$\mathcal{S} \subseteq \mathcal{U}$ is drawn with a given sampling design, and 
the inclusion probabilities are given by $\pi_u = P(u \in \mathcal{S})$ and 
assumed to be strictly positive. To guarantee a designed-unbiased estimator, 
denote $w_u = \pi_u^{-1}$ the design weight associated to unit $u$, such that 
$\pi_u > 0$ for all $u \in \mathcal{U}$. Terms $w_u$ are usually interpreted 
as the number of population units represented by the $u$th sample unit. For a 
known $N$, the weighted estimator, or Horvitz-Thompson (HT) estimator 
$$\bar{y}_{w} = \frac{1}{N} \sum_{u \: \in \: \mathcal{S}} w_u \: y_u.$$

Properties of the HT estimator follow from randomisation theory. The inverse 
probability weighting (IPW) is a very commom technique in design-based 
estimation methods, and also popular in missing data problems. For a further 
discussion about weighted estimators see \cite{Thompson1997} and \cite{Lohr2010}.

## Probability of success

To compare randomisation approaches with model-based methods consider the same 
problem of estimating a mean, but now our goal is to estimate the complete-data 
likelihood function, which is the product of the joint probability distribution 
for each unit $u$. Let us denote the probability distribution of variable $Y$ as 
$$ f\left( y \right) = P\left(Y=y\right) = (1-\pi)^{(1-y)} \pi^y,$$
where $y \in \{0,1\}$. The logarithm of the distribution $f$ is $\log f(y) = y \log \pi + (1-y) \log (1-\pi)$, and the derivative of function $f$ with respect to $\boldsymbol{\pi}$ is $\frac{\partial}{\partial \boldsymbol{\pi}} f(y) = f(y) \times \frac{\partial}{\partial \boldsymbol{\pi}} \log f(y)$, hence the derivative of the logarithm is $\frac{\partial}{\partial \boldsymbol{\pi}} \log f(y) = \frac{y}{\pi} - \frac{(1-y)}{1-\pi}$.


```r
N <- 15000          # population size
n <- 50             # sample size
p <- 0.05           # probability of being a case
y <- rbinom(N,1,p)  # finite pop. sampled sample from Y ~ Bernoulli(p)
U <- 1:N            # set of indexes for population
S <- sample(U, n)   # sample of size n from U
table(y[S])/n
```

```
## 
## 0 
## 1
```

## Individual-level likelihood function (complete unobserved data)

The complete-data likelihood assumes $Y$ observed for all units in 
$\mathcal{U}$ then the joint distribution for all observations is given by 
$$f\left(\boldsymbol{y} \right) \: = \: \prod_{u \in \mathcal{U}} 
    f\left(y_u \right) \: = \: \prod_{u \in \mathcal{U}} \pi^{y_u} 
    \: = \: \pi^{\sum_{u \in \mathcal{U}} y_u} (1-\pi)^{\sum_{u \in 
    \mathcal{U}} (1-y_u)} \: = \: \pi^{N_1} (1-\pi)^{N-N_1},$$
with $y_u$ being the indicator $Y$ for the $u$th unit. Note that the 
population number of cases $(Y=1)$ can be written as $N_1 = \sum_{u \in \mathcal{U}} y_u$. 
For a given population counts $\boldsymbol{N_{y}} = (N_0, N_1)$, 
$L_u\left( \boldsymbol{\pi} \right) = L_u\left( \boldsymbol{\pi}; \boldsymbol{N_{y}} \right) = f\left( \boldsymbol{y} \right)$. Then the complete-data log-likelihood function can be expressed as 

\begin{eqnarray}
\label{cohort_loglik}
\ell\left(\boldsymbol{\pi} \right) & = & \log f\left(\boldsymbol{y}\right) 
\nonumber \\
        & = & \sum_{u \in \mathcal{U}} \log f\left(y_u\right) \nonumber \\
        & = & N_1 \log \pi + (N-N_1) \log (1-\pi). \nonumber \\
\end{eqnarray}
One may note that the score function for an observed population counts 
$\boldsymbol{N_y}$ can be written as 
$S \left(\boldsymbol{\pi}\right) = \frac{\partial}{\partial \pi} \ell(\pi) = 
\frac{1}{L_t\boldsymbol{\theta}} \: \frac{\partial L_t\left(\boldsymbol{\theta}\right)}
{\partial \boldsymbol{\theta}} = \frac{1}{L_t\boldsymbol{\theta}} 
\prod_{u \in \mathcal{U}} \frac{\partial}{\partial \boldsymbol{\theta}} 
f_{\boldsymbol{\theta}}(y_u)$.


## Weighted (pseudo) likelihood



population \in \mathbb{R}$ Suposse X a population of interest follawing a 




## Likelihood plot simple example
Here we show some simulation results...


```
## Error in loadNamespace(name): there is no package called 'webshot'
```



## Likelihood theory for infinite population
Supose a finite pop $f$, and a simple random sample of size $n$ from $f$, then...  

Following likelihood theory the maximum likelihood estimator is given by $S \left(\boldsymbol{\pi}\right) = \frac{\partial}{\partial \pi} \ell(\pi) = \boldsymbol{0}$, when the maximum is exists...  

Also consider $J \left(\boldsymbol{\pi}\right) = - \frac{\partial}{\partial \pi} S(\pi)$ the observed Fisher information, and the expected Fisher information is $I \left(\boldsymbol{\pi}\right) = E \left[ J \left(\boldsymbol{\pi} \right) \right]$.  

If regularity conditions are satisfied then $I \left(\boldsymbol{\pi}\right) = Var \left[ S\left(\boldsymbol{\pi}\right) \right]$.  

By the weak law of large numbers $ $. Also the central limit theorem holds for .... to compare to clt for finite pop. $\sum_{i=1}^N$ but a sample $n$ is observed... or $\sum_{i=1}^n$ and connect to missing....  

## Parameters of superpopulation and survey population
In Godambe and Thompson (1986???) (.. other references???...) the authors discuss the simultaneous estimation of superpopulation model and survey (finite population) parameters. The estimation equation theory applied to this problem is shown to give optimal, which means design unbiased and minimal joint mean squared error with respect to the sampling design and the model.  

### Superpopulation model
Assume the simple linear model $Y = \beta \: X + \mathcal{E}$, where $X=x$ is a realization of a random variable (...), $E \left( \mathcal{E} \right) = 0$ and $Var \left( \mathcal{E} \right) = \sigma^2$. We call $\beta$ the superpopulation parameter of interest, and assuming that we observe a finite population of 
size $N$, $(x_i, y_i)$ for $i=1, \ldots, N$, denote $\beta_N$ the least square estimation of $\beta$.  

Assume that $\boldsymbol{y} = (y_1, \ldots, y_N)$ is generated from a distribution $ $ where is known to be a member of a class  $C = \left\{ \xi \right\}$. $C$ is then called a *superpopulation model*. Estimating equantion theory theory for $\boldsymbol{\theta}$ can be used to show that, for a given function $g(\boldsymbol{y}, \boldsymbol{\theta})$ for all $g(\cdot) \in C$, $E_{\xi} \left[ g(\boldsymbol{y}, \boldsymbol{\theta}(\xi)) \right]=\boldsymbol{0}$. 

If we have complete population information, then by the law of large numbers $\beta_N$ converges to the superpopulation parameter $\beta$ as $N$ increases. But we do not observe all Elements in $\mathcal{U}$, so the aim is also to estimate $\boldsymbol{\beta_N}$.

### Finite population
Chapter 3 Godambe and Thompson (1986)...



## Mathematical Results 

Obs. 1: $E \left[ e^{ix} \right] = \sum_{k=0}^n \dfrac{(ix)^k}{k!} E[X^k],\hspace{2cm} t \rightarrow 0$  

Obs. 2: (Taylos expansion) (Billingsley 3rd ed., pg. 345) Define a $X$ random variable in a probability space... . If $E \left( X^2 \right)$ is finite, then 
$$ \phi \left( \boldsymbol{t} \right) = 1 + i \boldsymbol{t}^\top E \left( \boldsymbol{X} \right) + \frac{ i^2  }{2} \boldsymbol{t}^\top E \left( \boldsymbol{X}^\top \boldsymbol{X} \right) \boldsymbol{t} + o \left( \boldsymbol{t}^\top \boldsymbol{t} \right),\hspace{4cm} \boldsymbol{t} \rightarrow \boldsymbol{0}.$$
...to check Cramer-Wold device/theorem...

Obs. 3: Lemma (Importance Sampling convergence) Assume that conditions for the Strong Law of Large Numbers (SLLN) holds for a colection of $X_1, X_2, \ldots, X_n$ independent and identically distributed random variables (defined in a continuous space $\mathcal{X}$) from $g$ (coninuous). Also consider a function $f$ and assume that $x_1, x_2, \ldots , x_n$ are  sampled from $g$, so
$$ n^{-1} \sum_{i=1}^n \frac{f(x_i)}{g} \xrightarrow[]{a.s} E \left[ f \right].$$
Proof:   
It follow from the SLLN that 
\begin{eqnarray}
\label{importance_sampling}
n^{-1} \sum_{i=1}^n \frac{f(x_i)}{g} 
    & \xrightarrow[]{a.s} & E \left[ \frac{f(x_i)}{g} \right] \nonumber \\
    & = & \int_{\mathcal{X}} \frac{f(x_i)}{g} g \nonumber \\
    & = & E \left[ f \right].
\end{eqnarray}
The $\mathcal{X}$ discrete case follows straight from the definition of discrete function. 
