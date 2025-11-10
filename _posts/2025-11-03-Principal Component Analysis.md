---
title: Principal Component Analysis
subtitle: Mathematical Deivation
layout: default
date: 2025-11-03
keywords: PCA, dimensionality reduction
published: false
---

In this blog post I will outline the mathematical derivation of principal component analysis (PCA). (Much) more detail can be found in various books, the one I used as the primary source of information is _An Introduction to Multivariate Statistical Analysis_ by T. W. Anderson.

## Derivation

Suppose that we have $X: \Omega \xrightarrow{} \mathbb{R}^p$, i.e. a $p$-dimensional random variable with $Var(X) = \Sigma$. For simplicity, we assume $E(X) = 0$, otherwise we can take a different random variable and define it as $Y = X - E(X)$.

Principal component analysis poses the following problem. What is the linear combination of the $p$-components of the random variable X with maximum variance. That is, which $\beta \in \mathbb{R}^p$ maximizes $Var(\beta^\prime X)$?

By itself, this is an ill-posed problem since we can simply let $\beta$ go to $\infty$. Therefore we induce an additional restriction $\beta^\prime \beta = 1$. Rewritting $Var(\beta^\prime X) = \beta^\prime Var(X) \beta = \beta^\prime \Sigma \beta$ we have the following constrained optimization problem: $$ \max_\beta \beta^\prime \Sigma \beta \hspace{1cm} \text{subject to  } \beta^\prime \beta = 1. $$
Now let $$ \mathscr{L} (\beta, \lambda) = \beta^\prime \Sigma \beta - \lambda (\beta^\prime \beta - 1) $$ where $\lambda$ is a Lagrange multiplier. The vector of derivatives is $$ \frac{\partial \mathscr{L}}{\partial \beta} = 2\Sigma \beta - 2 \lambda \beta. $$ 
Setting this expression to 0 we get $$ (\Sigma - \lambda I)\beta_1 = 0, $$
that is our first principal component $\beta_1$ is an eigenvector of the variance-covariance matrix $\Sigma$ with some corresponding eigenvalue $\lambda_1$. Therefore the objective function is $\beta_1^\prime \Sigma \beta_1 = \lambda_1 \beta_1^\prime \beta_1 = \lambda_1$. Since we wish to maximize this expression, we choose $\lambda_1 = \lambda_{max} (\Sigma)$. We have found the first principal component, which we may denote $Z_1 = \beta_1 X$.

Now let's find a normalized combination $\beta^\prime X$ that has maximum variance and is uncorrelated with $Z_1$. This means that (since $Var(Z_1) = 1, Var(\beta^\prime X) = 1$) $$ Corr(\beta^\prime X, Z_1) = Cov(\beta^\prime X, \beta_1^\prime X) = \beta^\prime Cov(X,X) \beta_1 = \beta^\prime \Sigma \beta_1 = \lambda_1 \beta^\prime \beta_1 = 0. $$

Therefore $\beta^\prime X$ is orthogonal to $Z_1$ in statistical sense as well as a pure geometric one (for $\lambda_1 \neq 0$). Let us set up the second lagrangian
$$ \mathscr{L}_2 (\beta, \lambda, \mu_1) = \beta^\prime \Sigma \beta - \lambda (\beta^\prime \beta - 1) - 2 \mu_1 \beta^\prime \Sigma \beta_1, $$
where we now have two lagrange multipliers, $\lambda$ and $\mu_1$. The partial derivative is
$$ \frac{\partial \mathscr{L}_2}{\partial \beta} = 2 \Sigma \beta - 2\lambda \beta - 2 \mu_1 \Sigma \beta_1, $$
and we set this to 0. We can left multiply this with $\beta_1^\prime$ to get
$$ 0 = 2\beta_1^\prime \Sigma \beta - 2 \lambda \beta_1^\prime \beta - 2 \mu_1 \beta_1^\prime \Sigma \beta_1 = -2 \mu_1 \lambda_1. $$
Therefore, $\mu_1 = 0$ and $\beta$ must satisfy $(\Sigma - \lambda I) \beta = 0$ again. Let $\lambda_{2}(\Sigma)$ be the second largest eigenvalue among $\lambda_{(1)}, \lambda_{(2)}, \dots, \lambda_{(p)}$ of $\Sigma$. Assuming that we don't have repeating eigenvalues, we then select the second largest eigenvalue and label the corresponding eigenvector as the second principal component $\beta_2$.

Similar procedure yields that the $k$th principal component is the eigenvector of $\Sigma$ with the $k$th largest eigenvalue.

If we let $B = [\beta_1 \dots \beta_p]$ and $ \Lambda = diag(\lambda_1, \lambda_2, \dots, \lambda_p) $, the equations $\Sigma \beta_r = \lambda_r \beta_r$ can be written in matrix form as $$ \Sigma B = B \Lambda, $$ and the equations $\beta_r^\prime \beta_r = 1$ and $\beta_r^\prime \beta_s = 0$ for $r\neq s$ can be written as $$ B^\prime B = I. $$
Thus, principal components are defined by an orthogonal transformation of the original random vector $X$.

## Some properties of PCA

1. Let the $p$ component random variable $X$ have $Var(X) = \Sigma$. Then there exists an orthogonal linear transformation $Z = B^\prime X$ such that $Var(Z) = \Lambda$ and $\Lambda = diag(\lambda_1, \lambda_2, \dots, \lambda_p)$ where $\lambda_1 \geq \lambda_2 \geq \dots \geq \lambda_p \geq 0$ are the eigenvalues of $\Sigma$.

Theorem

Given a symmetric matrix $\Sigma$, there exists an orthogonal matrix $C$ such that $C^\prime \Sigma C = D = diag(d_1, d_2, \dots, d_p)$. If $\Sigma$ is positive definite we have $d_i > 0$, if it is positive semidefinite then $d_i \geq 0$.

2. Orthogonal transformation of a random vector $X$ leaves the trace and determinant of the variance-covariance matrix invariant. That is $$ \sum_{i=1}^p Var(Z_i) tr(\Sigma_Z) = tr(Var(Z)) = tr(B \Sigma B^\prime) = tr(\Sigma B^\prime B) = tr(\Sigma) = \sum_{i=1}^p Var(X_i), $$
and $$ |\Sigma_z| = |B \Sigma B^\prime| = |\Sigma| \cdot |B B^\prime| = |\Sigma|. $$

Interpretation of this fact (theorem) may be that if we wish to reduce the dimensionality of our random variable, but wish to preserve as much variance as possible, we should take the first $r$ principal components and discard the rest (where $r < p$).

3. The transformation $Z = B^\prime X$ is a rotation of the coordinate axes. If we take a look at $$ x^\prime \Sigma x = x^\prime B^\prime \Lambda B x = z^\prime \Lambda z = \sum_{i=1}^p \frac{z_i^2}{1/\lambda_i} $$ and set it to a constant, i.e. $\sum_{i=1}^p \frac{z_i^2}{1/\lambda_i} = C$ we recognize an ellipse equation. The i-th principal component defines an ellipse axis with length $\sqrt{C/\lambda_i}$. 

![PCA ellipse](../assets/images/PCA_elipse.png)

## Inference

Consider now a sample of size $N$, $x_1, \dots, x_n$ iid drawn from some multivariate $p$-dimensional ($p<N$) distribution $X$ with $Var(X) = \Sigma$. This matrix $\Sigma$ has $p$ different characteristic roots, i.e. eigenvalues $\lambda_1, \dots, \lambda_p$ and eigenvectors $\beta_1, \dots, \beta_p$ as before.

In case that $X \sim \mathcal{N} (\mu, \Sigma)$, the set of maximum likelihood estimators of $\lambda_i$ and $\beta_i$ consists of the roots $l_1 > \dots > l_p$ of $$ |\hat{\Sigma} - lI| = 0 $$ and the corresponding mutually orthogonal vectors $b_1, \dots, b_p$ ($b_i^\prime b_i = 1$ and $b_i^\prime b_j = 0$) satisfy $$ (\hat{\Sigma} - l_i I) b_i = 0, $$
where $\hat{\Sigma}$ is the maximum likelihood estimate of $\Sigma$. Usually, we denote $\hat{\Sigma}$ by $S = \frac{1}{N} \sum_{i=1}^N (x_i - \overline{x})(x_i - \overline{x})^\prime$.

For general distributions, this matrix $S$ is not the MLE, but it remains a consistent estimator of $\Sigma$ under mild moment conditions. The PCA eigenvalues and eigenvectors of $S$ can still be analyzed asymptotically.

The most general result of the asymptotic eigenvalue distribution is $$ \sqrt{N} (l_i - \lambda_i) \xrightarrow{d} \mathcal{N} (0, \sigma_i^2) $$ where $$ \sigma_i^2 = Var((\beta_i^\prime X)^2). $$ The distribution of the eigenvectors is more complicated. If $X \sim \mathcal{N} (0, \Sigma)$, there are nicer closed form solutions: $$ \sqrt{N} (l_i - \lambda_i) \xrightarrow{d} \mathcal{N} (0, 2\lambda_i^2) $$ and $$ \sqrt{N} (b_i - \beta_i) \xrightarrow{d} \mathcal{N} (0, \sum_{k\neq i} \frac{\lambda_i \lambda_k}{(\lambda_i - \lambda_k)^2} \beta_k \beta_k^\prime). $$

To illustrate the asymptotic distributions, let's do a Monte Carlo simulation on a 2-dimensional probability distribution. We take a bivariate t-distribution, that is $$ X = \frac{Z}{\sqrt{W / \nu}}, \hspace{0.5cm} Z \sim \mathcal{N} (0, \Sigma), \hspace{0.5cm}, W \sim \chi_\nu^2, $$ where $\nu$ denotes the degrees of freedom. Set $$ \Sigma = \begin{pmatrix} 1 & 0.8 \\ 0.8 & 1 \end{pmatrix} $$ and $\nu = 5$.