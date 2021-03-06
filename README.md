---
title: "README"
author: "Zachary McCaw"
date: "2018-09-11"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Introduction](#introduction)
* [Example Data](#example-data)
* [Estimation](#estimation)

# Introduction

## Overview

Suppose that a multivariate normal outcome vector is observed for each subject. This package considers the general case where separate regression models are specified for each element of the outcome vector, and the outcome covariance matrix is left unstructured. Utilities are provided for model estimation and for inference on the regression parameters. 

## Model

Suppose an outcome vector $y_{i}\in\mathbb{R}^{k}$ is observed for each subject. Let $y_{ij}$ denote the $j$th element of $y_{i}$. A regression model for the subject-specific mean of the $j$th element is given by:

$$
\mu_{ij} = x_{ij}'\beta_{j}
$$

Conditional on covariates, the outcome follows a multivariate normal distribution with unstructured covariance:

$$
y_{i}\big|(x_{i1},\cdots,x_{ik}) \sim N \left(
\begin{array}{c}
x_{i1}'\beta_{1} \\
\vdots \\
x_{ik}'\beta_{k}
\end{array}
\right), \left(
\begin{array}{c c c c}
\Sigma_{11} & \Sigma_{12} & \cdots & \Sigma_{1k} \\
\Sigma_{21} & \Sigma_{22} & & \vdots \\
\vdots & & \ddots & \\
\Sigma_{k1} & \cdots & & \Sigma_{kk}
\end{array}
\right)
$$

Suppose $\beta_{j}$ is the covariate of interest. Partition $\beta_{j}=(\beta_{1,j},\beta_{2,j})$. Consider the hypothesis $\beta_{1,j} = \beta_{1,j}^{\dagger}$. Let $\eta_{1,j}$ denote the collected nuisance regression parameters. The score test of $H_{0}:\beta_{1,j} = \beta_{1,j}^{\dagger}$ takes the form:

$$
T_{S} = U_{\beta_{1,j}}(\beta_{1,j}^{\dagger},\eta_{1,j})'I_{\beta_{1,j}\beta_{1,j}\big|\eta_{1,j}}^{-1}
U_{\beta_{1,j}}(\beta_{1,j}^{\dagger},\eta_{1,j})
$$

# Example Data

Below, data is simulated for $n=10^{3}$ subjects. The outcome is trivariate normal $y_{i}\in\mathbb{R}^{3}$. The covariance structure is exchangeable with diagonal one and off diagonal $\rho=0.5$. Separate regression models are specified for each component of $y_{i}$. The mean of the first, second, and third components of the outcome depends on a design matrix with an intercept and two, three, or four independent, standard normal covariates, respectively. 


```r
library(MNR);
set.seed(100);
# Observations
n = 1e3;
## Design matrices
X1 = cbind(1,matrix(rnorm(2*n),nrow=n));
colnames(X1) = c("int",paste0("x0",seq(1:2)));
X2 = cbind(1,matrix(rnorm(3*n),nrow=n));
colnames(X2) = c("int",paste0("x1",seq(1:3)));
X3 = cbind(1,matrix(rnorm(4*n),nrow=n));
colnames(X3) = c("int",paste0("x2",seq(1:4)));
X = list(X1,X2,X3);
# Target Parameter
b1 = c(-1,0.1,-0.1);
b2 = c(1,-0.1,0.1,0);
b3 = c(0,0.1,-0.1,0.1,-0.1);
b = list(b1,b2,b3);
# Exchangeable covariance structure
S = array(0.5,dim=c(3,3)) + 0.5*diag(3);
# Generate data
Y = rMNR(X=X,b=b,S=S);
```

## Data Formatting

The outcome `Y` is expected as a numeric matrix, with observations as rows. Covariates are supplied as a *list* of numeric matrices, one for each column of `Y`. 

# Estimation

Below, parameters are estimated for the trivariate normal data generated above. The fitted model is of class `mnr`. The show method provides a table of estimated regression coefficients, along with standard errors, and Wald-type confidence intervals and p-values. The fitted model `M` contains the following components:

* Regression coefficients, extracted using `coef(M)` or `M@Coefficients`.
* Outcome covariance matrix, extracted using `vcov(M,type="Outcome")` or `M@Covariance`.
* Regression coefficient information matrix, extracted using `vcov(M,type="Information")` or `M@Information`.
* Outcome residual matrix, extracted using `resid(M)` or `M@Residuals`. 


```r
# Model fitting
M = fit.mnr(Y=Y,X=X,eps=1e-8);
cat("Show method for fitted model:\n");
show(M);
cat("\n");
cat("Comparison of estimated coefficients with the truth:\n");
Coeff = coef(M);
show(data.frame(Coeff[,c("Outcome","Coeff")],"Est"=round(Coeff$Point,digits=3),"Truth"=unlist(b)));
```

```
## Objective increment:  6.41 
## Objective increment:  0.00253 
## Objective increment:  1.16e-06 
## Objective increment:  3.79e-10 
## 3 update(s) performed before tolerance limit.
## 
## Show method for fitted model:
##    Outcome Coeff    Point     SE       L       U         p
## 1       y1   int -1.00000 0.0312 -1.0600 -0.9420 3.55e-227
## 2       y1   x01  0.14500 0.0248  0.0962  0.1930  5.33e-09
## 3       y1   x02 -0.09390 0.0261 -0.1450 -0.0427  3.21e-04
## 4       y2   int  0.99900 0.0308  0.9390  1.0600 1.08e-230
## 5       y2   x11 -0.02940 0.0244 -0.0771  0.0184  2.28e-01
## 6       y2   x12  0.12000 0.0266  0.0680  0.1720  6.27e-06
## 7       y2   x13  0.01060 0.0263 -0.0408  0.0621  6.85e-01
## 8       y3   int  0.00247 0.0303 -0.0569  0.0618  9.35e-01
## 9       y3   x21  0.11400 0.0262  0.0622  0.1650  1.48e-05
## 10      y3   x22 -0.11000 0.0259 -0.1610 -0.0593  2.11e-05
## 11      y3   x23  0.07000 0.0261  0.0188  0.1210  7.37e-03
## 12      y3   x24 -0.07820 0.0258 -0.1290 -0.0276  2.44e-03
## 
## Comparison of estimated coefficients with the truth:
##    Outcome Coeff    Est Truth
## 1       y1   int -1.003  -1.0
## 2       y1   x01  0.145   0.1
## 3       y1   x02 -0.094  -0.1
## 4       y2   int  0.999   1.0
## 5       y2   x11 -0.029  -0.1
## 6       y2   x12  0.120   0.1
## 7       y2   x13  0.011   0.0
## 8       y3   int  0.002   0.0
## 9       y3   x21  0.114   0.1
## 10      y3   x22 -0.110  -0.1
## 11      y3   x23  0.070   0.1
## 12      y3   x24 -0.078  -0.1
```

# Inference

Below, hypotheses are tested on the model fit to the trivariate normal data generated above. The column of `Y` which is of interest is specified using `j`. The hypothesis test is specified using a logical vector `L`. Elements of `L` that are constrained under the null are set to `TRUE`, those which are estimated under the null are set to `FALSE`. At least one element of `L` must be `TRUE` (a test must be specified), and at least one element of `L` must be `FALSE` (the null model must be estimable). 

* The first assesses $H_{0}:\beta_{13}=0$, which is false. 
* The second test assesses $H_{0}:\beta_{24}=0$, which is true. 
* The third test assesses $H_{0}:\beta_{32} = \cdots = \beta_{35} = 0$, which is false. 
* The fouth test assesses $H_{0}:\beta_{32} = \beta_{34} = 0.1$, which is true. 


```r
cat("Test b13 = 0, which is FALSE:\n");
Score.mnr(Y=Y,j=1,X=X,L=c(F,F,T));
cat("\n");
cat("Test b24 = 0, which is TRUE:\n");
Score.mnr(Y=Y,j=2,X=X,L=c(F,F,F,T));
cat("\n");
cat("Test b32 = ... = b35 = 0, which is FALSE:\n");
Score.mnr(Y=Y,j=3,X=X,L=c(F,T,T,T,T));
cat("\n");
cat("Test b32 = b34 = 0.1, which is TRUE:\n");
Score.mnr(Y=Y,j=3,X=X,b10=c(0.1,0.1),L=c(F,T,F,T,F));
```

```
## Test b13 = 0, which is FALSE:
##        Score           df            p 
## 12.572560188  1.000000000  0.000391452 
## 
## Test b24 = 0, which is TRUE:
##     Score        df         p 
## 0.1637477 1.0000000 0.6857293 
## 
## Test b32 = ... = b35 = 0, which is FALSE:
##        Score           df            p 
## 5.088565e+01 4.000000e+00 2.358459e-10 
## 
## Test b32 = b34 = 0.1, which is TRUE:
##     Score        df         p 
## 1.5734197 2.0000000 0.4553405
```

## Repeated Score Test

Consider testing for association between a specified column of `Y`, and each column of a matrix `G`. The function `rScore.mnr` accelerates association testing by recycling the same null model for each hypothesis test. This function is motivated by genetic association testing, where `G` represents genotype at different loci across the genome. The *genotype* matrix `G` may contain missing values, although the outcome `Y` and design matrices `X` still should not. 


```r
# Genotype matrix
G = replicate(2000,rbinom(n=1000,size=2,prob=0.25));
storage.mode(G) = "numeric";
# Introduce missingness
G[sort(sample(length(G),size=0.01*length(G)))] = NA;
# Repeated Score test
doMC::registerDoMC(cores=4);
R = rScore.mnr(Y=Y,G=G,X=X,report=T,parallel=T);
# Estimated size
mean(R<=0.05);
```


