---
title: "README"
author: "Zachary McCaw"
date: "2018-07-14"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Setting](#setting)
* [Example Data](#example-data)
* [Bivariate Outcome Regression](#bivariate-outcome-regression)
* [Multivariate Outcome Regression](#multivariate-outcome-regression)

# Setting

Suppose that a continuous response vector $y_{i}\in\mathbb{R}^{k}$ is observed for each of $n$ subjects. Partition the response vector as $y_{i} = (t_{i},s_{i})$. The first component $t_{i}\in\mathbb{R}^{1}$ is the (scalar) outcome of primary interest, or *target* outcome. The second component $s_{i}$ is the (potentially vectoral) outcome of secondary interest, or *surrogate* outcome. For example, suppose $t_{i}$ is the gold-standard measurement, and $s_{i}$ a collection of one or more surrogate measurements. 

Regression models are specified for the target $\mu_{T,i} = x_{i}'\beta$ and surrogate $\mu_{S,i} = \Xi_{i}\alpha$ means. Here $x_{i}$ is a $p\times 1$ vector of covariates for the target mean, and $\Xi_{i}$ is a $(k-1)\times q$ matrix of covariates for the surrogate mean. Collect $x_{i}$ and $\Xi_{i}$ into the subject-specific design matrix $Z_{i}$. The overall regression model is:
$$
\mu_{i} \equiv \binom{\mu_{T,i}}{\mu_{S,i}} = \left(
\begin{array}{c c}
x_{i}' & 0 \\
0 & \Xi_{i}
\end{array}
\right)
\binom{\beta}{\alpha} \equiv Z_{i}\gamma
$$

Conditional on covariates, the response vector $y_{i}\in\mathbb{R}^{k}$ is taken as multivariate normal:

$$
y_{i}|Z_{i} \sim N \binom{\mu_{T,i}}{\mu_{S,i}}, \left(
\begin{array}{c c}
\Sigma_{TT} & \Sigma_{TS} \\
\Sigma_{ST} & \Sigma_{SS}
\end{array}
\right)
$$

Finally, partition the regression coefficient for the target mean as $\beta = (\beta_{1},\beta_{2})$. This package provides score tests for assessing $H_{0}:\beta_{1}=\beta_{10}$. The advantage of jointly modeling the target and surrogate outcomes is improved power for performing inference on the target mean model.

# Example Data
## Bivariate Outcome
The list `D.bnr` contains simulated data for $n=1000$ subjects with bivariate normal outcomes. The target `yt` and surrogate `ys` outcomes have unit variances and correlation $\rho=0.5$. Each outcome depends on a model matrix containing an intercept and three independent, standard normal covariates. The model matrix for the target outcome is `Zt`, and the model matrix for the surrogate outcome is `Zs`. `Beta` contains the regression coefficients for the target mean model, and `Alpha` contains the regression coefficients for the surrogate mean model. 


```r
library(MNR);
# Extract components
yt = D.bnr$yt;
ys = D.bnr$ys;
Zt = D.bnr$Zt;
Zs = D.bnr$Zs;
cat("Target outcome:","\n");
head(yt);
cat("\n");
cat("Surrogate outcome:","\n");
head(ys);
cat("\n");
cat("Target model matrix:","\n");
head(Zt);
cat("\n");
cat("Surrogate model matrix:","\n");
head(Zs);
cat("\n");
cat("Target regression parameters:","\n");
show(D.bnr$Beta);
cat("\n");
cat("Surrogate regression parameters:","\n");
show(D.bnr$Alpha);
```

```
## Target outcome: 
##          1          2          3          4          5          6 
## -2.1775941 -0.6433499 -1.2380953 -1.1208658  0.6387681 -0.9404964 
## 
## Surrogate outcome: 
##           1           2           3           4           5           6 
## -0.00438639  1.15376520  1.25420024  0.14000159  0.28943876  0.88898029 
## 
## Target model matrix: 
##   int.t          x1         x2          x3
## 1     1  0.21475678  0.4998932 -0.08250593
## 2     1  1.59638109 -0.8656795  1.66117669
## 3     1  0.09539529  0.4924589 -1.97274622
## 4     1 -1.25801656 -1.3526159 -2.14845755
## 5     1  0.06733360 -0.3943881 -1.14182522
## 6     1 -0.78741717  0.7711460  2.98557356
## 
## Surrogate model matrix: 
##   int.s         xi1         xi2         xi3
## 1     1 -0.60635079 -0.01501802  0.21116190
## 2     1 -1.25536475  0.48508411  0.63642953
## 3     1  1.08003701 -0.37448487  0.09266153
## 4     1  0.08325608 -1.18135648  0.01246761
## 5     1 -0.09917261  1.19255878 -0.90633463
## 6     1  2.49662596  0.57356246 -1.08234428
## 
## Target regression parameters: 
##   b0   b1   b2   b3 
## -1.0  0.1  0.1  0.0 
## 
## Surrogate regression parameters: 
##   a0   a1   a2   a3 
##  1.0 -0.1 -0.1 -0.1
```
## Trivariate Outcome
The list `D.mnr` contains a simulated data set for $n=1000$ subjects with trivariate normal outcomes. The target `yt` and surrogate `Ys` outcomes have unit variances, and an exchangeable correlation structure, with $\rho=0.5$. The target and first surrogate outcomes each depend on three independent, standard normal covariates. The second surrogate outcome depends on four such covariates. The model matrix for the target outcome is `Zt`. The list `Ls` contains two model matrices. The first, `Xi1` is an $1000\times 4$ model matrix for the first surrogate outcome, consistening of an intercept and the three covariates. The second, `Xi2` is an $1000\times 5$ model matrix for the second surrogate outcome, consisting of an intercept and the four covariates. The overall surrogate model matrix `Zs` is generated within the testing function. By default, each surrogate outcome is permitted its own regression parameters. To fit a *parallel coefficient model*, in which the surrogates share a common set of regression coefficients, specify `parallel=T` during hypothesis testing. Fitting a parallel coefficient model requires that each design matrix has the same columns. `Beta` contains the regression coefficients for the target mean model, and `Alpha` contains the regression coefficients for the surrogate mean models.


```r
library(MNR);
# Extract components
yt = D.mnr$yt;
Ys = D.mnr$Ys;
Zt = D.mnr$Zt;
Ls = D.mnr$Ls;
cat("Target outcome:","\n");
head(yt);
cat("\n");
cat("Surrogate outcomes:","\n");
head(Ys);
cat("\n");
cat("Target model matrix:","\n");
head(Zt);
cat("\n");
cat("Surrogate model matrices:","\n");
head(Ls$Xi1);
cat("\n");
head(Ls$Xi2);
cat("Target regression parameters:","\n");
show(D.mnr$Beta);
cat("\n");
cat("Surrogate regression parameters","\n");
show(D.mnr$Alpha);
```

```
## Target outcome: 
##          1          2          3          4          5          6 
## -1.1355092 -1.1555698 -3.1820548 -2.3106878 -0.5271279 -1.5898073 
## 
## Surrogate outcomes: 
##        y.s1      y.s2
## 1 1.0859936 1.8317890
## 2 1.1687834 2.1326281
## 3 0.3840986 1.0835291
## 4 0.5095097 0.1996144
## 5 1.0211076 2.2989638
## 6 1.9449832 2.0919847
## 
## Target model matrix: 
##   int.t         x1           x2          x3
## 1     1 -0.9157211 -0.967320838 -1.06790563
## 2     1 -1.6354771 -0.008087443  0.09711701
## 3     1 -0.2941834 -0.147777850 -0.87567080
## 4     1 -0.1428488  1.778802601  0.63141124
## 5     1  1.1652948  0.468497648  0.71073700
## 6     1  0.6457990 -0.227475132  2.03778239
## 
## Surrogate model matrices: 
##      int.s1       xi11       xi12          xi13
## [1,]      1 -1.4395537 -0.3740681  1.5130474139
## [2,]      1  0.1674263 -0.2899386  2.6596479915
## [3,]      1  0.8630031 -0.5530195  0.9662611896
## [4,]      1  1.5602511  0.3050545  0.1524413326
## [5,]      1  0.7642962 -1.1030871  0.2749556706
## [6,]      1  1.1860374  0.3805394 -0.0007527277
## 
##      int.s2        xi21       xi22        xi23       xi24
## [1,]      1  1.34363436 -1.1195984 -0.00130081 -0.6959728
## [2,]      1 -0.42858546 -1.0723851 -0.44431810  0.4985890
## [3,]      1 -0.87064367 -2.0306128 -0.08852820 -0.2751421
## [4,]      1 -0.07328819 -1.7535199 -1.11290174 -0.8265508
## [5,]      1 -2.00706806 -0.8361248 -0.03223253  0.2806040
## [6,]      1  0.71148498  1.1797593 -2.23657957 -0.4730092
## Target regression parameters: 
##   b0   b1   b2   b3 
## -1.0  0.1  0.1  0.0 
## 
## Surrogate regression parameters 
##  a10  a11  a12  a13  a20  a21  a22  a23  a24 
##  1.0 -0.1 -0.1 -0.1  2.0  0.1  0.1  0.1  0.1
```

## Formatting Assumptions
The target outcome `yt` is a numeric vector of length $n$. In the case of a single surrogate, `ys` is a numeric vector of length $n$. In the case of multiple surrogates, `Ys` is an $n\times(k-1)$ numeric matrix. The target model `Zt` is an $n\times p$ numeric matrix. In the case of a single surrogate, `Zs` is an $n\times q$ numeric matrix. In the case of multiple surrogates, $Ls$ is a list of numeric matrices, each with $n$ rows. In the case of a parallel coefficient model, each matrix in $Ls$ must additionally have the same number of columns $q$. All factors and interactions must have been expanded in advance. None of the inputs should contain missing data. 

# Bivariate Outcome Regression
The case of a single surrogate is implemented separately from the case of multiple surrogates because the fitting procedure is faster in the former. The following demonstrates various score tests possible using the example data. The first is an overall test of $H_{0}:\beta_{1}=\beta_{2}=\beta_{3}=0$, which is false. The second assesses $H_{0}:\beta_{1}=\beta_{2}=0$, which is again false, treating $\beta_{3}$ as a nuisance. The third considers $H_{0}:\beta_{1}=0.1$, which is in fact true, treating $\beta_{2}$ and $\beta_{3}$ as nuisances. The final considers $H_{0}:\beta_{3}=0$, which is again true, treating $\beta_{1}$ and $\beta_{2}$ as nuisances. All models include an intercept.


```r
# Extract components
yt = D.bnr$yt;
ys = D.bnr$ys;
Zt = D.bnr$Zt;
Zs = D.bnr$Zs;
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(Score.bnr(yt,ys,Zt,Zs,L=c(F,T,T,T)),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(Score.bnr(yt,ys,Zt,Zs,L=c(F,T,T,F)),digits=2);
cat("\n","Joint score test of b1 = 0.1, treating b2 and b3 as nuisances","\n");
signif(Score.bnr(yt,ys,Zt,Zs,L=c(F,T,F,F),b10=c(0.1)),digits=2);
cat("\n","Individual score test of b3 = 0, treating b1 and b2 as nuisances","\n");
signif(Score.bnr(yt,ys,Zt,Zs,L=c(F,F,F,T)),digits=2);
```

```
## Joint score test of b1 = b2 = b3 = 0 
##   Score      df       p 
## 1.7e+01 3.0e+00 5.9e-04 
## 
##  Joint score test of b1 = b2 = 0, treating b3 as a nuisance 
##   Score      df       p 
## 1.7e+01 2.0e+00 1.8e-04 
## 
##  Joint score test of b1 = 0.1, treating b2 and b3 as nuisances 
## Score    df     p 
##   1.6   1.0   0.2 
## 
##  Individual score test of b3 = 0, treating b1 and b2 as nuisances 
## Score    df     p 
##  0.19  1.00  0.67
```

# Trivariate Outcome Regression
Hypothesis testing for two or more surrogate outcomes is analogous to the bivariate case. The same hypothesis tests discussed above are repeated for the trivariate outcome model. 


```r
# Extract components
yt = D.mnr$yt;
Ys = D.mnr$Ys;
Zt = D.mnr$Zt;
Ls = D.mnr$Ls;
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(Score.mnr(yt,Ys,Zt,Ls,L=c(F,T,T,T)),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(Score.mnr(yt,Ys,Zt,Ls,L=c(F,T,T,F)),digits=2);
cat("\n","Joint score test of b1 = 0.1, treating b2 and b3 as nuisances","\n");
signif(Score.mnr(yt,Ys,Zt,Ls,L=c(F,T,F,F),b10=c(0.1)),digits=2);
cat("\n","Individual score test of b3 = 0, treating b1 and b2 as nuisances","\n");
signif(Score.mnr(yt,Ys,Zt,Ls,L=c(F,F,F,T)),digits=2);
```

```
## Joint score test of b1 = b2 = b3 = 0 
##   Score      df       p 
## 5.1e+01 3.0e+00 5.4e-11 
## 
##  Joint score test of b1 = b2 = 0, treating b3 as a nuisance 
##   Score      df       p 
## 5.1e+01 2.0e+00 1.0e-11 
## 
##  Joint score test of b1 = 0.1, treating b2 and b3 as nuisances 
## Score    df     p 
##  0.44  1.00  0.51 
## 
##  Individual score test of b3 = 0, treating b1 and b2 as nuisances 
## Score    df     p 
##  0.38  1.00  0.54
```
