# README
Zachary McCaw  
`r Sys.Date()`  

# Package Vignette




# Contents

* [Setting](#setting)
* [Example Data](#example-data)
* [Bivariate Outcome Regression](#bivariate-outcome-regression)
* [Multivariate Outcome Regression](#multivariate-outcome-regression)

# Setting

Suppose that a continuous response vector $y_{i}\in\mathbb{R}^{k}$ is observed for each of $n$ subjects. Partition the response vector as $y_{i} = (t_{i},s_{i})$, where $t_{i}\in\mathbb{R}^{1}$ is the outcome of primary interest, or *target* outcome, and $s_{i}\in\mathbb{R}^{k-1}$ is an outcome of secondary interest. For example, $t_{i}$ may denote a gold-standard measurement, and $s_{i}$ a collection of *surrogate* measurements. In addition, for each subject, a matrix $Z_{i}$ of covariates is observed. Conditional on covariates, the response vector is taken to follow a multivariate normal distribution:

$$
y_{i}|Z_{i} \sim N \binom{\mu_{T,i}}{\mu_{S,i}}, \left(
\begin{array}{c c}
\Sigma_{TT} & \Sigma_{TS} \\
\Sigma_{ST} & \Sigma_{SS}
\end{array}
\right)
$$

Let $x_{T,i}$ denote a $p \times 1$ vector of covariates relevant to the target outcome $t_{i}$, and $X_{S,i}$ a $(k-1) \times q$ matrix of covariates relevant to the secondary outcome $s_{i}$. Specify regression models for the target and secondary means as:

$$
\binom{\mu_{T,i}}{\mu_{S,i}} = \left(
\begin{array}{c c}
x_{T,i}' & 0 \\
0 & X_{S,i}
\end{array}
\right)
\binom{\beta}{\alpha}
$$

Finally, partition the regression coefficient for the target mean as $\beta = (\beta_{1},\beta_{2})$. Suppose interest lies in testing the null hypothesis that a subset of the target regression coefficients agree with a pre-specified value, i.e. $H_{0}:\beta_{1} = \beta_{10}$. Score tests for evaluating this hypothesis are provided in the case of a single and of multiple secondary outcomes. 

# Example Data
## Bivariate Outcome
The list `D.bnr` contains a simulated data set for $n=1000$ subjects with bivariate normal outcomes. The target and surrogate outcomes have unit variances and correlation $\rho=0.5$. Each outcome depends on a model matrix containing an intercept and three independent, standard normal covariates. `Beta` contains the regression coefficients linking the target outcome `y.t` to the target model matrix `Z.t`. Likewise, `Alpha` contains the regression coefficients linking the surrogate outcome `y.s` to the surrogate model matrix `Z.s`.


```r
library(MNR);
# Target outcome
y.b.t = D.bnr$y.t;
# Surrogate outcome
y.b.s = D.bnr$y.s;
# Target model matrix
Z.b.t = D.bnr$Z.t;
# Surrogate model matrix
Z.b.s = D.bnr$Z.s;
# Beta
print(D.bnr$Beta);
# Alpha
print(D.bnr$Alpha);
```

```
##   b0   b1   b2   b3 
## -1.0  0.1  0.1  0.0 
##   a0   a1   a2   a3 
##  1.0 -0.1 -0.1 -0.1
```
## Trivariate Outcome
The list `D.mnr` contains a simulated data set for $n=1000$ subjects with trivariate normal outcomes. The target and surrogate outcomes have unit variances, and an exchangeable correlation structure, with $\rho=0.5$. The target and first surrogate outcomes depend on three independent standard normal covariates, while the second surrogate outcome depends on four such covariates. `Beta` contains the regression coefficients linking the target outcome `y.t` to the target model matrix `Z.t`. Likewise, `Alpha` contains the regression coefficients linking the surrogate outcomes `y.s` to the surrogate model matrices in `L.s`. In the case of multiple surrogates, `L.s` is a *list* of design matrices, with one matrix per surrogate outcome. By default, each surrogate mean is afford its own regression parameters. To fit a **parallel coefficient model**, in which the surrogate means share a common set of regression coefficients, set `parallel=T`. Fitting a parallel coefficient model requires that each design matrix has the same columns. 


```r
library(MNR);
# Target outcome
y.m.t = D.mnr$y.t;
# Surrogate outcome
y.m.s = D.mnr$y.s;
# Target model matrix
Z.m.t = D.mnr$Z.t;
# Surrogate model matrices
Z.m.s = D.mnr$L.s;
# Beta
print(D.mnr$Beta);
# Alpha
print(D.mnr$Alpha);
```

```
##   b0   b1   b2   b3 
## -1.0  0.1  0.1  0.0 
##  a10  a11  a12  a13  a20  a21  a22  a23  a24 
##  1.0 -0.1 -0.1 -0.1  2.0  0.1  0.1  0.1  0.1
```

# Bivariate Outcome Regression
The case of a single secondary outcome is implemented separately from the case of multiple secondary outcomes because the fitting procedure is faster in the former. The following demonstrates various score tests possible using the example data. The first is an overall test of $H_{0}:\beta_{1}=\beta_{2}=\beta_{3}=0$, which is false. The second assesses $H_{0}:\beta_{1}=\beta_{2}=0$, which is again false, treating $\beta_{3}$ as a nuisance. The third considers $H_{0}:\beta_{1}=0.1$, which is in fact true, while treating $\beta_{2}$ and $\beta_{3}$ as nuisances. The final considers $H_{0}:\beta_{3}=0$, which is again true, treating $\beta_{1}$ and $\beta_{2}$ as nuisances. All models include an intercept. 


```r
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(Score.bnr(y.b.t,y.b.s,Z.b.t,Z.b.s,L=c(F,T,T,T)),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(Score.bnr(y.b.t,y.b.s,Z.b.t,Z.b.s,L=c(F,T,T,F)),digits=2);
cat("\n","Joint score test of b1 = 0.1, treating b2 and b3 as nuisances","\n");
signif(Score.bnr(y.b.t,y.b.s,Z.b.t,Z.b.s,L=c(F,T,F,F),b0=c(0.1)),digits=2);
cat("\n","Individual score test of b3 = 0, treating b1 and b2 as nuisances","\n");
signif(Score.bnr(y.b.t,y.b.s,Z.b.t,Z.b.s,L=c(F,F,F,T)),digits=2);
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
Hypothesis testing for two or more surrogate outcomes is analogous to the bivariate case. The same hypothesis tests considered above are repeated for the trivariate outcome model. 


```r
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(Score.mnr(y.m.t,y.m.s,Z.m.t,Z.m.s,L=c(F,T,T,T)),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(Score.mnr(y.m.t,y.m.s,Z.m.t,Z.m.s,L=c(F,T,T,F)),digits=2);
cat("\n","Joint score test of b1 = 0.1, treating b2 and b3 as nuisances","\n");
signif(Score.mnr(y.m.t,y.m.s,Z.m.t,Z.m.s,L=c(F,T,F,F),b0=c(0.1)),digits=2);
cat("\n","Individual score test of b3 = 0, treating b1 and b2 as nuisances","\n");
signif(Score.mnr(y.m.t,y.m.s,Z.m.t,Z.m.s,L=c(F,F,F,T)),digits=2);
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
