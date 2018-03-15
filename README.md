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
The list `D.bvr` contains a simulated data set for $n=1000$ subjects with bivariate normal outcomes. The target and surrogate outcomes have unit variances and correlation $\rho=0.5$. Each outcome depends on a design matrix containing three independent standard normal covariates. `Beta` contains the regression coefficients used to generate the subject-specific target mean $\mu_{T,i}$ from the target design matrix `D.t`. Likewise, `Alpha` contains the regression coefficients used to generate the subject-specific surrogate mean $\mu_{S,i}$ from the surrogate design matrix `D.s`. 


```r
library(MNR);
# Target mean
y.t = D.bvr$y.t;
# Surrogate mean
y.s = D.bvr$y.s;
# Target design
D.t = D.bvr$D.t;
# Surrogate design
D.s = D.bvr$D.s;
# Beta
print(D.bvr$Beta);
# Alpha
print(D.bvr$Alpha);
```

```
##   b0   b1   b2   b3 
## -1.0  0.1  0.1  0.0 
##   a0   a1   a2   a3 
##  1.0 -0.1 -0.1 -0.1
```

# Bivariate Outcome Regression
The case of a single secondary outcome is implemented separately from the case of multiple secondary outcomes because the fitting procedure is simplified, and thereby faster, in the former. The following demonstrates various score tests possible using the example data. The first is an overall test of $H_{0}:\beta_{1}=\beta_{2}=\beta_{3}=0$, which is false. The second assesses whether $H_{0}:\beta_{1}=\beta_{2}=0$ holds, which is again false, treating $\beta_{3}$ as a nuisance. The third considers $H_{0}:\beta_{1}=0.1$, which is in fact true, while treating $\beta_{2}$ and $\beta_{3}$ as nuisances. The final considers $H_{0}:\beta_{3}=0$, which is again true, treating $\beta_{1}$ and $\beta_{2}$ as nuisances


```r
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(Score.bnr(y.t,y.s,D.t,D.s,L=c(T,T,T)),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(Score.bnr(y.t,y.s,D.t,D.s,L=c(T,T,F)),digits=2);
cat("\n","Joint score test of b1 = 0.1, treating b2 and b3 as nuisances","\n");
signif(Score.bnr(y.t,y.s,D.t,D.s,L=c(T,F,F),b0=c(0.1)),digits=2);
cat("\n","Individual score test of b3 = 0, treating b1 and b2 as nuisances","\n");
signif(Score.bnr(y.t,y.s,D.t,D.s,L=c(F,F,T)),digits=2);
```

```
## Joint score test of b1 = b2 = b3 = 0 
##   Score      df       p 
## 2.7e+01 3.0e+00 6.6e-06 
## 
##  Joint score test of b1 = b2 = 0, treating b3 as a nuisance 
##   Score      df       p 
## 2.7e+01 2.0e+00 1.6e-06 
## 
##  Joint score test of b1 = 0.1, treating b2 and b3 as nuisances 
## Score    df     p 
##  0.58  1.00  0.45 
## 
##  Individual score test of b3 = 0, treating b1 and b2 as nuisances 
##  Score     df      p 
## 0.0031 1.0000 0.9600
```
