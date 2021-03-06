# independence (R package)

## Fast Rank-Based Independence Tests
Performs three ranking-based nonparametric tests for the independence of two continuous variables:
(1) the classical Hoeffding's D test; 
(2) a refined variant of it, named R;
(3) the Bergsma-Dassios T* sign covariance.
The first test is consistent assuming an absolutely continuous bivariate distribution, i.e., the population coefficient D=0 iff the variables are independent. The latter two are consistent under no restriction on the distribution.
All three statistics are computed in time O(n log n) given n iid paired samples. The computation of R and T* uses a new algorithm, following work of Even-Zohar and Leng (2019), see https://arxiv.org/abs/2010.09712 and references therein.
## Install
From CRAN
```
install.packages("independence")
library(independence)
```
From GitHub
```
library(devtools)
install_github("chaim-e/independence")
library(independence)
```
## Examples
```
tau.star.test(rnorm(1000000),rnorm(1000000))
hoeffding.D.test(rnorm(1000000),rnorm(1000000))
hoeffding.refined.test(rnorm(1000000),rnorm(1000000))
```
### Why you should use the refined Hoeffding's D or tau*
Example (I)
```
set.seed(12345)
f <- function(a,b) ifelse(a>b, pmin(b,a/2), pmax(b,(a+1)/2))
x <- runif(300)
y <- f(x, runif(300))
hoeffding.D.test(x,y)$p.value
# [1] 0.4589397
hoeffding.refined.test(x,y)$p.value
# [1] 2.138784e-10
tau.star.test(x,y)$p.value
# [1] 1.053099e-08
```
Example (II)
```
set.seed(13579)
f <- function(a,b) ifelse(a>0.5, pmin(b,1-0.5/a), pmax(b,0.5/(1-a)))
x <- runif(50)
y <- f(x, runif(50))
hoeffding.D.test(x,y)$p.value
# [1] 0.4157275
hoeffding.refined.test(x,y)$p.value
# [1] 8.04458e-08
tau.star.test(x,y)$p.value
# [1] 1.009128e-05
```
Example (III)
```
set.seed(54321)
f <- function(a,b) (a+b-b%%(2^ceiling(log(pmin(a,1-a),2))))%%1
x <- runif(1000)
y <- f(x, runif(1000))
hoeffding.D.test(x,y)$p.value
# [1] 0.5391983
hoeffding.refined.test(x,y)$p.value
# [1] 6.539453e-09
tau.star.test(x,y)$p.value
# [1] 8.865218e-07
```
