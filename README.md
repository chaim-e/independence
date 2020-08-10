# independence (R package)

## Fast Rank-Based Independence Tests

Install:
```
library(devtools)
install_github("chaim-e/independence")
```
Examples:
```
tau.star.test(rnorm(1000),rnorm(1000))
hoeffding.D.test(rnorm(1000),rnorm(1000))
hoeffding.refined.test(rnorm(1000),rnorm(1000))
```
