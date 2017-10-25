





[![version](http://www.r-pkg.org/badges/version/personalized)](https://cran.r-project.org/package=personalized)
[![Build Status](https://travis-ci.org/jaredhuling/personalized.svg?branch=master)](https://travis-ci.org/jaredhuling/personalized)
[![codecov](https://codecov.io/gh/jaredhuling/personalized/branch/master/graph/badge.svg)](https://codecov.io/gh/jaredhuling/personalized)
[![](http://cranlogs.r-pkg.org/badges/personalized)](http://cran.rstudio.com/web/packages/personalized/index.html)


# Overview of 'personalized'

The 'personalized' package is designed for the analysis of data where the effect of a treatment or intervention may vary for different patients. It can be used for either data from randomized controlled trials or observational studies and is not limited specifically to the analysis of medical data.

The personalized package provides estimation methods for subgroup identification under the framework of [Chen et al (2017)](http://onlinelibrary.wiley.com/doi/10.1111/biom.12676/abstract). It also provides routines for valid estimation of the subgroup-specific treatment effects.


<img src="vignettes/usage_overview-1.png" width="100%" />

## Documentation

[Documentation](http://jaredhuling.org/personalized/)

- [Tutorial / Vignette](http://jaredhuling.org/personalized/articles/usage_of_the_personalized_package.html) - tutorial of statistical methodology and usage of the package

- [Function help files](http://jaredhuling.org/personalized/reference/index.html)


# Installing the 'personalized' package


Install from CRAN using:

```r
install.packages("personalized")
```

or install the development version using the **devtools** package:

```r
devtools::install_github("jaredhuling/personalized")
```

or by cloning and building using `R CMD INSTALL`

# Quick Usage Overview

Load the package:

```r
library(personalized)
```



### Create a propensity score model 
(it should be a function which inputs covariates and treatments and returns propensity score):

```r
prop.func <- function(x, trt)
{
    # fit propensity score model
    propens.model <- cv.glmnet(y = trt,
                               x = x, family = "binomial")
    pi.x <- predict(propens.model, s = "lambda.min",
                    newx = x, type = "response")[,1]
    pi.x
}
```

### Fit a model to estimate subgroup:

```r
subgrp.model <- fit.subgroup(x = x, y = y,
                             trt = trt,
                             propensity.func = prop.func,
                             loss   = "sq_loss_lasso",
                             nfolds = 5)              # option for cv.glmnet
```

### Display estimated subgroups and variables selected which determine the subgroups:

```r
summary(subgrp.model)
```

```
## family:  gaussian 
## loss:    sq_loss_lasso 
## method:  weighting 
## propensity 
## function: propensity.func 
## 
## Average Outcomes:
##                 Recommended Ctrl    Recommended Trt
## Received Ctrl  -5.6613 (n = 117) -23.6003 (n = 114)
## Received Trt  -20.9028 (n = 132)  -7.3036 (n = 137)
## 
## Ctrl effect among recommended Ctrl   Trt effect among recommended Trt 
##                  15.2414 (n = 249)                  16.2968 (n = 251) 
## 
## Benefit score quantiles: 
##        0%       25%       50%       75%      100% 
## -14.15602  -3.58120   0.04648   3.51676  14.78106 
## 
## 9 out of 50 variables selected in total by the lasso (cross validation criterion).
## 
##     Estimate
## Trt   0.3389
## V2    1.3120
## V11  -0.8576
## V17  -0.3681
## V32   0.2421
## V35   0.3570
## V39  -0.1401
## V40   0.0275
## V45   0.0945
## V50  -0.0422
```

### Use repeated train and test splitting to estimate subgroup treatment effects:

```r
val.model <- validate.subgroup(subgrp.model, B = 100,
                               method = "training_test",
                               train.fraction = 0.75)
```

### Display estimated subgroup treatment effects:

```r
print(val.model, digits = 2, sample.pct = TRUE)
```

```
## family:  gaussian 
## loss:    sq_loss_lasso 
## method:  weighting 
## 
## validation method:  training_test_replication 
## iterations:  100 
## 
## Average Test Set Outcomes:
##                         Recommended Ctrl            Recommended Trt
## Received Ctrl  -9.18 (SE = 6.83, 20.74%) -19.97 (SE = 6.14, 25.81%)
## Received Trt  -16.35 (SE = 5.66, 24.18%) -13.03 (SE = 5.67, 29.26%)
## 
## Ctrl effect among recommended Ctrl   Trt effect among recommended Trt 
##           7.17 (SE = 10.5, 44.93%)           6.94 (SE = 8.96, 55.07%)
```

Visualize subgroup-specific treatment effect estimates across training/testing iterations:

```r
plot(val.model)
```

<img src="vignettes/vis_val-1.png" width="75%" style="display: block; margin: auto;" />

### Investigate the marginal characteristics of the two estimated subgroups

Here we only display covariates with a significantly different mean value (at level 0.05)

```r
summarize.subgroups(subgrp.model, p.value = 0.05)
```

```
##     Avg (recom Ctrl) Avg (recom Trt) Ctrl - Trt pval Ctrl - Trt
## V1          0.266492       -0.057653   0.324146       2.147e-01
## V2         -1.899413        1.870326  -3.769739       6.746e-55
## V3          0.033790        0.123058  -0.089268       7.378e-01
## V4          0.134217        0.217828  -0.083611       7.619e-01
## V5         -0.082638       -0.300738   0.218100       3.866e-01
## V6          0.173496       -0.029785   0.203281       4.567e-01
## V7         -0.002880       -0.132264   0.129384       6.261e-01
## V8          0.088475       -0.062135   0.150610       5.740e-01
## V9         -0.038066       -0.242125   0.204059       4.478e-01
## V10        -0.119139        0.013597  -0.132736       6.219e-01
## V11         1.239408       -1.121597   2.361005       1.325e-21
## V12         0.200841        0.104925   0.095917       7.264e-01
## V13         0.070726       -0.265149   0.335875       2.235e-01
## V14        -0.014771       -0.183628   0.168857       5.129e-01
## V15         0.112392        0.206887  -0.094494       7.140e-01
## V16        -0.225112       -0.029797  -0.195315       4.836e-01
## V17         0.956636       -0.652966   1.609602       3.814e-09
## V18        -0.161601        0.102548  -0.264149       3.517e-01
## V19        -0.114320       -0.065134  -0.049186       8.524e-01
## V20        -0.301426       -0.094245  -0.207180       4.323e-01
## V21         0.151313        0.099241   0.052072       8.421e-01
## V22         0.036973        0.025914   0.011059       9.684e-01
## V23        -0.061897        0.383565  -0.445462       8.684e-02
## V24         0.059496       -0.083219   0.142715       5.985e-01
## V25        -0.124776       -0.333858   0.209082       4.444e-01
## V26        -0.171079       -0.126225  -0.044854       8.672e-01
## V27         0.179029        0.198470  -0.019441       9.426e-01
## V28         0.017898        0.030955  -0.013057       9.612e-01
## V29         0.006661        0.179831  -0.173170       5.230e-01
## V30        -0.021850        0.042430  -0.064280       8.090e-01
## V31         0.427910       -0.060331   0.488241       7.183e-02
## V32        -0.499087        0.174704  -0.673791       9.574e-03
## V33        -0.161684       -0.138999  -0.022685       9.342e-01
## V34        -0.031119       -0.286260   0.255141       3.502e-01
## V35        -0.717008        0.464186  -1.181194       1.356e-05
## V36         0.017353        0.051908  -0.034555       8.970e-01
## V37        -0.086710       -0.128919   0.042209       8.697e-01
## V38         0.003970       -0.117545   0.121515       6.553e-01
## V39         0.309017       -0.355246   0.664264       1.700e-02
## V40        -0.428111       -0.024037  -0.404074       1.226e-01
## V41         0.171776       -0.028166   0.199942       4.722e-01
## V42        -0.257506       -0.114546  -0.142960       5.872e-01
## V43        -0.293105        0.378708  -0.671814       1.706e-02
## V44        -0.119815        0.036825  -0.156640       5.348e-01
## V45        -0.274946        0.042420  -0.317366       2.423e-01
## V46         0.053626       -0.031523   0.085149       7.542e-01
## V47        -0.027894        0.004777  -0.032671       9.030e-01
## V48         0.109658        0.169670  -0.060012       8.180e-01
## V49        -0.086639       -0.083815  -0.002824       9.920e-01
## V50         0.049473       -0.279293   0.328767       2.440e-01
##     SE (recom Ctrl) SE (recom Trt)
## V1           0.1874         0.1815
## V2           0.1469         0.1536
## V3           0.1800         0.1966
## V4           0.1937         0.1963
## V5           0.1780         0.1780
## V6           0.2019         0.1836
## V7           0.1869         0.1883
## V8           0.1827         0.1956
## V9           0.1901         0.1898
## V10          0.1935         0.1868
## V11          0.1753         0.1577
## V12          0.1975         0.1898
## V13          0.1933         0.1965
## V14          0.1758         0.1887
## V15          0.1750         0.1891
## V16          0.2053         0.1883
## V17          0.1989         0.1800
## V18          0.2111         0.1891
## V19          0.1844         0.1893
## V20          0.1913         0.1813
## V21          0.1924         0.1767
## V22          0.1896         0.2045
## V23          0.1818         0.1853
## V24          0.1928         0.1902
## V25          0.1931         0.1932
## V26          0.1999         0.1786
## V27          0.1930         0.1887
## V28          0.1797         0.1993
## V29          0.1861         0.1970
## V30          0.1923         0.1835
## V31          0.1986         0.1839
## V32          0.1715         0.1941
## V33          0.1986         0.1899
## V34          0.1975         0.1883
## V35          0.1825         0.1973
## V36          0.1918         0.1853
## V37          0.1699         0.1930
## V38          0.1972         0.1875
## V39          0.1906         0.2016
## V40          0.1791         0.1902
## V41          0.1997         0.1933
## V42          0.1939         0.1779
## V43          0.2091         0.1872
## V44          0.1756         0.1810
## V45          0.1870         0.1963
## V46          0.1871         0.1972
## V47          0.1857         0.1933
## V48          0.1769         0.1914
## V49          0.2010         0.1953
## V50          0.1952         0.2033
```

## Accessing Help Files for Main Functions of `personalized`

Access help files for the main functions of the `personalized` package:

```r
?fit.subgroup
?validate.subgroup
```


