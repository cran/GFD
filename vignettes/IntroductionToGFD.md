<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{An Introduction to GFD}
-->


## Introduction

This vignette documents the use of the `GFD` function for the analysis of general factorial designs. The `GFD` function
calculates the Wald-type statistic (WTS), the ANOVA-type statistic (ATS) and a permuted Wald-type statistic (WTPS). 
These test statistics can be used for general factorial designs (crossed or nested) with an arbitrary number of 
factors, unequal covariance matrices among groups and unbalanced data even for small sample sizes.

## Data Example (crossed design)

For illustration purposes, we will use the data set `pizza` which is included in the `GFD` package. We first load the
`GFD` package and the data set.


```r
library(GFD)
data(pizza)
```

The objective of the study was to see how the delivery time in minutes would be affected by
three different factors: whether thick or thin crust was ordered (factor A), whether Coke was
ordered with the pizza or not (factor B), and whether or not garlic bread was ordered as a
side (factor C). 


```r
head(pizza)
```

```
##   Crust Coke Bread Driver  Hour Delivery
## 1  thin  yes   yes      M 20.87       14
## 2 thick  yes    no      M 20.78       21
## 3  thin   no    no      M 20.75       18
## 4  thin   no   yes      F 20.60       17
## 5 thick   no    no      M 20.70       19
## 6 thick   no   yes      M 20.95       17
```

This is a three-way crossed design, where each factor has two levels. We will now analyze this design with the `GFD` function. 
The `GFD` function takes as arguments:
* `formula`: A formula consisting of the outcome variable on the left hand side of a \~ operator and the factor
variables of interest on the right hand side. An interaction term must be specified.
* `data`: A data.frame, list or environment containing the variables in `formula`.
* `nperm`: The number of permutations. Default value is 10000.
* `alpha`: The significance level, default is 0.05.


```r
set.seed(1234)
model1 <- GFD(Delivery ~ Crust * Coke * Bread, data = pizza, nperm = 10000, alpha = 0.05)
summary(model1)
```

```
## Call: 
## Delivery ~ Crust * Coke * Bread
## 
## Descriptive:
##   Crust Coke Bread n Means Variances Lower 95 % CI Upper 95 % CI
## 1  thin   no    no 2  19.0       2.0      14.69735      23.30265
## 5  thin   no   yes 2  17.5       0.5      15.34867      19.65133
## 3  thin  yes    no 2  17.5       4.5      11.04602      23.95398
## 7  thin  yes   yes 2  15.0       2.0      10.69735      19.30265
## 2 thick   no    no 2  19.5       0.5      17.34867      21.65133
## 6 thick   no   yes 2  18.0       2.0      13.69735      22.30265
## 4 thick  yes    no 2  21.5       0.5      19.34867      23.65133
## 8 thick  yes   yes 2  18.5       0.5      16.34867      20.65133
## 
## Wald-Type Statistic (WTS):
##                  Test statistic df      p-value p-value WTPS
## Crust                     11.56  1 0.0006738585       0.0089
## Coke                       0.36  1 0.5485062355       0.5613
## Crust:Coke                 6.76  1 0.0093223760       0.0286
## Bread                     11.56  1 0.0006738585       0.0073
## Crust:Bread                0.04  1 0.8414805811       0.8153
## Coke:Bread                 1.00  1 0.3173105079       0.3457
## Crust:Coke:Bread           0.04  1 0.8414805811       0.8212
## 
## ANOVA-Type Statistic (ATS):
##                  Test statistic df1      df2    p-value
## Crust                     11.56   1 4.699248 0.02121110
## Coke                       0.36   1 4.699248 0.57625702
## Crust:Coke                 6.76   1 4.699248 0.05122842
## Bread                     11.56   1 4.699248 0.02121110
## Crust:Bread                0.04   1 4.699248 0.84984482
## Coke:Bread                 1.00   1 4.699248 0.36598284
## Crust:Coke:Bread           0.04   1 4.699248 0.84984482
```

The output consists of three parts: `model1$Descriptive` gives an overview of the descriptive statistics: The number of observations, 
mean and variance as well as confidence intervals (based on quantiles of the t-distribution) are displayed for each factor level combination.
`model1$WTS` contains the results for the Wald-type test: The test statistic, degree of freedom and p-values based on the asymptotic $\chi^2$ distribution
and the permutation procedure, respectively, are displayed. Note that the $\chi^2$ approximation is very liberal for small sample sizes and therefore the WTPS is 
recommended for such situations. Finally, `model1$ATS` contains the corresponding results based on the ATS. This test statistic tends to rather
conservative decisions in the case of small sample sizes and is even asymptotically only an approximation, thus not providing an asymptotic level $\alpha$ test.

We find a significant influence of the factors Crust and Bread. The WTS and WTPS also
suggest a significant interaction between the factors Crust and Coke at 5% level, which is only
borderline significant when using the ATS.

## Data example (nested design)

Nested designs can also be analyzed using the `GFD` function. We consider the data set `curdies` from the `GFD` package:


```r
data("curdies")
set.seed(987)
nested <- GFD(dugesia ~ season + season:site, data = curdies)
summary(nested)
```

```
## Call: 
## dugesia ~ season + season:site
## 
## Descriptive:
##         season site n     Means  Variances Lower 95 % CI Upper 95 % CI
## 1 SUMMER          4 6 0.4190947  0.4615290   -0.25954958     1.0977390
## 2 SUMMER          5 6 0.2290862  0.3148830   -0.33146759     0.7896401
## 3 SUMMER          6 6 0.1942443  0.0729142   -0.07549781     0.4639864
## 4 WINTER          1 6 2.0494375  4.0647606    0.03543415     4.0634408
## 5 WINTER          2 6 4.1819078 35.6801853   -1.78509515    10.1489107
## 6 WINTER          3 6 0.6782063  0.1910970    0.24151987     1.1148927
## 
## Wald-Type Statistic (WTS):
##             Test statistic df    p-value p-value WTPS
## season            5.415180  1 0.01996239       0.0001
## season:site       5.200991  4 0.26728919       0.3154
## 
## ANOVA-Type Statistic (ATS):
##             Test statistic      df1      df2    p-value
## season            5.415180 1.000000 6.447707 0.05593278
## season:site       1.382224 1.217424 6.447707 0.29278958
```

The aim of the study was to describe basic patterns of variation in a small flatworm, Dugesia, in the Curdies
River, Western Victoria. Therefore, worms were sampled at two different seasons and three
different sites within each season. For our analyses we consider both factors as fixed (e.g.,
some sites may only be accessed in summer). In this setting, both WTS and WTPS detect a significant influence of the season whereas
the ATS, again, only shows a borderline significance at 5% level. The effect of the site is not
significant.

## Plotting

The `GFD` package is equipped with a plotting function, displaying the calculated means along with $(1-\alpha)$ confidence intervals.
The `plot` function takes a `GFD` object as an argument. In addition, the factor of interest may be specified. If this argument is 
omitted in a two- or higher-way layout, the user is asked to specify the factor for plotting. Furthermore, additional graphical parameters
can be used to customize the plots. The optional argument `legendpos` specifies the position of the legend in higher-way layouts.


```r
plot(model1, factor = "Crust:Coke:Bread", legendpos = "center", main = "Delivery time of pizza", xlab = "Bread")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

```r
plot(model1, factor = "Crust:Coke", legendpos = "topleft", main = "Two-way interaction", xlab = "Coke", col = 3:5, pch = 17)
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-2.png)

```r
plot(nested, factor = "season:site", xlab = "site")
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-3.png)

## optional GUI

The `GFD` package is equipped with an optional graphical user interface, which is based on `RGtk2`. The GUI may be started in `R` (if `RGtk2` is installed) using the
command `calculateGUI()`. 


```r
calculateGUI()
```

The user can specify the data location
(either directly or via the "load data" button), the formula, the number of permutations and
the significance level. Additionally, one can specify whether or not headers are
included in the data file, and which separator and character symbols are used for decimals
in the data file. The GUI also provides a plotting option, which generates a new window
for specifying the factors to be plotted (in higher-way layouts) along with a few plotting
parameters.
