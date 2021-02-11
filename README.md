
<!-- README.md is generated from README.Rmd. Please edit that file -->

# staggered

<!-- badges: start -->

<!-- badges: end -->

This packages computes the efficient estimator for settings with
randomized treatment timing, based on the theoretical results in [Roth
and Sant’Anna (2021)](https://arxiv.org/pdf/2102.01291.pdf). If units
are randomly (or quasi-randomly) assigned to begin treatment at
different dates, the efficient estimator can potentially offer
substantial gains over methods that only impose parallel trends. The
package also allows for calculating the estimator of [Callaway and
Sant’Anna
(2020)](https://www.sciencedirect.com/science/article/pii/S0304407620303948?dgcid=author)
and the simple-difference-in-means as special cases.

## Installation

You can install the package from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jonathandroth/staggered")
```

## Example

We now illustrate how to use the package by re-creating some of the
results in the application section of [Roth and Sant’Anna
(2021)](https://arxiv.org/pdf/2102.01291.pdf). Our data contains a
balanced panel of police officers in Chicago who were randomly given a
procedural justice training on different dates.

### Loading the package and the data

We first load the staggered package as well as some auxiliary packages
for modifying and plotting the results.

``` r
library(staggered) #load the staggered package
library(dplyr) #load dplyr for data manipulation
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(ggplot2) #load ggplot2 for plotting the results
library(purrr)

df <- staggered::pj_officer_level_balanced #load the officer data
```

We modify the data so that the time dimension is named t, the period of
treatment is named g, the outcome is named y, and the individual
identifiers are named i.

``` r
df <- df %>% rename(t = period, y = complaints, g = first_trained, i = uid)
```

### Simple aggregate parameters

We now can call the function calculate\_adjusted\_estimator\_and\_se to
calculate the efficient estimator. With staggered treatment timing,
there are several ways to aggregate treatment effects across cohorts and
time periods. The following block of code calculates the simple,
calendar-weighted, and cohort-weighted average treatment effects (see
Section 3.2 of [Roth and Sant’Anna
(2021)](https://arxiv.org/pdf/2102.01291.pdf) for more about different
aggregation schemes).

``` r
#Calculate efficient estimator for the simple weighted average
staggered(df = df, estimand = "simple")
#>       estimate          se   se_neyman
#> 1 -0.001126981 0.002115194 0.002119248
```

``` r
#Calculate efficient estimator for the cohort weighted average
staggered(df = df, estimand = "cohort") 
#>       estimate          se   se_neyman
#> 1 -0.001084689 0.002261011 0.002264876
```

``` r
#Calculate efficient estimator for the calendar weighted average
staggered(df = df, estimand = "calendar")
#>      estimate         se   se_neyman
#> 1 -0.00187198 0.00255863 0.002561472
```

### Event-study Parameters

We can also calculate an \`\`event-study’’ that computes the
average-treatment effect at each lag since treatment.

``` r
#Calculate event-study coefficients for the first 24 months (month 0 is instantaneous effect)
eventPlotResults <- staggered(df = df, estimand = "eventstudy", eventTime = 0:23)
  
eventPlotResults %>% head()
#>        estimate          se   se_neyman eventTime
#> 1  3.083575e-04 0.002645327 0.002650957         0
#> 2  2.591678e-03 0.002614563 0.002621513         1
#> 3 -4.872562e-05 0.002622640 0.002623634         2
#> 4  2.043434e-03 0.002715695 0.002720467         3
#> 5  2.977076e-03 0.002653917 0.002659630         4
#> 6  7.979656e-04 0.002721784 0.002727140         5
```

``` r
#Create event-study plot from the results of the event-study
eventPlotResults %>% 
    mutate(ymin_ptwise = estimate + 1.96*se,
           ymax_ptwise = estimate - 1.96*se)%>%
  ggplot(aes(x=eventTime, y =estimate)) +
  geom_pointrange(aes(ymin = ymin_ptwise, ymax = ymax_ptwise))+ 
  geom_hline(yintercept =0) +
  xlab("Event Time") + ylab("Estimate") +
  ggtitle("Effect of Procedural Justice Training on Officer Complaints")
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

### Other Estimators

If instead of the plug-in efficient estimator, one wishes to calculate
the simple difference-in-means or Callaway and Sant’Anna estimator, one
can specify the argument beta=0 or beta=1, respectively, to the
functions above.
