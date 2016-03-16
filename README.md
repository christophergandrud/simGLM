# simGLM

> Simulate and plot quantities of interest from generalised linear
    models using [King, Tomz, and Wittenburg (2000)](http://www.jstor.org/stable/2669316).
    Currently only supports normal linear and logistic regression models.

Christopher Gandrud

[![Build Status](https://travis-ci.org/christophergandrud/simGLM.svg?branch=master)](https://travis-ci.org/christophergandrud/simGLM)

## Example: normal linear regression


```r
library(car) # contains data
library(simGLM)
library(ggplot2) # only needed for adding additional arguments outside of sim_glm

# Estimate model
m1 <- lm(prestige ~ education + type, data = Prestige)

# Create fitted values
fitted_prestige <- expand.grid(education = 6:16, typewc = 1)

# Simulate and plot
sim_glm(obj = m1, newdata = fitted_prestige, x_coef = 'education') +
        ylab('Predicted Job Prestige\n') + xlab('\nYears of Education')
```

```
## typeprof fitted at 0.
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

## Example: logistic regression


```r
# Download data
URL <- 'http://www.ats.ucla.edu/stat/data/binary.csv'
Admission <- read.csv(URL)
Admission$rank <- as.factor(Admission$rank)

# Estimate model
m2 <- glm(admit ~ gre + gpa + rank, data = Admission, family = 'binomial')

# Create fitted values
fitted_admit <- expand.grid(gre = seq(220, 800, by = 10), gpa = c(2, 4), 
                            rank4 = 1)

# Simulate and plot
sim_glm(obj = m2, newdata = fitted_admit, model = 'logit', x_coef = 'gre', 
        group_coef = 'gpa')
```

```
## rank2 fitted at 0.
```

```
## rank3 fitted at 0.
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

## Examples: `bayesglm`

`sim_glm` also works with estimates made using the `bayesglm` funciton in the [arm](https://cran.r-project.org/web/packages/arm/index.html) package.


```r
library(arm)
```


```r
# Estimate model
m2 <- bayesglm(admit ~ gre + gpa + rank, data = Admission, 
               family = binomial(link = 'logit'))

# Simulate and plot
sim_glm(obj = m2, newdata = fitted_admit, model = 'logit', x_coef = 'gre', 
        group_coef = 'gpa')
```

```
## rank2 fitted at 0.
```

```
## rank3 fitted at 0.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

## Install

To install the development version of **simGLM** use:


```r
ghit::install_github('christophergandrud/simGLM')
```
