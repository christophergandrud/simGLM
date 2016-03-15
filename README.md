# simGLM

> Simulate and plot quantities of interest from generalised linear
    models. Currently only supports normal linear and logistic regression models.

Christopher Gandrud



## Example: normal linear regression


```r
library(car) # Contains data
m1 <- lm(prestige ~ education + type,
         data = Prestige)

fitted_prestige <- expand.grid(education = 6:16, typewc = 1)

sim_glm(obj = m1, newdata = fitted_prestige, x_coef = 'education') +
    ylab('Predicted Job Prestige\n') + xlab('\nYears of Education')
```

```
## typeprof fitted at 0.
```

```
## Error in eval(expr, envir, enclos): could not find function "ylab"
```

## Example: logistis regression


```r
URL <- 'http://www.ats.ucla.edu/stat/data/binary.csv'
Admission <- read.csv(URL)
Admission$rank <- as.factor(Admission$rank)

m2 <- glm(admit ~ gre + gpa + rank,
          data = Admission, family = 'binomial')

fitted_admit <- expand.grid(gre = seq(220, 800, by = 10), gpa = c(2, 4),
                              rank4 = 1)
sim_glm(obj = m2, newdata = fitted_admit, model = 'logit',
        x_coef = 'gre', group_coef = 'gpa')
```

```
## rank2 fitted at 0.
```

```
## rank3 fitted at 0.
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)