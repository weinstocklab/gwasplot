# QQ plot for p-values

QQ plot for p-values

## Usage

``` r
qqplot(
  pvalues,
  should.thin = TRUE,
  thin.obs.places = 2,
  thin.exp.places = 2,
  xlab = expression(paste("Expected (", -log[10], " p-value)")),
  ylab = expression(paste("Observed (", -log[10], " p-value)")),
  draw.conf = TRUE,
  conf.points = 1000,
  conf.col = "lightgray",
  conf.alpha = 0.05,
  already.transformed = FALSE,
  point.size = 1.5,
  point.alpha = 0.7,
  point.color = "#3366FF",
  ...
)
```

## Arguments

- pvalues:

  A numeric vector of p-values or a list of numeric vectors of p-values.

- should.thin:

  Logical indicating whether to thin the points for plotting. Default is
  TRUE.

- thin.obs.places:

  Number of decimal places to round the observed p-values for thinning.
  Default is 2.

- thin.exp.places:

  Number of decimal places to round the expected p-values for thinning.
  Default is 2.
