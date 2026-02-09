# Pull top hits from a GWASFormatter or a data.frame/tibble

Pull top hits from a GWASFormatter or a data.frame/tibble

## Usage

``` r
select_top_hits(x, threshold = 5e-08, ...)

# S3 method for class 'GWASFormatter'
select_top_hits(x, threshold = 5e-08, ...)

# S3 method for class 'data.frame'
select_top_hits(x, threshold = 5e-08, ...)
```

## Arguments

- x:

  A GWASFormatter object, data.frame, or tibble.

- threshold:

  The p-value threshold to filter the top hits. Default is 5e-8.

- ...:

  Additional arguments (unused).

## Value

For GWASFormatter, a tibble of filtered hits; for data.frame/tibble, a
filtered data.frame/tibble.

## Methods (by class)

- `select_top_hits(GWASFormatter)`: Method for GWASFormatter objects

- `select_top_hits(data.frame)`: Method for data.frame/tibble objects
