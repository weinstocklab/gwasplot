# Volcano plot for GWAS results

Volcano plot for GWAS results

## Usage

``` r
volcano(x, phenotype_label = NULL, ...)

# S3 method for class 'GWASFormatter'
volcano(x, phenotype_label = NULL, ...)

# S3 method for class 'data.frame'
volcano(x, phenotype_label = NULL, ...)

# S3 method for class 'tbl_df'
volcano(x, phenotype_label = NULL, ...)
```

## Arguments

- x:

  A GWASFormatter object or a data.frame/tibble containing GWAS results.

- phenotype_label:

  Optional label for the plot title.

- ...:

  Additional arguments passed to methods.

## Value

A ggplot2 object with the volcano plot.

## Methods (by class)

- `volcano(GWASFormatter)`: Method for GWASFormatter objects

- `volcano(data.frame)`: Method for data.frame/tibble objects

- `volcano(tbl_df)`: Alias for data.frame method
