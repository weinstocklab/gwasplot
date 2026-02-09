# Append lfsr results from ashr to a data frame

This function appends lfsr results from the ashr package to a data frame
containing GWAS summary statistics.

## Usage

``` r
append_ashr_results(x, ...)
```

## Arguments

- x:

  A data frame containing GWAS summary statistics with columns "BETA"
  and "SE".

- ...:

  Additional arguments passed to the method.

## Value

A data frame with additional columns for lfsr, qvalue, pm, and psd.
