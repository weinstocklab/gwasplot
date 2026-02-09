# Rename columns in the dataset to standard names

Rename columns in the dataset to standard names

## Usage

``` r
reformat_names(ds, format)
```

## Arguments

- ds:

  A data frame or tibble or duckdb tbl containing the GWAS summary
  statistics.

- format:

  The format of the input dataset (either "saige" or "regenie").

## Value

An object of the same class as \`ds\` with renamed columns.
