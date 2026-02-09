# Filter variants

Filter variants in a GWASFormatter object based on a whitelist or
blacklist of variants.

## Usage

``` r
filter_variants(x, subset = NULL, exclude = NULL, ...)
```

## Arguments

- x:

  A GWASFormatter object.

- subset:

  A file path to a whitelist of variants in Parquet format. Only
  variants in this file will be kept.
