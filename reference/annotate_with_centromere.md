# Annotate data with centromere information

Annotate data with centromere information

## Usage

``` r
annotate_with_centromere(x, ...)

# S3 method for class 'data.frame'
annotate_with_centromere(x, chrom_col = "CHROM", pos_col = "POS", ...)

# S3 method for class 'tbl_df'
annotate_with_centromere(x, ...)

# S3 method for class 'GWASFormatter'
annotate_with_centromere(x, ...)
```

## Arguments

- x:

  A data frame or tibble containing variant data.

- ...:

  Additional arguments passed to methods.

- chrom_col:

  Name of the column containing chromosome information. Default is
  "CHROM".

- pos_col:

  Name of the column containing position information. Default is "POS".

## Value

A data frame with the centromere information.

## Methods (by class)

- `annotate_with_centromere(data.frame)`: Annotate a data frame or
  tibble with centromere information

- `annotate_with_centromere(tbl_df)`: Alias for the data.frame method

- `annotate_with_centromere(GWASFormatter)`: Centromere annotation
  method for GWASFormatter objects
