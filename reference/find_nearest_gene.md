# Find the nearest gene for variants

Find the nearest gene for variants

## Usage

``` r
find_nearest_gene(x, ...)

# S3 method for class 'GWASFormatter'
find_nearest_gene(x, threshold = 1e+05, ...)

# S3 method for class 'tbl_df'
find_nearest_gene(
  x,
  threshold = 1e+05,
  chrom_col = "CHROM",
  pos_col = "POS",
  id_col = "ID",
  ...
)
```

## Arguments

- x:

  A gwas object or tibble containing variant data.

- ...:

  Additional arguments passed to methods.

- threshold:

  The distance threshold to consider a gene as nearest. Default is 1e5.

- chrom_col:

  Name of the column containing chromosome information. Default is
  "CHROM".

- pos_col:

  Name of the column containing position information. Default is "POS".

- id_col:

  Name of the column containing variant IDs. Default is "ID".

## Value

The input object with gene annotations added.

## Methods (by class)

- `find_nearest_gene(GWASFormatter)`: Find the nearest gene for each
  variant in a gwas object

- `find_nearest_gene(tbl_df)`: Find the nearest gene for each variant in
  a tibble
