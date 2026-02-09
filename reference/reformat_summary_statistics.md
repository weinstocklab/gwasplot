# Reformat GWAS summary statistics from regenie or SAIGE

This function takes a file path to a parquet or CSV file of GWAS summary
statistics output from regenie or SAIGE and reformats it to have a
standard set of columns.

## Usage

``` r
reformat_summary_statistics(file_path, use_cache = FALSE, read_only = FALSE)
```

## Arguments

- file_path:

  Path to the GWAS summary statistics file (parquet or CSV).

- use_cache:

  Logical. If TRUE, uses the cached table 'summary_stats' instead of
  reading the file again.

- read_only:

  Logical. If TRUE, opens the database connection in read-only mode.

## Value

An R6 object of class \`GWASFormatter\` containing the reformatted
summary statistics.
