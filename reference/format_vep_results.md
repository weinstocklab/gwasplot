# Format VEP API results into a tidy data.frame

Format VEP API results into a tidy data.frame

## Usage

``` r
format_vep_results(vep_results, flag_pick = FALSE)
```

## Arguments

- vep_results:

  Raw results from VEP API (can be data.frame or list)

- flag_pick:

  Logical. If TRUE, filter for transcripts with PICK flag

## Value

Formatted data.frame
