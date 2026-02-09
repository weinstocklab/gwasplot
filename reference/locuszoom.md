# LocusZoom-style plot for GWAS results

LocusZoom-style plot for GWAS results

## Usage

``` r
locuszoom(
  x,
  locus_chr,
  locus_start,
  locus_end,
  include_ccres = FALSE,
  ccre_biosample = NULL,
  ccre_cell_type = NULL,
  valid_biotype = "protein_coding",
  ...
)

# S3 method for class 'GWASFormatter'
locuszoom(
  x,
  locus_chr,
  locus_start,
  locus_end,
  include_ccres = FALSE,
  ccre_biosample = NULL,
  ccre_cell_type = NULL,
  valid_biotype = "protein_coding",
  ...
)

# S3 method for class 'data.frame'
locuszoom(
  x,
  locus_chr,
  locus_start,
  locus_end,
  include_ccres = FALSE,
  ccre_biosample = NULL,
  ccre_cell_type = NULL,
  valid_biotype = "protein_coding",
  ...
)
```

## Arguments

- x:

  A GWASFormatter object or a data.frame/tibble containing GWAS results.

- locus_chr:

  Chromosome of the locus to plot (e.g., 'chr1').

- locus_start:

  Start position of the locus (integer).

- locus_end:

  End position of the locus (integer).

- include_ccres:

  Logical, whether to include ENCODE SCREEN cCREs track. Default is
  FALSE.

- ccre_biosample:

  Optional biosample ID for biosample-specific cCRE filtering (e.g.,
  "GM12878_ENCDO000AAK").

- ccre_cell_type:

  Optional cell type name for cCRE filtering (e.g., "GM12878"). Will be
  converted to biosample ID.

- ...:

  Additional arguments passed to methods.

## Value

A ggplot2 object with the locuszoom-style plot.

## Methods (by class)

- `locuszoom(GWASFormatter)`: Method for GWASFormatter objects

- `locuszoom(data.frame)`: Method for data.frame/tibble objects
