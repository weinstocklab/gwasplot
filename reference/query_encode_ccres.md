# Query ENCODE SCREEN cCREs for a genomic region

Query ENCODE SCREEN cCREs for a genomic region

## Usage

``` r
query_encode_ccres(
  locus_chr,
  locus_start,
  locus_end,
  assembly = "grch38",
  biosample = NULL,
  cell_type = NULL
)
```

## Arguments

- locus_chr:

  Chromosome (e.g., "chr1")

- locus_start:

  Start position

- locus_end:

  End position

- assembly:

  Genome assembly ("grch38" or "mm10")

- biosample:

  Optional biosample ID for cell-type-specific filtering (e.g.,
  "GM12878_ENCDO000AAK")

- cell_type:

  Optional cell type name for filtering (e.g., "GM12878"). If provided,
  will be converted to biosample ID.

## Value

Data frame with cCRE information
