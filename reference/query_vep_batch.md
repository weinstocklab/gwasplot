# Query VEP API for a batch of variants

Query VEP API for a batch of variants

## Usage

``` r
query_vep_batch(
  variant_batch,
  flag_pick = FALSE,
  include_hgvs = TRUE,
  include_domains = FALSE,
  include_regulatory = FALSE,
  verbose = TRUE
)
```

## Arguments

- variant_batch:

  Character vector of VCF-formatted variants

- flag_pick:

  Include flag_pick to mark the selected transcript

- include_hgvs:

  Include HGVS nomenclature

- include_domains:

  Include protein domains

- include_regulatory:

  Include regulatory consequences

## Value

data.frame with VEP results
