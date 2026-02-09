# Query Open Targets Platform API for Locus-to-Gene (L2G) predictions

This function retrieves L2G predictions from credible sets containing
the specified variant. It replaces the deprecated V2G functionality from
the old Open Targets Genetics API.

## Usage

``` r
query_ot_api_v2g(variant_id = "19_44908822_C_T", pageindex = 0, pagesize = 20)
```

## Arguments

- variant_id:

  A string representing the variant ID in format "CHR_POS_REF_ALT"
  (e.g., "19_44908822_C_T") or an rsID (e.g., "rs123456").

- pageindex:

  Pagination index (currently unused in new API structure).

- pagesize:

  Pagination size (currently unused in new API structure).

## Value

A list containing L2G predictions and associated data from credible
sets.

## Note

This function has been migrated from the deprecated Open Targets
Genetics API to the new Open Targets Platform API. The functionality has
changed from direct variant-to-gene mapping to locus-to-gene predictions
via credible sets.

## Examples

``` r
if (FALSE) { # \dontrun{
  # Query by variant ID
  result <- query_ot_api_v2g("19_44908822_C_T")
  
  # Query by rsID  
  result <- query_ot_api_v2g("rs123456")
} # }
```
