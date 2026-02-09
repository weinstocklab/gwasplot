# Annotate variants with functional consequences using Ensembl VEP API

This function queries the Ensembl Variant Effect Predictor (VEP) REST
API to annotate a list of variants with their functional consequences.
It processes variants in batches to respect API rate limits and returns
detailed consequence information including transcript effects, protein
changes, and regulatory impacts.

## Usage

``` r
annotate_variants_ensembl(
  variants,
  batch_size = 200,
  flag_pick = TRUE,
  include_hgvs = TRUE,
  include_domains = FALSE,
  include_regulatory = FALSE,
  sleep_time = 1,
  verbose = FALSE
)
```

## Arguments

- variants:

  A character vector containing variant information in the format
  "chr_pos_ref_alt" (e.g., "chr21_26960070_G_A" or "21_26960070_G_A").

- batch_size:

  Integer specifying the number of variants to process per API call.
  Default is 200 (max 1000). Smaller batches are more reliable for large
  datasets.

- flag_pick:

  Logical. If TRUE, only report the transcript with the PICK flag.
  Default is TRUE to get the single best transcript per variant.

- include_hgvs:

  Logical. If TRUE, include HGVS nomenclature in results. Default is
  TRUE.

- include_domains:

  Logical. If TRUE, include protein domain information. Default is
  FALSE.

- include_regulatory:

  Logical. If TRUE, include regulatory feature consequences. Default is
  FALSE.

- sleep_time:

  Numeric. Time in seconds to wait between API calls to respect rate
  limits. Default is 1 second. Increase if you encounter rate limiting.

- verbose:

  Logical. If TRUE, print progress messages. Default is TRUE.

## Value

A tibble containing variant consequences with the following columns:

- ID - Original variant string

- CHROM - Chromosome (extracted from variant)

- POS - Position (extracted from variant)

- REF - Reference allele (extracted from variant)

- ALT - Alternate allele (extracted from variant)

- most_severe_consequence - Most severe consequence for the variant

- gene_id - Ensembl gene ID

- gene_symbol - Gene symbol

- transcript_id - Ensembl transcript ID

- consequence_terms - All consequence terms (comma-separated)

- impact - Impact level (HIGH, MODERATE, LOW, MODIFIER)

- protein_position - Position in protein sequence

- amino_acids - Reference and alternate amino acids

- codons - Reference and alternate codons

- existing_variation - Known variant IDs (e.g., rsIDs)

- hgvsc - HGVS coding sequence nomenclature (if include_hgvs=TRUE)

- hgvsp - HGVS protein nomenclature (if include_hgvs=TRUE)

- domains - Protein domains affected (if include_domains=TRUE)

## Details

The function handles variant input in the format "chr_pos_ref_alt":

- Input format: "chr21_26960070_G_A" or "21_26960070_G_A"

- Automatic rate limiting to respect Ensembl's 15 requests/second limit

- Batch processing for efficient handling of large variant lists

- Error handling and retry logic for failed requests

## Note

- Respects Ensembl API rate limits (15 requests/second)

- Large variant lists are automatically batched

- Requires internet connection to Ensembl REST API

- For very large datasets (\>10,000 variants), consider using local VEP
  installation

- Only supports human variants (homo_sapiens)

## References

McLaren et al. (2016). The Ensembl Variant Effect Predictor. Genome
Biology 17, 122. doi:10.1186/s13059-016-0974-4

## See also

<https://rest.ensembl.org/documentation/info/vep_region_post>
<https://github.com/Ensembl/ensembl-vep>

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Single variant
variants <- "chr21_26960070_G_A"
consequences <- annotate_variants_ensembl(variants)

# Example 2: Multiple variants
variants <- c(
  "chr21_26960070_G_A",
  "21_26965148_G_A",
  "chrX_155066068_C_T"
)
consequences <- annotate_variants_ensembl(variants)

# Example 3: With additional options
consequences <- annotate_variants_ensembl(
  variants,
  flag_pick = TRUE,
  include_domains = TRUE,
  batch_size = 100
)
} # }
```
