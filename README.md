# Overview

**gwasplot** provides fast and efficient tools for visualizing and annotating GWAS 
summary statistics. It offers functions to reformat data, create plots, and annotate 
top hits with genomic features. 

It is designed to work with output from GWAS performed on WGS data with many rare variants, using [duckdb](https://duckdb.org/) to manipulate data. 

For more information, visit the [GitHub repository](https://github.com/weinstockj/gwasplot) or the [project website](https://weinstockj.github.io/gwasplot/).

## Roadmap
 - add LD into locuszoom plots, perhaps referencing this [setup](https://github.com/weinstockj/LD_REFERENCE_PANEL).
 - build out OpenTargets API functionality for annotations
 - add in helper functions to faciliate fine-mapping with other packages (e.g., SuSIE)

# Installation

Install the development version of **gwasplot** from GitHub:

```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install gwasplot from GitHub
remotes::install_github("weinstockj/gwasplot")
```


## Input formats

`gwasplot` expects output from either regenie or saige in csv or parquet format (parquet is generally recommended). 

Regenie output looks like this:

```r
dplyr::glimpse(df)
Rows: 5
Columns: 12
$ CHROM     <chr> "chr3", "chr3", "chr3", "chr3", "chr3"
$ POS       <int> 90000011, 90000261, 90000274, 90000324, 90000366
$ ID        <chr> "chr3-90000011-T-C", "chr3-90000261-GATT-G", "chr3-900002…"
$ ALLELE0   <chr> "T", "GATT", "T", "C", "T"
$ ALLELE1   <chr> "C", "G", "A", "A", "G"
$ A1FREQ    <dbl> 1.37775e-04, 2.91090e-04, 2.38259e-04, 1.62637e-04, 4.66158e…
$ N         <int> 482670, 482669, 482669, 482669, 482669
$ BETA      <dbl> -0.00157451, -0.06463220, -0.06452010, 0.02796820, -0.021063…
$ SE        <dbl> 0.0730887, 0.0479369, 0.0540183, 0.0646778, 0.1213880
$ CHISQ     <dbl> 0.000464077, 1.817850000, 1.426620000, 0.186990000, 0.030108…
$ LOG10P    <dbl> 0.00752913, 0.75063100, 0.63391900, 0.17689500, 0.06437000
```

gwastools will read in data from either regenie or saige and standarize the column names,
and then create a table called `summary_stats` in a `duckdb` database. 

# Key Features

## Data Quality Control and Filtering

**gwasplot** includes robust filtering capabilities to remove variants from problematic genomic regions that can lead to spurious associations:

### Exclude Difficult Regions
The `exclude_difficult_regions()` function removes variants from regions with known mapping difficulties:

```r
# Remove variants from multiple types of difficult regions
gwas_filtered <- exclude_difficult_regions(gwas_stats, 
                                         beds_exclude = c("hg19diff", "UCSC_unusual", "GRC_exclusions"))
```

Available region types include:
- **hg19diff**: Regions with assembly differences between genome builds
- **UCSC_unusual**: Unusual regions identified by UCSC (gaps, heterochromatin, etc.)
- **GRC_exclusions**: Genome Reference Consortium exclusion regions
- **GIAB_difficult_regions**: Genome in a Bottle difficult-to-map regions

### Custom Variant Filtering
Filter variants using custom whitelists or blacklists:

```r
# Keep only variants in a specific region or set
filtered_gwas <- filter_variants(gwas_stats, subset = "my_variants.parquet")

# Exclude specific variants
filtered_gwas <- filter_variants(gwas_stats, exclude = "bad_variants.parquet")
```

## Gene Annotation

Annotate variants with nearest protein-coding genes:

```r
# Add gene annotations
annotated_gwas <- find_nearest_gene(gwas_stats, threshold = 1e5)  # 100kb window
```

Additional specialized annotations:
- **CHIP genes**: `annotate_with_chip_genes()` - flags variants in clonal hematopoiesis genes
- **Centromeres**: `annotate_with_centromere()` - identifies variants in centromeric regions  
- **Immunoglobulin loci**: `annotate_with_immunoglobulin()` - flags variants in Ig heavy/light chain regions

## Visualization

### Manhattan Plots
```r
manhattan(gwas = gwas_stats, output_prefix = "my_manhattan")
```

### QQ Plots
```r
qqplot_save(gwas = gwas_stats, output_prefix = "my_qqplot")
```

### LocusZoom Plots
Create detailed regional association plots with gene tracks:

```r
# Plot a specific locus with gene structure
locuszoom(gwas_stats, locus_chr = "chr2", locus_start = 25000000, locus_end = 25500000)
```

### Volcano Plots
For post-GWAS analysis (e.g., after empirical Bayes shrinkage):

```r
# Automatically uses lfsr if available, otherwise p-values
volcano(gwas_data, phenotype_label = "My Phenotype")
```

# Usage
Below are some sample usage examples:

## Basic Workflow

Reformat GWAS Summary Statistics
Use the reformat_summary_statistics function to read and reformat GWAS summary statistics from a parquet or CSV file:

```r
library(gwasplot)

# Reformat summary statistics from a file
gwas_stats <- reformat_summary_statistics("path/to/your/file.parquet")
print(gwas_stats)

GWAS object
File path: ../concatenated_results.parquet
Detected format: regenie
Data names: CHROM, POS, ID, ALLELE0, ALLELE1, A1FREQ, N, BETA, SE, CHISQ, LOG10P, phenotype
Data preview:
# A tibble: 5 × 9
  CHROM      POS REF   ALT      AF_ALT     BETA  LOG10P PVALUE ID
  <chr>    <dbl> <chr> <chr>     <dbl>    <dbl>   <dbl>  <dbl> <chr>
1 chr3  90000011 T     C     0.000138  -0.00157 0.00753  0.983 chr3_90000011_T_C
2 chr3  90000261 GATT  G     0.000291  -0.0646  0.751    0.178 chr3_90000261_GA…
3 chr3  90000274 T     A     0.000238  -0.0645  0.634    0.232 chr3_90000274_T_A
4 chr3  90000324 C     A     0.000163   0.0280  0.177    0.665 chr3_90000324_C_A
5 chr3  90000366 T     G     0.0000466 -0.0211  0.0644   0.862 chr3_90000366_T_G
```

## Complete Analysis Pipeline

```r
# Load and reformat data
gwas_stats <- reformat_summary_statistics("path/to/gwas_results.parquet")

# Quality control: remove problematic regions for WGS GWAS
gwas_clean <- gwas_stats %>%
  exclude_difficult_regions(
    beds_exclude = c(
      "hg19diff", 
      "UCSC_unusual", 
      "GRC_exclusions"
      )
    )


manhattan(gwas_clean, output_prefix = "my_analysis")
qqplot(gwas_clean, output_prefix = "my_analysis")

# Get top hits and annotate with genes
top_hits <- gwas_clean %>%
  select_top_hits(threshold = 5e-8) %>%
  find_nearest_gene()

# Plot specific loci
locuszoom(gwas_clean, locus_chr = "chr2", locus_start = 25000000, locus_end = 25500000)
```

## Contact

Email Josh Weinstock. 
See our other work [here](https://weinstocklab.org).