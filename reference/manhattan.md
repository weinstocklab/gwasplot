# Plot a Manhattan plot from a gwas object

Plot a Manhattan plot from a gwas object

## Usage

``` r
manhattan(gwas, output, lower_logp_threshold = 3, ...)
```

## Arguments

- gwas:

  A gwas object containing the data to plot.

- output:

  The output file name.

- lower_logp_threshold:

  The lower threshold for the -log10(p-value) to plot. Default is 3.0.

- ...:

  Additional arguments passed to \`ggsave\`.
