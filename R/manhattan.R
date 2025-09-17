
# Global variables to avoid NSE warnings in R CMD check
utils::globalVariables(c(
  "CHROM", "POS", "PVALUE", "CHROM_index", "chr_len", "tot", "POScum",
  "center", "aes"
))

# Helper function to create chromosome lookup table
create_chrom_lookup = function() {
  tibble::tibble(
    CHROM = c(glue::glue("chr{1:22}"), "chrX"),
    CHROM_index = 1:23
  )
}

# Helper function to prepare GWAS data for plotting
prepare_manhattan_data = function(gwas_data, chrom_lookup, lower_logp_threshold = 3.0) {
  scaling = 1e8
  
  # Compute chromosome positions and cumulative coordinates
  prepared = gwas_data %>%
    dplyr::select(CHROM, POS, PVALUE) %>%
    dplyr::inner_join(chrom_lookup, by = "CHROM") %>%
    # Compute chromosome size
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarise(chr_len = (max(POS) - min(POS)) / scaling) %>%
    dplyr::ungroup(.) %>%
    # Calculate cumulative start position of each chromosome
    dplyr::arrange(CHROM_index) %>%
    dplyr::mutate(tot = cumsum(chr_len) - chr_len) %>%
    # Add this info to the initial dataset
    dplyr::right_join(gwas_data %>% dplyr::select(CHROM, POS, PVALUE), by = "CHROM") %>%
    # Add a cumulative position of each SNP
    dplyr::arrange(CHROM_index, POS) %>%
    dplyr::mutate(POScum = POS / scaling + tot) %>%
    # Filter SNP to make the plot lighter
    dplyr::filter(-log10(PVALUE) > lower_logp_threshold)
  
  return(prepared)
}

# Helper function to calculate axis positions
calculate_axis_positions = function(prepared_data) {
  axis = prepared_data %>%
    dplyr::group_by(CHROM_index, CHROM) %>%
    dplyr::summarize(center = (max(POScum) + min(POScum)) / 2) %>% # integer overflow concern
    dplyr::ungroup(.) %>%
    dplyr::arrange(CHROM_index, center)
  
  return(axis)
}

# Helper function to create the Manhattan plot
create_manhattan_plot = function(prepared_data, axis_data, pvalue_threshold = 5e-8, y_max_cap = 300) {
  p = ggplot2::ggplot(prepared_data, aes(x=POScum, y=-log10(PVALUE))) +
    # Show all points
    ggrastr::rasterise(
      ggplot2::geom_point(aes(color=as.factor(CHROM_index)), alpha=0.8, size=.3),
      dev = "ragg_png",
      dpi = 300
     ) +
    ggplot2::scale_color_manual(values = c(rep(c("#7173C9", "#01035F"), 11), "#7173C9")) +
    # custom X axis:
    ggplot2::scale_x_continuous(
      label = stringr::str_replace(axis_data$CHROM, "chr", ""),
      breaks= axis_data$center, expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max(pmin(-log10(prepared_data$PVALUE), y_max_cap)), by = 4),
      expand = ggplot2::expansion(mult = c(0, 0), add = c(0, 1))
    ) +     # remove space between plot area and x axis
    # Custom the theme:
    ggplot2::geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "grey") +
    ggplot2::labs(y = expression(-log[10](pvalue))) +
    ggplot2::theme_bw(base_size = 7, base_family = "Helvetica") +
    ggplot2::theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 3),
      axis.text.y = element_text(size = 5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(size = .6)
    )
  
  return(p)
}

# Helper function to save the Manhattan plot
save_manhattan_plot = function(plot, output, ...) {
  log_info("Now rendering")
  
  ggplot2::ggsave(
    output,
    plot,
    ...,
    bg = "white"
  )
  
  log_info("done plotting.")
}

#' Plot a Manhattan plot from a gwas object
#' 
#' @param gwas A gwas object containing the data to plot.
#' @param output The output file name.
#' @param lower_logp_threshold The lower threshold for the -log10(p-value) to plot. Default is 3.0.
#' @param ... Additional arguments passed to `ggsave`.
#' @return NULL
#' @export 
manhattan = function(gwas, output, lower_logp_threshold = 3.0, ...) {
  UseMethod("manhattan")
}

#' @export
manhattan.tbl_df = function(gwas, output, lower_logp_threshold = 3.0, ...) {
  manhattan.data.frame(gwas, output, lower_logp_threshold, ...)
}

#' @export 
manhattan.data.frame = function(gwas, output, lower_logp_threshold = 3.0, ...) {
  log_info("Now preparing to plot")
  
  # Create chromosome lookup
  chrom_lookup = create_chrom_lookup()
  
  # Prepare data
  prepared = prepare_manhattan_data(gwas, chrom_lookup, lower_logp_threshold)
  
  # Calculate axis positions
  axis = calculate_axis_positions(prepared)
  
  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))
  
  # Create plot
  p = create_manhattan_plot(prepared, axis)
  
  # Save plot
  save_manhattan_plot(p, output, ...)
}

#' @export
manhattan.GWASFormatter = function(gwas, output, lower_logp_threshold = 3.0, ...) {
  log_info("Now preparing to plot")
  
  # Create chromosome lookup and copy to database connection
  chrom_lookup = create_chrom_lookup() %>%
    dplyr::copy_to(gwas$con, ., "chrom_lookup", overwrite = TRUE)
  
  # Prepare data
  prepared = prepare_manhattan_data(gwas$data, chrom_lookup, lower_logp_threshold)
  
  # Calculate axis positions and collect from database
  axis = calculate_axis_positions(prepared) %>%
    dplyr::collect(.)
  
  # Collect prepared data from database
  prepared = prepared %>%
    dplyr::collect(.)
  
  log_info(glue("Done preparing to plot {nrow(prepared)} SNPs."))
  
  # Create plot (no y_max_cap for GWASFormatter to avoid pmin issue with database)
  p = create_manhattan_plot(prepared, axis, y_max_cap = max(-log10(prepared$PVALUE)))
  
  # Save plot
  save_manhattan_plot(p, output, ...)
}
