#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2022-10-22
## Script to process allelic counts

## Load/install packages
packages = c("data.table", "argparser", "ggplot2", "cowplot", "ggsci", "scales", "patchwork", "gtools")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

theme_set(theme_cowplot())

## Parse arguments
parser = arg_parser("Generate summary table and plots")
parser = add_argument(parser, "--input", short = "-i", help = "Path to input file", nargs = Inf)
parser = add_argument(parser, "--output_plot", short = "-op", help = "Path to output plot", nargs = 1)
parser = add_argument(parser, "--output_table", short = "-ot", help = "Path to output table", nargs = 1)
argv = parse_args(parser)


#argv = list()
#argv$input = c("KG77_filter_counts.txt", "KG71_filter_counts.txt")
#argv$output = "summary_plot.pdf"

# Get input files and list over
dts = lapply(argv$input, function(x) {
  fread(x)
})
dt = rbindlist(dts)

# Add condition information
dt[, name := paste0(sample, " (", condition, ")")]

# Ordering
dt[, time_suffix := gsub("([0-9])", "", condition)]
dt[, time_prefix := as.numeric(gsub("([a-zA-Z:])", "", condition))]
setorder(dt, -time_suffix, time_prefix)


dt[, name := factor(name, levels = unique(name))]
dt[, filter := factor(filter, levels = c("input", 
                                         "hq_extracted",
                                         "prefix_match",
                                         "total_umis",
                                         "dedup_umis"))]
dt = dt[!is.na(filter)]

# Melt and plot
plt_cnt = ggplot(dt, aes(x = name, y = count, fill = filter)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(expand = c(0, 0), labels = scales::number_format()) +
  labs(y = "Number of reads", x = "Library ID", fill = "") +
  scale_fill_npg() +
  theme(strip.background = element_rect(fill = "#E8EDF0"),
        strip.text = element_text(size = 12, face = "bold"))

# Same plot for percentages
totals = dt[filter == "input", .(total = count), by = sample]
dt = merge(dt, totals)
dt[, fraction := count / total]

# Melt and plot
plt_perc = ggplot(dt, aes(x = name, y = fraction, fill = filter)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format()) +
  labs(y = "Percentage of input reads", x = "Library ID", fill = "") +
  scale_fill_npg() +
  theme(strip.background = element_rect(fill = "#E8EDF0"),
        strip.text = element_text(size = 12, face = "bold"))

# Combine
plt = (plt_cnt / plt_perc) + plot_layout(guides = "collect")


# Save plots
ggsave(argv$output_plot, plt, width = 14, height = 12)


# Write output table
dt_wide = merge(dcast(dt, sample + condition ~ filter, value.var = "count"),
                dcast(dt, sample + condition ~ filter, value.var = "fraction"),
                suffixes = c("", "_fraction"))
write.table(dt_wide, argv$output_table, col.names = TRUE, row.names = FALSE,  quote = FALSE, sep = "\t")