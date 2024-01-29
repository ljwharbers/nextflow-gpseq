#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2022-10-22
## Script to plot GPSeq Pizza plot

# Function for saving
# Function to save a plot in different file formats and show it within the notebook
save_and_plot <- function(x, bname, width, height, dpi = 300, 
                          use.cowplot = FALSE, ncol = 1, nrow = 1, base_height = 4, base_width = NULL, 
                          base_aspect_ratio = 1, plot = FALSE, onefile = TRUE, pointsize = 8,
                          output = c("cairo_ps", "cairo_pdf", "postscript", "png")) {
  
  # Create dir if not exists
    if (!dir.exists(dirname(bname))) {
    dir.create(dirname(bname), recursive = TRUE)
    }

    if ("png" %in% output) {
        save_plot(x, filename = file.path(paste0(bname, ".png")), ncol = 
                    ncol, nrow = nrow, base_height = base_height, base_width = base_width, 
                base_aspect_ratio = base_aspect_ratio, dpi = dpi)
    }
    if ("cairo_ps" %in% output) {
        while (!is.null(dev.list()))  invisible(dev.off())
        save_plot(x, filename = file.path(paste0(bname, "_cairo_ps.eps")), 
                ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                    base_width, base_aspect_ratio = base_aspect_ratio, device = cairo_ps)
    }
    if ("cairo_pdf" %in% output) {
        while (!is.null(dev.list()))  invisible(dev.off())
        save_plot(x, filename = file.path(paste0(bname, "_cairo_pdf.pdf")), 
                ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device = cairo_pdf)
    }
    if ("postscript" %in% output) {
        while (!is.null(dev.list()))  invisible(dev.off())
        save_plot(x, filename = file.path(paste0(bname, "_postscript.eps")), 
                ncol = ncol, nrow = nrow, base_height = base_height, base_width = 
                base_width, base_aspect_ratio = base_aspect_ratio, device = "ps")
    }
    while (!is.null(dev.list()))  invisible(dev.off())
 
    if (plot)
    print(x)
    dev.off()
}

## Load/install packages
packages = c("data.table", "argparser", "ggplot2", "dplyr", "paletteer", "cowplot")
invisible(sapply(packages, function(x) {
    require(x, character.only = TRUE)
    }))

## Parse arguments
parser = arg_parser("Generate summary table and plots")
parser = add_argument(parser, "--input", short = "-i", help = "Path to input rds", nargs = Inf)
parser = add_argument(parser, "--output", short = "-o", help = "Path to output plot folder", nargs = 1)
argv = parse_args(parser)

# Get list of all binsizes in the rds file
rds = readRDS(argv$input)
binsizes = rds[[1]]

# Remove chrom-wide
binsizes = binsizes[names(binsizes) != "chrom:wide"]

# Loop through binsizes and plot pizza plot
plot_list = lapply(binsizes, function(dt) {

    dt = dt[chrom %in% paste0("chr", c(1:22, "X"))]
    dt = dt[score >= 1 & score <= 2]

    # Get bins
    dt[, bin := seq_along(chrom)]
    dt[, end_cum := cumsum((end - start) + 1)]
    dt[, start_cum := c(1, end_cum[seq_along(end_cum) - 1] + 1)]

    # Randomly reorder bins
    set.seed(123)
    dt[, bin := sample(bin, size = .N), by = chrom]

    # Make chr_bounds
    chr_bounds = dt[, list(min = min(bin), max = max(bin), chrlen_bp = sum(end - start)), by = chrom]
    chr_bounds = chr_bounds %>% 
    mutate(mid = round(min + (max - min) / 2, 0),
            end_bp = cumsum(as.numeric(chrlen_bp)), 
            start_bp = end_bp - chrlen_bp, 
            mid_bp = round((chrlen_bp / 2) + start_bp, 0))

    # Plot
    plt = ggplot(dt, aes(x = bin, y = score, color = score)) +
                 geom_point(size = 0.5, alpha = 0.5) +
                 scale_color_viridis_c(option = "H", direction = 1) +
                 scale_y_reverse() +
                 geom_vline(data = chr_bounds, aes(xintercept = max), linewidth = 0.25, colour = "grey", linetype = "dashed") +
                 geom_text(data = chr_bounds, aes(x = mid, y = -Inf, label = chrom), size = 5, inherit.aes = F) +
                 coord_polar() +
                 labs(y = "GPSeq score", x = "", colour = "Log2 GPSeq Score") +
                 theme_void() + 
                 theme(legend.position = "right",
                 legend.direction = "vertical",
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank()) 
    
    ggsave(paste0(argv$output, "/pizza_plot_", dt[1, tag], ".png"), plt,
                  width = 12, height = 12)
})
