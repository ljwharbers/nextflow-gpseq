#!/usr/bin/env Rscript

## Author: Luuk Harbers
## Date: 2022-10-22
## Script to process allelic counts

## Load/install packages
packages = c("data.table", "argparser")
invisible(sapply(packages, function(x) suppressMessages(require(x, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))))

## Parse arguments
parser = arg_parser("Generate metadata table")
parser = add_argument(parser, "--runs", short = "-r", help = "List of run ids", nargs = Inf)
parser = add_argument(parser, "--conditions", short = "-c", help = "List of conditions", nargs = Inf)
parser = add_argument(parser, "--samples", short = "-s", help = "List of samples", nargs = Inf)
parser = add_argument(parser, "--filepaths", short = "-f", help = "List of bed filepaths", nargs = Inf)
parser = add_argument(parser, "--output", short = "-o", help = "Path to output file", nargs = 1)
argv = parse_args(parser)

# Check for equal lengths
lengths = sapply(list(argv$runs, argv$conditions, argv$samples, argv$filepaths), function(x) {
  length(x)
})
if(length(unique(lengths)) > 1) stop("Input argument are not of same length.\nPlease check that runs, conditions, samples, and filepaths are of identical lengths")

# Generate metadata table
dt = data.table(exid = argv$runs,
                cond = argv$conditions,
                libid = argv$samples,
                fpath = argv$filepaths)

# Write output
write.table(dt, argv$output, quote = F, row.names = F, col.names = T, sep = "\t")
