#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# License: MIT - Copyright (c) 2020 Gabriele Girelli
# ------------------------------------------------------------------------------

# UTILITIES ====================================================================

version = "v0.0.9"
if ("--version" %in% commandArgs(trailingOnly=TRUE)) {
    cat(sprintf("GPSeq-RadiCal %s\n\n", version))
    quit()
}
if ("--debug-info" %in% commandArgs(trailingOnly=TRUE)) {
    pd = as.data.frame(installed.packages())[,
        c("Package", "Version", "Priority")]
    rownames(pd) = NULL
    cat(sprintf("GPSeq-RadiCal %s\n", version))
    print(R.version)
    print(pd)
    quit()
}

# DEPENDENCIES =================================================================

depends_on = c("argparser", "data.table", "logging",
    "outliers", "pbapply", "rtracklayer")
tmp = lapply(depends_on, function(x) {
    expr = x %in% names(installed.packages()[, "Package"])
    if (!expr) {
        cat(sprintf(" Missing '%s' package.\n", x),
            sprintf("Install it with 'install.packages(\"%s\")'.\n", x))
    }
    stopifnot(expr)
})
pbapply::pboptions(type="timer")

# FUNCTIONS ====================================================================

chrom_to_chrom_id = function(chrom, nchrom=24, hetero=c("X", "Y")) {
    # Convert a chromosome name (e.g., "chr1", "chrX") to a numerical ID.
    if (grepl(":", chrom)) {
        return(floor(as.numeric(gsub(":", ".",
            substr(chrom, 4, nchar(chrom))))))
    } else {
        chrom = unlist(strsplit(chrom, "_", fixed=T))[1]
        chrom_id = substr(chrom, 4, nchar(chrom))
        if (chrom_id %in% hetero)
            chrom_id = nchrom - which(rev(hetero) == chrom_id) + 1
        chrom_id = as.numeric(chrom_id)
        stopifnot(!is.na(chrom_id))
        return(chrom_id)
    }
}

add_chrom_id = function(data, key = "chrom", nchrom=24, hetero=c("X", "Y")) {
    # Add chromosome ID to a data.table.
    # key should be the name of the column with chromosome names.
    stopifnot(data.table::is.data.table(data))
    stopifnot("chrom" %in% colnames(data))

    cid_table = data.table::data.table(
        chrom = as.character(unique(data$chrom)))
    cid_table$chrom_id = unlist(lapply(
        cid_table$chrom, FUN=chrom_to_chrom_id, nchrom, hetero))
    data.table::setkeyv(cid_table, "chrom")
    data.table::setkeyv(data, key)

    return(data[cid_table, , nomatch = 0])
}

assert = function(expr, msg) {
    if (!expr) logging::logerror(msg)
    stopifnot(expr)
}

btag2extremes = function(btag) {
    if (grepl(":", btag)) {
        extremes = as.numeric(unlist(strsplit(btag, ":")))
    } else {
        extremes = rep(as.numeric(btag), 2)
    }
    extremes = data.table::as.data.table(t(extremes))
    colnames(extremes) = c("size", "step")
    extremes$tag = btag
    return(extremes)
}

bstring2specs = function(bstring) {
    taglist = unlist(strsplit(bstring, ","))
    binspecs = data.table::rbindlist(lapply(taglist, btag2extremes))
    return(binspecs)
}

scientific_with_signif_digits = function(x) {
    x_string = sprintf("%d", x)
    integer_part = ifelse("."  %in% x_string,
        strsplit(x_string, ".", fixed=T)[[1]], x_string)
    for (i in seq(nchar(integer_part))) {
        base = as.numeric(substr(integer_part, 1, i))
        exponent = nchar(integer_part) - i
        if (base * (10**exponent) == x) {
            return(sprintf(paste0("%.", i-1, "e"), x))
        }
    }
    return(as.character(x))
}

mk_genome_wide_bins = function(brid, bspecs, cinfo, args) {
    bins = cinfo[, .(start=seq(start, end, by=bspecs[brid, step]),
        size=end), by=chrom]
    bins[, end := start + bspecs[brid, size] - 1]
    if (args$elongate_ter_bin) {
        bins = data.table::rbindlist(list(bins[end<=size, .(chrom, start, end)],
            bins[end > size, .(start=min(start), end=min(end)), by=chrom]))
    } else {
        bins = data.table::rbindlist(list(bins[end<=size, .(chrom, start, end)],
            bins[end > size, .(start=min(start), end=size[1]), by=chrom]))
    }
    nchrom = as.numeric(strsplit(args$chrom_tag, ":")[[1]][1])
    hetero = unlist(strsplit(strsplit(args$chrom_tag, ":")[[1]][2], ","))
    bins = add_chrom_id(bins, "chrom", nchrom, hetero)
    bins[, chrom := reorder(chrom, chrom_id)]
    bins[, chrom_id := NULL]
    bins = bins[order(chrom, start)]
    bins$tag = bspecs[brid, paste0(
        scientific_with_signif_digits(size), ":",
        scientific_with_signif_digits(step))]
    bins[, start := floor(start)]
    bins[, end := ceiling(end)]
    return(bins)
}

mk_roi_centered_bins = function(brid, bspecs, rois) {
    bins = data.table::copy(rois)
    bins[, start := (start+end)/2]
    bins[, end := start]
    half_width = ceiling(bspecs[brid, size]/2)
    bins[, start := start - half_width]
    bins[, end := end + half_width]
    bins$tag = bspecs[brid,
        paste0(scientific_with_signif_digits(size), ":rois")]
    return(bins)
}

import_gpseq_bed = function(brid, bmeta) {
    assert("fpath" %in% colnames(bmeta),
        "Missing 'fpath' column in metadata file.")
    brmeta = data.table::copy(bmeta[brid])
    o = data.table::as.data.table(rtracklayer::import.bed(brmeta$fpath))
    data.table::setnames(o, "seqnames", "chrom")
    o[, c("width", "strand", "name") := NULL]
    brmeta[, fpath := NULL]
    o[, condition := sprintf("cid_%d", brid)]
    return(o)
}

parse_bed_meta = function(bbmeta, args) {
    logging::loginfo(sprintf("Reading %d bed files.", nrow(bbmeta)))
    bd = data.table::rbindlist(pbapply::pblapply(seq_len(nrow(bbmeta)),
        import_gpseq_bed, bbmeta, cl=args$threads))
    bd[, end := start]
    logging::loginfo(sprintf("Dcasting bed data."))
    dups = bd[, .N, by=c("chrom", "start", "end", "condition")][N != 1]
    if (0 != nrow(dups)) {
        logging::logwarn(sprintf(
            "Found %d duplicated regions, using minimum signal.", nrow(dups)))
        print(dups)
    }
    bd = data.table::dcast(bd, chrom+start+end~condition,
        value.var="score", fill=0, fun.aggreg=min)
    if (3 <= args$export_level) {
        logging::loginfo("Exporting dcasted input bed...")
        saveRDS(bd, file.path(args$exp_output_folder, "input_bed.rds"))
    }
    return(bd)
}

otag2specs = function(otag) {
    outlier_specs = unlist(strsplit(otag, ":"))
    outlier_specs = outlier_specs[1:min(2, length(outlier_specs))]
    outlier_specs = data.table::as.data.table(t(outlier_specs))
    colnames(outlier_specs) = c("method", "threshold")
    outlier_specs[, method := tolower(method)]
    outlier_specs[, threshold := as.numeric(threshold)]
    assert(outlier_specs$method %in% outlier_methods,
        sprintf("Unrecognized outlier method '%s'.", outlier_specs$method))
    outlier_specs$alpha = .5  # Default
    outlier_specs$lim = 1.5   # Default
    if (outlier_specs$method == "iqr") {
        outlier_specs[, lim := threshold]
    } else {
        outlier_specs[, alpha := threshold]
    }
    outlier_specs[, threshold := NULL]
    return(outlier_specs)
}

export_output = function(odata, odir, format, suffix, rm_tag=FALSE) {
    assert("tag" %in% colnames(odata), "Cannot find required 'tag' column.")
    assert(format %in% c("tsv", "tsv.gz", "csv", "csv.gz", "rds"),
        sprintf("Unrecognized output format '%s'.", format))
    opath_base = sprintf("%s.bins_%s", suffix, gsub(":", "_", odata[1, tag]))
    columns_to_export = colnames(odata)
    if (rm_tag) {
        columns_to_export = columns_to_export["tag" != columns_to_export]
    }
    if ("rds" == format) {
        saveRDS(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.rds", opath_base)))
    }
    if (format %in% c("tsv", "tsv.gz")) {
        data.table::fwrite(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.%s", opath_base, format)),
            sep="\t", na="NA", quote=F)
    }
    if (format %in% c("csv", "csv.gz")) {
        data.table::fwrite(odata[, .SD, .SDcols=columns_to_export],
            file.path(odir, sprintf("%s.%s", opath_base, format)),
            sep=",", na="NA", quote=F)
    }
}

calc_condition_outliers_stats = function(x, specs) {
    x = x[0 != x]
    is_outlier = outliers::scores(x,
        type=specs$method, prob=1-specs$alpha, lim=specs$lim)
    ostats = as.numeric(summary(x[is_outlier]))
    ostats = c(ostats, sum(is_outlier), length(x))
    return(ostats)
}

rm_condition_outliers = function(x, specs) {
    is_outlier = outliers::scores(x[0 != x],
        type=specs$method, prob=1-specs$alpha, lim=specs$lim)
    x[which(0 != x)[is_outlier]] = 0
    return(x)
}

calc_bed_outliers_stats = function(bd, args) {
    logging::loginfo(sprintf(
        "Calculating outlier stats. [%s]", args$bed_outlier_tag))
    outlier_stats = bd[, lapply(.SD,
        calc_condition_outliers_stats, otag2specs(args$bed_outlier_tag)),
        .SDcols=args$cond_cols]
    outlier_stats = data.table::data.table(t(outlier_stats))
    colnames(outlier_stats) = c(names(summary(1)), "nout", "ntot")
    outlier_stats_opath = file.path(args$exp_output_folder, "outlier_stats.tsv")
    logging::loginfo(sprintf(
        "Exporting outlier stats to '%s'.", outlier_stats_opath))
    data.table::fwrite(outlier_stats, outlier_stats_opath, sep="\t")
}

rm_bed_outliers = function(bd, args) {
    logging::loginfo(sprintf(
        "Removing outliers. [%s]", args$bed_outlier_tag))
    bd[, c(args$cond_cols) := lapply(.SD,
        rm_condition_outliers, otag2specs(args$bed_outlier_tag)),
        .SDcols=args$cond_cols]
    bd = bd[0 != apply(bd[, .SD, .SDcols=args$cond_cols], MARGIN=1, FUN=sum)]
    if (3 <= args$export_level) {
        logging::loginfo(
            "Exporting dcasted input bed after outlier removal...")
        saveRDS(bd, file.path(args$exp_output_folder, "clean_bed.rds"))
    }
    return(bd)
}

bin_bed_data = function(bbins, bd, args, site_universe=NULL) {
    logging::loginfo(sprintf("Binning... [%s]", bbins[1, tag]))
    binned = data.table::rbindlist(pbapply::pblapply(split(bd, bd$chrom),
        bin_chromosome, bbins, args, site_universe
        ))[order(tag, chrom, start, cid)]
    binned[, start := floor(start)]
    binned[, end := ceiling(end)]
    return(binned)
}

bin_chromosome = function(
    bbd, bbins, args, site_universe=NULL) {
    selected_chromosome = bbd[1, chrom]
    bbins2 = data.table::copy(bbins)[selected_chromosome==chrom]
    data.table::setkeyv(bbins2, bed3_colnames)

    ovlps = data.table::foverlaps(bbins2, bbd)[!is.na(tag)]
    ovlps[, c(bed3_colnames[2:3]) := .(NULL, NULL)]
    data.table::setnames(ovlps,
        paste0("i.", bed3_colnames[2:3]), bed3_colnames[2:3])

    nreads = ovlps[, lapply(.SD, sum, na.rm=T), by=c(bed3_colnames, "tag"),
        .SDcols=args$cond_cols][order(tag, chrom, start)]
    nreads = data.table::melt(nreads, id.vars=c(bed3_colnames, "tag"))
    data.table::setnames(nreads, c("variable", "value"), c("cid", "nreads"))
    nreads[, cid := match(cid, args$cond_cols)]
    data.table::setkeyv(nreads, c(bed3_colnames, "tag", "cid"))

    if ("universe" == args$site_domain) {
        assert(!is.null(site_universe),
            "Missing site universe data with site domain 'universe'.")
        nsites = data.table::foverlaps(
            site_universe, bbins2[chrom==selected_chromosome]
            )[!is.na(start), .(
                tag=bbins2[1, tag], cid=seq_along(args$cond_cols), nsites=.N
            ), by=bed3_colnames]
    } else {
        if ("union" == args$site_domain) {
            nsites = ovlps[, lapply(.SD, function(x) length(x)),
                by=c(bed3_colnames, "tag"), .SDcols=args$cond_cols
                ][order(tag, chrom, start)]
        } else {
            nsites = ovlps[, lapply(.SD, function(x) sum(0 != x)),
                by=c(bed3_colnames, "tag"), .SDcols=args$cond_cols
                ][order(tag, chrom, start)]
        }
        nsites = data.table::melt(nsites, id.vars=c(bed3_colnames, "tag"))
        data.table::setnames(nsites, c("variable", "value"), c("cid", "nsites"))
        nsites[, cid := match(cid, args$cond_cols)]
    }
    nsites[is.na(nsites), nsites := 0]
    data.table::setkeyv(nsites, c(bed3_colnames, "tag", "cid"))

    combined = nreads[nsites]
    combined[, lib_nreads := as.numeric(args$total_lib_nreads)[cid]]
    combined[, chrom_nreads := as.numeric(args$total_chrom_nreads[
        selected_chromosome==chrom, .SD, .SDcols=args$cond_cols])[cid]]

    return(combined)
}

export_binned_bed_data = function(binned, odir, format="rds") {
    export_output(binned, odir, format, "binned")
}

export_masked_data = function(masked, odir, format="rds") {
    export_output(masked, odir, format, "binned_masked")
}

estimate_centrality = function(bbind, normalize_by) {
    bbind = bbind[,
        .(pres=nreads/nsites/get(sprintf("%s_nreads", normalize_by))),
        by=c(bed3_colnames, "tag", "cid")]
    centr = bbind[order(chrom, start, cid), .(
            score=sum(pres[2:.N]/pres[1:(.N-1)])
        ), by=c(bed3_colnames, "tag")]
    return(centr)
}

export_estimated_centrality = function(centr, odir, format="rds") {
    export_output(centr, odir, format, "estimated")
}

calc_rescale_factor_by_lib = function(estmd, specs) {
    estmd[!is.na(score),
        outlier_status := outliers::scores(estmd[!is.na(score), score],
            type=specs$method, prob=1-specs$alpha, lim=specs$lim)]
    lowest = estmd[FALSE == outlier_status, min(score, na.rm = T)]
    highest = estmd[FALSE == outlier_status, max(score-lowest, na.rm = T)]
    estmd[, outlier_status := NULL]
    return(c(lowest, highest))
}

rescale_by_lib = function(estmd, specs) {
    rescaling_factors = calc_rescale_factor_by_lib(estmd, specs)
    estmd[, score := .(score - rescaling_factors[1])]
    estmd[, score := .(score / rescaling_factors[2])]
    estmd[, score := 2**score]
    return(estmd)
}

rescale_by_chr = function(estmd, specs) {
    data.table::rbindlist(by(estmd, estmd$chrom, rescale_by_lib, specs))
}

export_rescaled_centrality = function(rscld, odir, format="tsv.gz") {
    export_output(rscld, odir, format, "rescaled", rm_tag=TRUE)
}

apply_intersection_site_domain = function(bd, args) {
    logging::loginfo("Retaining only sites shared across conditions...")
    n_condition_empty = rowSums(0 == bd[, .SD, .SDcols=args$cond_cols])
    bd = bd[0 == n_condition_empty]
    logging::loginfo(sprintf("Retained %d/%d (%.1f%%) common sites.",
        nrow(bd), length(n_condition_empty),
        nrow(bd)/length(n_condition_empty)*100))
    if (3 <= args$export_level) {
        logging::loginfo(
            "Exporting dcasted input bed after site intersection...")
        saveRDS(bd, file.path(args$exp_output_folder,
            "clean_bed.intersected.rds"))
    }
    return(bd)
}

mask_dcasted_bed = function(bd, args) {
    assert(file.exists(args$mask_bed),
        sprintf("Cannot find mask bed file '%s'.", args$mask_bed))
    mask = data.table::as.data.table(
        rtracklayer::import.bed(args$mask_bed))[, .(chrom=seqnames, start, end)]
    assert(all(mask[, start <= end]), sprintf(
        "Mask not conforming to end >= start condition on row: %d",
        which(mask[, start > end])))
    data.table::setkeyv(mask, bed3_colnames)
    data.table::setkeyv(bd, bed3_colnames)
    masked = data.table::foverlaps(bd, mask)
    masked[!is.na(start), c(args$cond_cols) := lapply(.SD, function(x) {
            x = 0
            return(x)
        }), .SDcols=args$cond_cols]
    masked = masked[, bed3_colnames[2:3] := .(NULL, NULL)]
    data.table::setnames(masked,
        paste0("i.", bed3_colnames[2:3]), bed3_colnames[2:3])
    if (3 <= args$export_level) {
        logging::loginfo("Exporting masked dcasted input bed...")
        saveRDS(masked, file.path(args$exp_output_folder, "masked_bed.rds"))
    }
    return(masked)
}

mask_binned_track = function(bbins, mask) {
    if ("chrom:wide" == bbins[1, tag]) {
        logging::logwarn("Skipped masking for chromosome-wide bins.")
        return(bbins)
    }
    data.table::setkeyv(bbins, bed3_colnames)
    masked = unique(data.table::foverlaps(bbins, mask)[,
        .(tag, nreads, nsites, lib_nreads, chrom_nreads,
            mask_overlaps=.N, mask_overlapped=!is.na(start)),
        by=c(bed3_colnames[1], paste0("i.", bed3_colnames[2:3]),  "cid")])
    data.table::setnames(masked,
        paste0("i.", bed3_colnames[2:3]), bed3_colnames[2:3])
    masked[!(mask_overlapped), mask_overlaps := 0]
    masked[, mask_overlapped := NULL]
    masked[0 < mask_overlaps, c("nreads", "nsites") := .(0, 0)]
    masked[, mask_overlaps := NULL]

    return(masked)
}

mask_binned = function(binned, args) {
    assert(file.exists(args$mask_bed),
        sprintf("Cannot find mask bed file '%s'.", args$mask_bed))
    mask = data.table::as.data.table(
        rtracklayer::import.bed(args$mask_bed))[, .(chrom=seqnames, start, end)]
    assert(all(mask[, start <= end]), sprintf(
        "Mask not conforming to end >= start condition on row: %d",
        which(mask[, start > end])))
    data.table::setkeyv(mask, bed3_colnames)
    binned = pbapply::pblapply(binned, mask_binned_track, mask, cl=args$threads)
    if (1 <= args$export_level) {
        logging::loginfo(sprintf("Exporting binned bed data..."))
        tmp = lapply(binned, export_masked_data, args$exp_output_folder)
    }
    return(binned)
}

rescale_estimated = function(estimated, args) {
    logging::loginfo(sprintf("Rescaling estimates... [%s]", args$normalize_by))
    if ("chrom" == args$normalize_by) {
        if (args$chromosome_wide) logging::logwarn(
            "Skipped rescaling by chromosome for chromosome-wide bins.")
        rescaled = pbapply::pblapply(estimated, function(estmd) {
            if ("chrom:wide" == estmd[1, tag]) return(estmd)
            estmd = rescale_by_chr(estmd, otag2specs(args$score_outlier_tag))
        }, cl=args$threads)
    }
    if ("lib" == args$normalize_by) {
        rescaled = pbapply::pblapply(estimated, function(estmd) {
            estmd = rescale_by_lib(estmd, otag2specs(args$score_outlier_tag))
        }, cl=args$threads)
    }
    logging::loginfo(sprintf("Exporting rescaled centrality..."))
    tmp = lapply(rescaled, export_rescaled_centrality, args$exp_output_folder)
    return(rescaled)
}

process_experiment = function(bbmeta, bins, args) {
    exid = bbmeta[1, exid]
    logging::loginfo(sprintf("Processing experiment '%s'.", exid))
    args$exp_output_folder = file.path(args$output_folder)
    #dir.create(args$exp_output_folder)

    assert(2 <= nrow(bbmeta),
        sprintf("Provide at least two bed files. [%s]", exid))

    logging::loginfo("Storing metadata.")
    data.table::fwrite(bbmeta,
        file.path(args$exp_output_folder, "bed.metadata.tsv"), sep="\t")

    args$cond_cols = sprintf("cid_%d", seq_len(nrow(bbmeta)))

    # Read bed files -----------------------------------------------------------

        bd = parse_bed_meta(bbmeta, args)

    # Get outlier stats --------------------------------------------------------

        if (0 != nchar(args$bed_outlier_tag)) {
            calc_bed_outliers_stats(bd, args)
        } else {
            logging::loginfo(sprintf("Skipped outlier removal."))
        }

    # Clean bed outliers -------------------------------------------------------

        if (0 != nchar(args$bed_outlier_tag)) bd = rm_bed_outliers(bd, args)

    # Apply intersection site domain -------------------------------------------

        if ("intersection" == args$site_domain) {
            bd = apply_intersection_site_domain(bd, args)
        } else if ("universe" == args$site_domain) {
            logging::loginfo(sprintf(
                "Reading site bed file '%s'", args$site_bed))
            site_universe = data.table::as.data.table(
                rtracklayer::import.bed(args$site_bed))[,
                    .(chrom=seqnames, start, end=start)]
        }

    # Calculate normalization factors ------------------------------------------

        logging::loginfo(sprintf("Calculating normalization factors."))
        args$total_lib_nreads = bd[, lapply(.SD, sum), .SDcols=args$cond_cols]
        args$total_chrom_nreads = bd[, lapply(.SD, sum),
            by=chrom, .SDcols=args$cond_cols]

    # Masking bed --------------------------------------------------------------

        if (!is.na(args$mask_bed)) bd = mask_dcasted_bed(bd, args)

    # Assign to bins -----------------------------------------------------------

        bin_tags = bins[, unique(tag)]
        binned = by(bins, bins$tag, bin_bed_data, bd, args, site_universe)
        if (1 <= args$export_level) {
            logging::loginfo(sprintf("Exporting binned bed data..."))
            tmp = lapply(binned, export_binned_bed_data, args$exp_output_folder)
        }

    # Calculate centrality -----------------------------------------------------

        logging::loginfo(sprintf("Estimating centrality..."))
        estimated = pbapply::pblapply(
            binned, estimate_centrality, args$normalize_by, cl=args$threads)
        if (1 <= args$export_level) {
            logging::loginfo(sprintf("Exporting estimated centrality..."))
            tmp = lapply(estimated,
                export_estimated_centrality, args$exp_output_folder)
        }

    # Rescale estimates --------------------------------------------------------

        if (0 == nchar(args$score_outlier_tag)) {
            logging::loginfo(sprintf("Skipped rescaling."))
            logging::loginfo(sprintf("Exporting estimated centrality..."))
            tmp = lapply(estimated, export_estimated_centrality,
                args$exp_output_folder, format="tsv.gz")
            saveRDS(estimated, file.path(args$exp_output_folder,
                "gpseq-radical.out.rds"))
            return(estimated)
        } else {
            rescaled = rescale_estimated(estimated, args)
            saveRDS(rescaled, file.path(args$exp_output_folder,
                "gpseq-radical.out.rds"))
            return(rescaled)
        }
}

# COMMON PARAMETERS ============================================================

bed3_colnames = c("chrom", "start", "end")
outlier_methods = c("z", "t", "chisq", "iqr", "mad")

# INPUT ========================================================================

parser = argparser::arg_parser("
Provide the path to a 4-columns tabulation-separated metadata file, containing
the sequencing run ID, a condition description, the library ID, and full
absolute path to each BED file (possibly gzipped) from the GPSeq experiment.
The bed files should be reported in order of condition strength, with the top
rows being *weaker* than bottom ones. Also, the first line should contain the
column headers: exid, cond, libid, and path. An example metadata file would be
is available at
https://github.com/ggirelli/GPSeq-RadiCal/blob/master/example_meta.tsv.

Fore more details, use the '--more-help' option.
", name="gpseq-radical.R")

parser = argparser::add_argument(parser, arg="bmeta_path",
    help="Path to bed metadata tsv file.")
parser = argparser::add_argument(parser, arg="output_folder",
    help="Path to output directory. Stops if it already exists.")

parser = argparser::add_argument(parser, arg="--cinfo-path",
    help="Path to bed file with chromsome sizes. Queries UCSC if not provided.",
    default=NULL)
parser = argparser::add_argument(parser, arg="--ref-genome", help=paste0(
        "Used when --cinfo-path is not provided to query UCSC. ",
        "If the provided reference genome is not found, reverts to 'hg38'."))

parser = argparser::add_argument(parser, arg="--bin-tags",
    help="Comma-separated bin tags. Use --more-help for more details.",
    default="1e6:1e5,1e5:1e4")
parser = argparser::add_argument(parser, arg="--bin-bed",
    help="Path to bed with regions on which to build bins.")

parser = argparser::add_argument(parser, arg="--bed-outlier-tag",
    help="Method:threshold for input bed outlier removal.",
    default="chisq:0.01")
parser = argparser::add_argument(parser, arg="--score-outlier-tag",
    help="Method:threshold for output score outlier removal (rescaling).",
    default="iqr:1.5")

parser = argparser::add_argument(parser, arg="--normalize-by", help=paste0(
        "Whether rescaling should be performed library-wise ('lib') ",
        "or chromosome-wise ('chrom')."),
    default="lib")

parser = argparser::add_argument(parser, arg="--site-domain",
    help="Site domain method. Use --more-help for more detail.",
    default="universe")
parser = argparser::add_argument(parser, arg="--site-bed", help=paste0(
        "Path to bed file with recognition site locations for the enzyme ",
        "in use. This is required when using '--site-domain universe'."),
    default=NULL)

parser = argparser::add_argument(parser, arg="--mask-bed",
    help="Path to bed file with regions to be masked.",
    default=NULL)

parser = argparser::add_argument(parser, arg="--chrom-tag",
    help="Two values, column (:) separated: the number of chromosomes, and a
    string with comma-separated heterosome names.",
    default="24:X,Y")

parser = argparser::add_argument(parser, arg="--threads", help=paste0(
        "Number of threads for parallelization. ",
        "Optimal when using at least one core per bed file."),
    default=1)
parser = argparser::add_argument(parser, arg="--export-level",
    help="Limits the amount of output. Use --more-help for more details",
    default=0)

parser = argparser::add_argument(parser, arg="--chromosome-base-delim",
    help="Delimi of chromosome base in chromosome name for patch recognition.",
    default="_")
parser = argparser::add_argument(parser,
    arg="--chromosome-strict-match", flag=TRUE,
    help="Use only chromosomes matching the expected names, discard patches.")

parser = argparser::add_argument(parser, arg="--chromosome-wide", flag=TRUE,
    help="Use this option to calculate also on chromosome-wide bins.")
parser = argparser::add_argument(parser, arg="--elongate-ter-bin", flag=TRUE,
    help=paste0("Use this option to elongate chromosome-terminal bins ",
        "and have equally-sized bins."))

parser = argparser::add_argument(parser, arg="--more-help", flag=TRUE,
    help="Show extended help page and exit")
parser = argparser::add_argument(parser, arg="--version", flag=TRUE,
    help="Show script version and exit")
parser = argparser::add_argument(parser, arg="--debug-info", flag=TRUE,
    help="Show debugging info and exit")

if ("--more-help" %in% commandArgs(trailingOnly=TRUE)) {
    cat("
Provide the path to a 4-columns tabulation-separated metadata file, containing
the sequencing run ID, a condition description, the library ID, and full
absolute path to each BED file (possibly gzipped) from the GPSeq experiment.
The bed files should be reported in order of condition strength, with the top
rows being *weaker* than bottom ones. Also, the first line should contain the
column headers: exid, cond, libid, and path. An example metadata file would be
is available at
https://github.com/ggirelli/GPSeq-RadiCal/blob/master/example_meta.tsv.

Binning can be specified by providing comma-separated bin 'tags', in the format
of 'bin_size:bin_step' or 'bin_size'. If 'bin_step' is not specified, the bins
are generated in a non-overlapping manner. Use the '--chromosome-wide' option to
estimated centrality also on chromosome-wide bins. Use the '--elongate-ter-bin'
option to elongate chromosome terminal bins to the specified bin size, the
default behaviour is for terminal bins' end to coincide with the chromosome end.

Using the '--bin-bed' option. it is also possible to provide a bed file with
regions of interest on which to build the bins. The bins are built with the bin
sizes specified in '--bin-tag' and centered on the midpoint of the regions of
interest. In other words, the bin step information is disregarded entirely. This
can be used, for example, to estimate radiality of FISH probe regions. When
using '--bin-bed', options '--cinfo-path', '--ref-genome', '--chromosome-wide',
and '--elongate-ter-bin' are ignored.

Outlier removal method can be specified with the 'method:threshold' format.
Available methods: z, t, chisq, iqr, and mad. The specified threshold is an
alpha threshold on the outlier score p-value for the z, t, chisq, and mad
methods, and should thus be 0 < threshold < 1. In the case of the iqr method,
outliers are identified by comparing their distance to the closest quartile (Q1
or Q3) and the product threshold*iqr. More details are available here:
https://www.rdocumentation.org/packages/outliers/versions/0.14/topics/scores

Outliers can be removed from input bed files by using the '--bed-outlier-tag'
option. To skip outlier removal use '--bed-outlier-tag \"\"'.

Moreover, centrality outliers are detected and used to rescale the final
estimates when using the '--score-outlier-tag' option. To skip score rescaling
use '--score-outlier-tag \"\"'. Rescaling can be performed in a
chromosome-wise or library-wise manner, by using the '--normalize-by chr' or
'--normalize-by lib' options. The '--normalize-by' is also used for centrality
calculation, even if rescaling is skipped.

The '--site-domain' option can be used to specify which recognition sites to
consider when estimating centrality. The 'separate' domain considers all
(non-outlier) sites with at least one mapped read in a condition. The 'union'
domain considers all (non-outlier) sites with at least one mapped read in at
least one condition. The 'intersection' domain retains only recognition sites
with at least on mapped read in all conditions. The default 'universe' domain considers
all known recognition site locations in the reference genome. When using
'--site-domain universe', a bed file with known recognition site locations must
be provided using the '--site-bed' option.

The output centrality estimation can be masked (masking performed after
binning) to removed known problematic regions. Use the '--mask-bed' option to
provide the path to a bed file with regions to be masked. A bin is masked if it
overlaps with any regions in the provided mask bed file.

The '--export-level' option is useful to limit the amount of output files
produced by the script. The default export level is set to 0.
 0: export only final results and statistics.
 1: export intermediate results.
 2: export supplementary (bins) files.
 3: export dcasted input (bed, clean bed, domained bed).

Finally, use the '--threads' option to specify how many threads to use for
parallelization. An optimal number of threads coincides with the number of
bed files in the input metadata file.
    \n")
    quit()
}

args = argparser::parse_args(parser)
assert(args$normalize_by %in% c("chrom", "lib"),
    sprintf("Unrecognized 'normalize_by' value: '%s'", args$normalize_by))
assert(args$site_domain %in% c("separate", "union", "intersection", "universe"),
    sprintf("Unrecognized 'site_domain' value: '%s'", args$site_domain))
if ("universal" == args$site_domain) {
    assert(!is.na(args$site_bed),
        "A 'site_bed' path must be specified when 'site_domain' is 'universe'.")
    assert(file.exists(args$site_bed),
        sprintf("Cannot find 'site_bed' file '%s'.", args$site_bed))
}

# RUN ==========================================================================

# Prelude ----------------------------------------------------------------------

    assert(!dir.exists(args$output_folder), paste0(
        sprintf("Output folder '%s' already exists. ", args$output_folder),
        "This can lead to mixed settings output and should be avoided. ",
        "Stopping."))
    dir.create(args$output_folder)

    logging::basicConfig()
    log_path = file.path(args$output_folder, "gpseq-radical.log")
    logging::loginfo(sprintf("This log will be stored at '%s'.", log_path))
    logging::addHandler(logging::writeToFile, file=log_path)

    logging::loginfo(sprintf("Created output folder '%s'.", args$output_folder))
    settings_path = file.path(args$output_folder, "gpseq-radical.opts.rds")
    saveRDS(args, settings_path)
    logging::loginfo(sprintf(
        "Exported input parameters to '%s'.", settings_path))

    data.table::setDTthreads(args$threads)
    if (args$threads != data.table::getDTthreads()) {
        logging::logwarn(sprintf("Changed from the requested %d to %d threads.",
            args$threads, data.table::getDTthreads()))
        args$threads = data.table::getDTthreads()
    }

# Parse metadata ---------------------------------------------------------------

    logging::loginfo(sprintf("Parsing metadata from '%s'.", args$bmeta_path))
    bmeta = data.table::fread(args$bmeta_path)
    assert(2 <= nrow(bmeta), "Provide at least two bed files.")
    logging::loginfo("Storing metadata.")
    data.table::fwrite(bmeta,
        file.path(args$output_folder, "bed.metadata.tsv"), sep="\t")

# Read chromosome info bed -----------------------------------------------------

    if (is.na(args$bin_bed)) {
        cinfo = NULL
        if (is.na(args$cinfo_path)) {
            logging::loginfo("Opening UCSC browser session...")
            ucsc = rtracklayer::browserSession("UCSC")
            rtracklayer::genome(ucsc) = args$ref_genome
            logging::loginfo(sprintf(
                "Querying UCSC for '%s' chromosome info...",
                rtracklayer::genome(ucsc)))
            cinfo = data.table::data.table(rtracklayer::getTable(
                rtracklayer::ucscTableQuery(ucsc,
                    table="chromInfo")))
            cinfo = cinfo[, .(start=1, end=size), by=chrom]
        } else {
            assert(file.exists(args$cinfo_path),
                sprintf("Cannot find chromosome info bed file '%s'.",
                    args$cinfo_path))
            logging::loginfo(sprintf(
                "Reading chromosome info from '%s'.", args$cinfo_path))
            cinfo = data.table::as.data.table(
                rtracklayer::import.bed(args$cinfo_path))
            data.table::setnames(cinfo, "seqnames", "chrom")
            cinfo[, c("width", "strand") := NULL]
            cinfo[, start := 1] # Added this because input chromsize format is different when calculating it from genome.fa.fai 
        }
        assert(!is.null(cinfo), "Failed to build or retrieve chromosome info.")
        assert(nrow(cinfo[start == 1]) == nrow(cinfo),
            "Chromosome start should be 1.")
    }

# Retain chromosomes according to chromosome tag -------------------------------

    if (!is.null(args$chrom_tag)) {
        if (exists("cinfo")) {        
        chrom_tag = unlist(strsplit(args$chrom_tag, ":"))
        chromosomes = paste0("chr", c(1:as.numeric(chrom_tag[1]),
            unlist(strsplit(chrom_tag[2], ","))))
        cinfo$chrom_base = unlist(lapply(as.character(cinfo$chrom),
            function(x) unlist(strsplit(x, args$chromosome_base_delim, fixed=T))[1]))
        if (args$chromosome_strict_match) {
            cinfo = cinfo[chrom %in% chromosomes]
        } else {
            cinfo = cinfo[chrom_base %in% chromosomes]
        }
        cinfo[, chrom_base := NULL]
    }
    }
 

# Build bins -------------------------------------------------------------------

    logging::loginfo(sprintf("Building bins."))
    bspecs = bstring2specs(args$bin_tags)
    if (0 == nrow(bspecs)) bins = data.table::data.table()
    if (!is.na(args$bin_bed)) {
        assert(file.exists(args$bin_bed), sprintf(
            "Cannot find bin bed file '%s'.", args$bin_bed))
        rois = data.table::as.data.table(rtracklayer::import.bed(args$bin_bed
            ))[, .(chrom=seqnames, start, end)]
        if (0 == nrow(rois)) {
            bins = data.table::data.table()
        } else if (0 < nrow(bspecs)) {
            bins = data.table::rbindlist(pbapply::pblapply(
                seq_len(nrow(bspecs)), mk_roi_centered_bins,
                bspecs, rois, cl=args$threads))
            bins = unique(bins)
        }
    } else {
        if (0 < nrow(bspecs)) {
            bins = data.table::rbindlist(pbapply::pblapply(
                seq_len(nrow(bspecs)), mk_genome_wide_bins,
                bspecs, cinfo, args, cl=args$threads))
        }
        if (args$chromosome_wide) {
            args$chromosome_wide_bins = data.table::copy(cinfo)
            args$chromosome_wide_bins[, tag := "chrom:wide"]
            bins = data.table::rbindlist(list(args$chromosome_wide_bins, bins))
        }
    }
    data.table::setkeyv(bins, bed3_colnames)
    assert(0 < nrow(bins), "No bins built. Stopping.")
    if (2 <= args$export_level) {
        logging::loginfo("Exporting bins...")
        saveRDS(bins, file.path(args$output_folder, "bins.rds"))
    }

# Process one experiment at a time ---------------------------------------------

    assert("exid" %in% colnames(bmeta), "Missing 'exid' column from metadata.")
    output = by(bmeta, bmeta$exid, process_experiment, bins, args)
    saveRDS(output, file.path(args$output_folder, "gpseq-radical.out.rds"))

# ------------------------------------------------------------------------------

logging::loginfo(sprintf("Done."))

################################################################################