#!/usr/bin/env Rscript
# cnv_heatmap.R
# Produce a CNV heatmap from a combined segmentation file.
# Rows = samples, x-axis = chromosome (genome coordinates), colour = seg.mean
# Optionally annotates rows with up to 3 metadata-driven group sidebars.
#
# Usage:
#   Rscript cnv_heatmap.R --seg combined_seg.tsv --output cnv_heatmap.pdf
#   Rscript cnv_heatmap.R --seg combined_seg.tsv --output cnv_heatmap.png \
#           --genome hg38 --cap 2 --no_sex_chr \
#           --metadata meta.csv --meta_id_col Sentrix_ID \
#           --group_col "Diagnosis,Gender,Batch" \
#           --group_palette "Set1,Set2,Dark2" \
#           --sort_by_group
#
# Expected segmentation columns (CBS / conumee style):
#   ID | chrom | loc.start | loc.end | num.mark | seg.mean
# Also accepts:
#   Sample | Chr | Start | End | Num_Probes | Seg_Mean  (auto-detected)
#
# --group_col    : comma-separated, up to 3 columns  e.g. "Diagnosis,Gender"
# --group_palette: comma-separated palettes, one per group col e.g. "Set1,Set2"
# --group_colors : semicolon-separated colour lists, one list per group col
#                  e.g. "#E41A1C,#377EB8;#4DAF4A,#984EA3"

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(grid)
})

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
option_list <- list(
  make_option(c("--seg"),           type = "character", help = "Combined segmentation file (TSV or CSV)"),
  make_option(c("--output"),        type = "character", default = "cnv_heatmap.pdf",
              help = "Output file (.pdf, .png, or .svg) [default: %default]"),
  make_option(c("--genome"),        type = "character", default = "hg19",
              help = "Reference genome: hg19 or hg38 [default: %default]"),
  make_option(c("--cap"),           type = "double",    default = 1.5,
              help = "Cap |seg.mean| at this value for colour scale [default: %default]"),
  make_option(c("--bin_kb"),        type = "integer",   default = 500L,
              help = "Bin size in kb for rasterising segments [default: %default]"),
  make_option(c("--sample_order"),  type = "character", default = NULL,
              help = "Optional file with sample IDs (one per line) to set row order"),
  make_option(c("--width"),         type = "double",    default = 20,
              help = "Plot width in inches [default: %default]"),
  make_option(c("--height"),        type = "double",    default = 0,
              help = "Plot height in inches; 0 = auto-scale [default: auto]"),
  make_option(c("--title"),         type = "character", default = "CNV Heatmap",
              help = "Plot title [default: %default]"),
  make_option(c("--dpi"),            type = "integer",   default = 300L,
              help = "DPI for png output [default: %default]"),
  make_option(c("--autoscale_h"),   action = "store_true", default = FALSE,
              help = "Auto-scale height to number of samples"),

  # --- Sex chromosome display ---
  make_option(c("--no_sex_chr"),    action = "store_true", default = FALSE,
              help = "Exclude chrX and chrY from the heatmap"),
  make_option(c("--no_chrX"),       action = "store_true", default = FALSE,
              help = "Exclude chrX only"),
  make_option(c("--no_chrY"),       action = "store_true", default = FALSE,
              help = "Exclude chrY only"),

  # --- Label options ---
  make_option(c("--label_col"),    type = "character", default = NULL,
              help = "Metadata column to use for y-axis sample labels (e.g. Sample_Name)"),
  make_option(c("--no_labels"),    action = "store_true", default = FALSE,
              help = "Hide all y-axis sample labels"),

  # --- CNV colour options ---
  make_option(c("--color_loss"),    type = "character", default = "#053061",
              help = "Colour for maximum loss [default: deep blue %default]"),
  make_option(c("--color_gain"),    type = "character", default = "#67001F",
              help = "Colour for maximum gain [default: deep red %default]"),
  make_option(c("--color_mid"),     type = "character", default = "#F7F7F7",
              help = "Colour for neutral (seg.mean = 0) [default: near-white %default]"),
  make_option(c("--color_na"),      type = "character", default = "#CCCCCC",
              help = "Colour for bins with no coverage [default: light grey %default]"),

  # --- Metadata / group annotation ---
  make_option(c("--metadata"),      type = "character", default = NULL,
              help = "CSV file with sample metadata for group annotation"),
  make_option(c("--group_col"),     type = "character", default = NULL,
              help = paste("Comma-separated metadata column(s) to colour samples by,",
                           "up to 3. e.g. \"Diagnosis,Gender,Batch\"",
                           "(required when --metadata is given)")),
  make_option(c("--meta_id_col"),   type = "character", default = NULL,
              help = paste("Column in --metadata whose values match sample IDs in",
                           "the segmentation file (e.g. Sentrix_ID).",
                           "Auto-detected if omitted.")),
  make_option(c("--sort_by_group"), action = "store_true", default = FALSE,
              help = "Sort samples by first group column, then alphabetically within group"),
  make_option(c("--group_palette"), type = "character", default = "Set1",
              help = paste("Comma-separated RColorBrewer palette(s), one per group column.",
                           "Single value is applied to all. [default: %default]")),
  make_option(c("--group_colors"),  type = "character", default = NULL,
              help = paste("Semicolon-separated colour specs, one per group column.",
                           "Named: \"Adult=#8DD3C7,Pediatric=#FFFFB3;...\"",
                           "Positional: \"#E41A1C,#377EB8;#4DAF4A,#984EA3\"",
                           "Use empty block to skip a column: \"Adult=#8DD3C7;\"",
                           "Unspecified groups fall back to --group_palette."))
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$seg)) stop("--seg is required")
if (!is.null(opt$metadata) && is.null(opt$group_col))
  stop("--group_col is required when --metadata is provided")

# ---------------------------------------------------------------------------
# Parse multi-value group options
# ---------------------------------------------------------------------------
group_cols <- if (!is.null(opt$group_col))
  trimws(strsplit(opt$group_col, ",")[[1]]) else character(0)

if (length(group_cols) > 3)
  stop("--group_col: maximum 3 columns allowed (got ", length(group_cols), ")")

n_group_cols <- length(group_cols)

# Parse palettes (one per group col, or recycle single value)
palette_list <- if (n_group_cols > 0) {
  palettes <- trimws(strsplit(opt$group_palette, ",")[[1]])
  if (length(palettes) == 1) rep(palettes, n_group_cols) else palettes
} else character(0)

if (length(palette_list) < n_group_cols)
  stop("--group_palette: supply 1 palette or one per --group_col (",
       n_group_cols, " needed, ", length(palette_list), " given)")

# Parse manual colours (semicolon-separated blocks, comma within each block)
# Supports named entries: "Adult=#8DD3C7,Pediatric=#FFFFB3" or positional: "#E41A1C,#377EB8"
color_blocks <- if (!is.null(opt$group_colors)) {
  lapply(strsplit(opt$group_colors, ";")[[1]], function(b) {
    entries <- trimws(strsplit(b, ",")[[1]])
    if (length(entries) == 0 || all(entries == "")) return(NULL)
    has_eq <- grepl("=", entries, fixed = TRUE)
    if (any(has_eq)) {
      parts <- strsplit(entries[has_eq], "=", fixed = TRUE)
      nms   <- trimws(sapply(parts, `[`, 1))
      vals  <- trimws(sapply(parts, `[`, 2))
      setNames(vals, nms)
    } else {
      entries
    }
  })
} else vector("list", n_group_cols)

if (length(color_blocks) < n_group_cols)
  color_blocks <- c(color_blocks,
                    vector("list", n_group_cols - length(color_blocks)))

# ---------------------------------------------------------------------------
# Chromosome sizes
# ---------------------------------------------------------------------------
chr_sizes_hg19 <- c(
  chr1  = 249250621, chr2  = 243199373, chr3  = 198022430,
  chr4  = 191154276, chr5  = 180915260, chr6  = 171115067,
  chr7  = 159138663, chr8  = 146364022, chr9  = 141213431,
  chr10 = 135534747, chr11 = 135006516, chr12 = 133851895,
  chr13 = 115169878, chr14 = 107349540, chr15 = 102531392,
  chr16 = 90354753,  chr17 = 81195210,  chr18 = 78077248,
  chr19 = 59128983,  chr20 = 63025520,  chr21 = 48129895,
  chr22 = 51304566,  chrX  = 155270560, chrY  = 59373566
)

chr_sizes_hg38 <- c(
  chr1  = 248956422, chr2  = 242193529, chr3  = 198295559,
  chr4  = 190214555, chr5  = 181538259, chr6  = 170805979,
  chr7  = 159345973, chr8  = 145138636, chr9  = 138394717,
  chr10 = 133797422, chr11 = 135086622, chr12 = 133275309,
  chr13 = 114364328, chr14 = 107043718, chr15 = 101991189,
  chr16 = 90338345,  chr17 = 83257441,  chr18 = 80373285,
  chr19 = 58617616,  chr20 = 64444167,  chr21 = 46709983,
  chr22 = 50818468,  chrX  = 156040895, chrY  = 57227415
)

chr_sizes <- switch(opt$genome,
  hg19 = chr_sizes_hg19,
  hg38 = chr_sizes_hg38,
  stop("--genome must be hg19 or hg38")
)

# ---------------------------------------------------------------------------
# Apply sex chromosome filters
# ---------------------------------------------------------------------------
exclude_chrs <- character(0)
if (isTRUE(opt$no_sex_chr)) {
  exclude_chrs <- c("chrX", "chrY"); message("Excluding chrX and chrY")
} else {
  if (isTRUE(opt$no_chrX)) { exclude_chrs <- c(exclude_chrs, "chrX"); message("Excluding chrX") }
  if (isTRUE(opt$no_chrY)) { exclude_chrs <- c(exclude_chrs, "chrY"); message("Excluding chrY") }
}
chr_sizes <- chr_sizes[!names(chr_sizes) %in% exclude_chrs]
chr_order  <- names(chr_sizes)

chr_offsets <- c(0, cumsum(as.numeric(chr_sizes))[-length(chr_sizes)])
names(chr_offsets) <- chr_order

# ---------------------------------------------------------------------------
# Read & normalise segmentation file
# ---------------------------------------------------------------------------
message("Reading segmentation file: ", opt$seg)
delim <- if (grepl("\\.csv$", opt$seg, ignore.case = TRUE)) "," else "\t"
seg   <- read_delim(opt$seg, delim = delim, show_col_types = FALSE, progress = FALSE)

col_map <- list(
  ID        = c("ID", "Sample", "sample", "id", "SAMPLE_ID", "SampleID"),
  chrom     = c("chrom", "Chr", "CHR", "chromosome", "Chromosome", "CHROM"),
  loc.start = c("loc.start", "Start", "START", "start", "chromStart"),
  loc.end   = c("loc.end",   "End",   "END",   "end",   "chromEnd"),
  seg.mean  = c("seg.mean",  "Seg_Mean", "SEG_MEAN", "mean", "Mean",
                "log2_ratio", "Log2Ratio", "value")
)

rename_col <- function(df, std_name, candidates) {
  hit <- intersect(candidates, names(df))[1]
  if (is.na(hit)) stop("Cannot find column for '", std_name,
                        "'. Expected one of: ", paste(candidates, collapse = ", "))
  if (hit != std_name) {
    names(df)[names(df) == hit] <- std_name
    message("  Mapped column '", hit, "' -> '", std_name, "'")
  }
  df
}

for (nm in names(col_map)) seg <- rename_col(seg, nm, col_map[[nm]])

seg <- seg %>%
  mutate(chrom     = as.character(chrom),
         loc.start = as.numeric(loc.start),
         loc.end   = as.numeric(loc.end),
         seg.mean  = as.numeric(seg.mean))

seg$chrom <- ifelse(startsWith(seg$chrom, "chr"), seg$chrom, paste0("chr", seg$chrom))
seg       <- seg %>% filter(chrom %in% chr_order)

all_seg_ids <- unique(seg$ID)
n_samples   <- length(all_seg_ids)
message(sprintf("  %d samples, %d segments found", n_samples, nrow(seg)))

# ---------------------------------------------------------------------------
# Read metadata & build per-column group maps
# ---------------------------------------------------------------------------
# group_maps  : list of named vectors (sample_id -> group label), one per col
# group_palettes: list of named vectors (group label -> hex), one per col
group_maps     <- vector("list", n_group_cols)
group_palettes <- vector("list", n_group_cols)
y_label_map    <- NULL

if (!is.null(opt$metadata) && (n_group_cols > 0 || !is.null(opt$label_col))) {
  message("Reading metadata: ", opt$metadata)
  meta <- read_csv(opt$metadata, show_col_types = FALSE, progress = FALSE)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)

  # ---------------------------------------------------------------------------
  # Auto-detect the ID column that best matches the segmentation sample IDs.
  # Priority order: user-specified > Sentrix_ID > other common names > first col
  # ---------------------------------------------------------------------------
  id_candidates <- c(
    opt$meta_id_col,                              # user-specified (may be NULL)
    "Sentrix_ID", "sentrix_id", "SENTRIX_ID",     # methylation array standard
    "Basename", "basename",                        # full IDAT basename
    "ID", "Sample_ID", "SampleID", "Sample",
    "sample", "id", "SAMPLE_ID"
  )
  id_candidates <- id_candidates[!is.null(id_candidates) & !is.na(id_candidates)]

  id_col <- NA_character_
  for (cand in id_candidates) {
    if (cand %in% names(meta)) {
      # Check if values in this column overlap with the seg IDs
      overlap <- sum(as.character(meta[[cand]]) %in% all_seg_ids)
      if (overlap > 0) {
        id_col <- cand
        break
      }
    }
  }

  # If no overlap found by priority, fall back to first column with any match
  if (is.na(id_col)) {
    for (cn in names(meta)) {
      if (sum(as.character(meta[[cn]]) %in% all_seg_ids) > 0) {
        id_col <- cn; break
      }
    }
  }

  # Last resort: first column
  if (is.na(id_col)) {
    id_col <- names(meta)[1]
    warning("No metadata column matched segmentation sample IDs; ",
            "using first column '", id_col, "' as sample ID column.")
  }

  message(sprintf("  Matching metadata to segments via column '%s'", id_col))
  meta[[id_col]] <- as.character(meta[[id_col]])

  # Check match rate
  n_matched   <- sum(meta[[id_col]] %in% all_seg_ids)
  n_unmatched <- nrow(meta) - n_matched
  message(sprintf("  %d / %d metadata rows matched to segmentation samples",
                  n_matched, nrow(meta)))
  if (n_unmatched > 0)
    message(sprintf("  %d metadata row(s) had no matching segment — ignored", n_unmatched))
  if (n_matched == 0)
    stop("No metadata rows matched any sample ID in the segmentation file.\n",
         "  Seg IDs (first 5): ", paste(head(all_seg_ids, 5), collapse = ", "), "\n",
         "  Metadata '", id_col, "' (first 5): ",
         paste(head(meta[[id_col]], 5), collapse = ", "))

  meta_matched <- meta[meta[[id_col]] %in% all_seg_ids, ]

  # Validate group columns exist
  bad_cols <- setdiff(group_cols, names(meta_matched))
  if (length(bad_cols) > 0)
    stop("Group column(s) not found in metadata: ",
         paste(bad_cols, collapse = ", "),
         "\nAvailable columns: ", paste(names(meta_matched), collapse = ", "))

  # Build one group_map and colour palette per group column
  # user_colors can be: NULL, unnamed vector (positional), or named vector (by label)
  build_palette <- function(groups, user_colors, palette_name) {
    n_grp <- length(groups)

    # Helper: generate n colours from the RColorBrewer palette
    pal_colors <- function(n) {
      if (!requireNamespace("RColorBrewer", quietly = TRUE))
        stop("Install 'RColorBrewer' or supply --group_colors")
      max_pal <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
      if (is.na(max_pal))
        stop("Unknown RColorBrewer palette: '", palette_name, "'")
      if (n <= max_pal) {
        RColorBrewer::brewer.pal(max(3, n), palette_name)[seq_len(n)]
      } else {
        colorRampPalette(RColorBrewer::brewer.pal(max_pal, palette_name))(n)
      }
    }

    if (!is.null(user_colors) && length(user_colors) > 0) {
      has_names <- !is.null(names(user_colors)) && any(nchar(names(user_colors)) > 0)
      if (has_names) {
        # Named colours: map by group label, fill gaps with palette
        result <- setNames(rep(NA_character_, n_grp), groups)
        for (g in groups) {
          if (g %in% names(user_colors)) result[g] <- user_colors[g]
        }
        missing <- groups[is.na(result)]
        if (length(missing) > 0) result[missing] <- pal_colors(length(missing))
        return(result)
      } else if (length(user_colors) >= n_grp) {
        # Positional colours (backward compatible)
        return(setNames(user_colors[seq_len(n_grp)], groups))
      }
    }

    setNames(pal_colors(n_grp), groups)
  }

  for (i in seq_len(n_group_cols)) {
    gc <- group_cols[i]
    meta_matched[[gc]] <- as.character(meta_matched[[gc]])

    gmap <- setNames(meta_matched[[gc]], meta_matched[[id_col]])
    group_maps[[i]] <- gmap

    groups <- sort(unique(gmap))   # alphabetical order for legend
    message(sprintf("  Column '%s': %d group(s): %s",
                    gc, length(groups), paste(groups, collapse = ", ")))

    group_palettes[[i]] <- build_palette(groups, color_blocks[[i]], palette_list[i])
  }

  # Build y-axis label map from --label_col
  if (!is.null(opt$label_col)) {
    if (!opt$label_col %in% names(meta_matched))
      stop("--label_col '", opt$label_col, "' not found in metadata.\n",
           "Available columns: ", paste(names(meta_matched), collapse = ", "))
    y_label_map <- setNames(as.character(meta_matched[[opt$label_col]]),
                             meta_matched[[id_col]])
    message(sprintf("  Using '%s' for y-axis labels", opt$label_col))
  }
}

# Label function for y-axis
y_label_fn <- if (!is.null(y_label_map)) {
  lm <- y_label_map
  function(x) ifelse(x %in% names(lm), lm[x], x)
} else {
  waiver()
}

# ---------------------------------------------------------------------------
# Bin the genome
# ---------------------------------------------------------------------------
bin_size <- opt$bin_kb * 1000L
message(sprintf("Binning genome into %d kb bins ...", opt$bin_kb))

bins <- data.frame(
  chrom     = rep(chr_order, ceiling(chr_sizes / bin_size)),
  bin_start = unlist(lapply(chr_sizes, function(sz) seq(1, sz, by = bin_size))),
  stringsAsFactors = FALSE
)
bins$bin_end    <- pmin(bins$bin_start + bin_size - 1, chr_sizes[bins$chrom])
bins$bin_idx    <- seq_len(nrow(bins))
bins$genome_mid <- chr_offsets[bins$chrom] + (bins$bin_start + bins$bin_end) / 2

# ---------------------------------------------------------------------------
# Assign segments to bins (overlap-weighted average)
# ---------------------------------------------------------------------------
message("Mapping segments to bins ...")

assign_segments <- function(seg_chr, bins_chr) {
  if (nrow(seg_chr) == 0 || nrow(bins_chr) == 0) return(NULL)
  result <- lapply(seq_len(nrow(bins_chr)), function(i) {
    bs   <- bins_chr$bin_start[i]
    be   <- bins_chr$bin_end[i]
    hits <- seg_chr[seg_chr$loc.start <= be & seg_chr$loc.end >= bs, , drop = FALSE]
    if (nrow(hits) == 0) return(NULL)
    ol  <- pmin(hits$loc.end, be) - pmax(hits$loc.start, bs)
    val <- sum(hits$seg.mean * ol) / sum(ol)
    data.frame(bin_idx    = bins_chr$bin_idx[i],
               genome_mid = bins_chr$genome_mid[i],
               seg.mean   = val,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

tile_list <- lapply(all_seg_ids, function(samp) {
  s    <- seg %>% filter(ID == samp)
  rows <- lapply(chr_order, function(ch)
    assign_segments(s %>% filter(chrom == ch), bins %>% filter(chrom == ch)))
  df <- do.call(rbind, rows)
  if (!is.null(df)) df$ID <- samp
  df
})
tiles <- bind_rows(tile_list)

# ---------------------------------------------------------------------------
# Sample ordering
# ---------------------------------------------------------------------------
if (!is.null(opt$sample_order) && file.exists(opt$sample_order)) {
  samp_order <- readLines(opt$sample_order)
  samp_order <- samp_order[samp_order %in% all_seg_ids]
  missing_s  <- setdiff(all_seg_ids, samp_order)
  if (length(missing_s) > 0)
    message("Samples not in order file (appended): ", paste(missing_s, collapse = ", "))
  samp_order <- c(samp_order, missing_s)

} else if (isTRUE(opt$sort_by_group) && length(group_maps) > 0 && !is.null(group_maps[[1]])) {
  # Sort by all group columns in the order they appear in the metadata
  # (first occurrence order, not alphabetical)
  grp_df <- data.frame(ID = all_seg_ids, stringsAsFactors = FALSE)
  for (i in seq_len(n_group_cols)) {
    gmap_i  <- group_maps[[i]]
    col_nm  <- paste0("G", i)
    labels  <- ifelse(all_seg_ids %in% names(gmap_i), gmap_i[all_seg_ids], "Unknown")
    # Factor levels follow the first-occurrence order in the metadata vector
    data_lvls <- unique(c(gmap_i, if (any(labels == "Unknown")) "Unknown"))
    grp_df[[col_nm]] <- factor(labels, levels = data_lvls)
  }
  sort_cols  <- paste0("G", seq_len(n_group_cols))
  grp_df     <- grp_df[do.call(order, grp_df[, sort_cols, drop = FALSE]), ]
  samp_order <- grp_df$ID

} else {
  samp_order <- sort(all_seg_ids)
}

# Reverse so first sample appears at top of heatmap
tiles$ID <- factor(tiles$ID, levels = rev(samp_order))

# ---------------------------------------------------------------------------
# Pre-compute group boundary rectangles
# Each contiguous run of the same group label gets one rect spanning the full
# x-axis width.  y coordinates are in ggplot discrete-scale integer space:
#   samp_order[k] (1 = top) maps to y integer = n_samples - k + 1
# ---------------------------------------------------------------------------
group_rects <- if (n_group_cols > 0) {
  n_s <- length(samp_order)
  lapply(seq_len(n_group_cols), function(i) {
    gmap_i <- group_maps[[i]]
    labels  <- ifelse(samp_order %in% names(gmap_i), gmap_i[samp_order], "Unknown")
    rle_res    <- rle(labels)
    run_ends   <- cumsum(rle_res$lengths)
    run_starts <- c(1L, head(run_ends, -1L) + 1L)
    data.frame(
      xmin  = -Inf,
      xmax  =  Inf,
      ymin  = (n_s - run_ends   + 1L) - 0.5,
      ymax  = (n_s - run_starts + 1L) + 0.5,
      stringsAsFactors = FALSE
    )
  })
} else NULL

# ---------------------------------------------------------------------------
# Cap seg.mean; build full sample × bin grid
# ---------------------------------------------------------------------------
cap              <- opt$cap
tiles$seg_capped <- pmin(pmax(tiles$seg.mean, -cap), cap)

all_bins    <- bins[, c("bin_idx", "genome_mid")]
all_samples <- data.frame(ID = factor(samp_order, levels = rev(samp_order)),
                           stringsAsFactors = FALSE)
full_grid   <- merge(all_samples, all_bins, by = NULL)
tiles_full  <- merge(full_grid, tiles[, c("ID", "bin_idx", "seg_capped")],
                     by = c("ID", "bin_idx"), all.x = TRUE)

# ---------------------------------------------------------------------------
# Chromosome boundary / label data
# ---------------------------------------------------------------------------
chr_bounds       <- data.frame(chrom = chr_order,
                                xmin  = chr_offsets,
                                xmax  = chr_offsets + as.numeric(chr_sizes),
                                stringsAsFactors = FALSE)
chr_bounds$xmid  <- (chr_bounds$xmin + chr_bounds$xmax) / 2
chr_bounds$label <- sub("chr", "", chr_bounds$chrom)
chr_bounds$shade <- ifelse(seq_len(nrow(chr_bounds)) %% 2 == 0, "#E8E8E8", "#F8F8F8")

# ---------------------------------------------------------------------------
# Build main heatmap
# ---------------------------------------------------------------------------
message("Building plot ...")

plot_height <- if (opt$height > 0) opt$height else max(4, 0.2 * n_samples + 2)

use_newscale <- requireNamespace("ggnewscale", quietly = TRUE)
label_size   <- max(7, min(11, 150 / n_samples))
has_sidebar  <- n_group_cols > 0

# ggplot2 >= 3.5.0 supports per-guide position; use it to keep the colorbar
# at the bottom while group-colour legends are collected on the right.
colorbar_guide <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
  guide_colorbar(barwidth       = 12, barheight = 0.7,
                 title.position = "top", title.hjust = 0.5,
                 ticks.colour   = "grey30", frame.colour = "grey30",
                 position       = "bottom")
} else {
  guide_colorbar(barwidth       = 12, barheight = 0.7,
                 title.position = "top", title.hjust = 0.5,
                 ticks.colour   = "grey30", frame.colour = "grey30")
}

if (use_newscale) {
  p_heat <- ggplot(tiles_full, aes(x = genome_mid, y = ID)) +
    geom_rect(data = chr_bounds,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = shade),
              inherit.aes = FALSE, alpha = 0.7) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    geom_tile(aes(fill = seg_capped), width = bin_size, height = 0.9)
} else {
  message("ggnewscale not available — chromosome bands omitted")
  p_heat <- ggplot(tiles_full, aes(x = genome_mid, y = ID)) +
    geom_tile(aes(fill = seg_capped), width = bin_size, height = 0.9)
}

p_heat <- p_heat +
  scale_fill_gradient2(
    low      = opt$color_loss,
    mid      = opt$color_mid,
    high     = opt$color_gain,
    midpoint = 0,
    limits   = c(-cap, cap),
    oob      = squish,
    na.value = opt$color_na,
    name     = "Seg. mean\n(log2 ratio)",
    guide    = colorbar_guide
  ) +
  geom_vline(data = chr_bounds[-1, ], aes(xintercept = xmin),
             colour = "grey40", linewidth = 0.3) +
  scale_x_continuous(breaks = chr_bounds$xmid, labels = chr_bounds$label,
                     expand = c(0, 0)) +
  labs(title = opt$title, x = "Chromosome", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x     = element_text(size = 10, colour = "grey20", face = "bold"),
    # Show y-axis labels on heatmap only when there is no sidebar
    axis.text.y     = if (isTRUE(opt$no_labels) || has_sidebar) element_blank()
                      else element_text(size = label_size, colour = "grey10"),
    axis.ticks.x    = element_blank(),
    panel.grid      = element_blank(),
    panel.border    = element_rect(colour = "grey30", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9,  face = "bold"),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5,
                                   margin = margin(b = 8)),
    plot.margin     = margin(10, 10, 10, 2)
  )

# Apply y-axis label mapping (--label_col)
if (!is.null(y_label_map)) {
  p_heat <- p_heat + scale_y_discrete(labels = y_label_fn)
}

# ---------------------------------------------------------------------------
# Overlay group boundary rectangles on the heatmap
# col 1 → thick solid black  (outer / Phase-level divisions)
# col 2 → medium dashed      (inner / SubClass-level divisions)
# col 3 → thin dotted        (finest divisions)
# ---------------------------------------------------------------------------
if (!is.null(group_rects)) {
  rect_styles <- list(
    list(color = "black",   linewidth = 1.1, linetype = "solid"),
    list(color = "grey25",  linewidth = 0.6, linetype = "dashed"),
    list(color = "grey50",  linewidth = 0.4, linetype = "dotted")
  )
  for (i in seq_len(n_group_cols)) {
    if (!is.null(group_rects[[i]]) && nrow(group_rects[[i]]) > 0) {
      rs <- rect_styles[[min(i, 3L)]]
      p_heat <- p_heat +
        geom_rect(data        = group_rects[[i]],
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  inherit.aes = FALSE,
                  fill        = NA,
                  color       = rs$color,
                  linewidth   = rs$linewidth,
                  linetype    = rs$linetype)
    }
  }
}

# ---------------------------------------------------------------------------
# Build group annotation sidebars (one per group column)
# ---------------------------------------------------------------------------
make_sidebar <- function(gmap, gcolors, col_title, samp_order,
                         show_sample_labels, label_size, n_samples,
                         y_label_fn = waiver()) {
  ann_df <- data.frame(
    ID    = factor(samp_order, levels = rev(samp_order)),
    Group = ifelse(samp_order %in% names(gmap), gmap[samp_order], "Unknown"),
    x     = 1L,
    stringsAsFactors = FALSE
  )

  # Ensure Unknown is in the palette
  if ("Unknown" %in% ann_df$Group && !"Unknown" %in% names(gcolors))
    gcolors <- c(gcolors, Unknown = "#AAAAAA")

  ann_df$Group <- factor(ann_df$Group, levels = names(gcolors))

  key_h <- unit(max(3, min(10, 180 / n_samples)) * 0.12, "cm")

  ggplot(ann_df, aes(x = x, y = ID, fill = Group)) +
    geom_tile(width = 1, height = 0.9) +
    scale_fill_manual(
      values   = gcolors,
      name     = col_title,
      na.value = "#AAAAAA",
      guide    = guide_legend(
        title.position = "top",
        title.hjust    = 0,
        keywidth       = unit(0.4, "cm"),
        keyheight      = key_h,
        ncol           = max(1, ceiling(length(gcolors) / 15))
      )
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(labels = y_label_fn) +
    labs(x = col_title) +
    theme_void() +
    theme(
      legend.position      = "right",
      legend.justification = c(0, 0.5),
      legend.title    = element_text(size = 10, face = "bold", hjust = 0),
      legend.text     = element_text(size = 9),
      # Only the leftmost sidebar carries sample labels
      axis.text.y     = if (show_sample_labels)
                          element_text(size = label_size, colour = "grey10",
                                       hjust = 1, margin = margin(r = 2))
                        else
                          element_blank(),
      axis.title.x    = element_text(size = 9, face = "bold", angle = -90,
                                      vjust = 0.5, margin = margin(t = 4)),
      plot.margin     = margin(10, 2, 10, if (show_sample_labels) 4 else 1)
    )
}

# ---------------------------------------------------------------------------
# Combine sidebars + heatmap with patchwork
# ---------------------------------------------------------------------------
final_plot <- p_heat

if (has_sidebar) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("Package 'patchwork' is required for group sidebars.\n",
         "Install with: install.packages('patchwork')")

  sidebar_plots <- lapply(seq_len(n_group_cols), function(i) {
    make_sidebar(
      gmap               = group_maps[[i]],
      gcolors            = group_palettes[[i]],
      col_title          = group_cols[i],
      samp_order         = samp_order,
      show_sample_labels = (i == 1) && !isTRUE(opt$no_labels),
      label_size         = label_size,
      n_samples          = n_samples,
      y_label_fn         = y_label_fn
    )
  })

  # Width fractions: each sidebar ~2.5%, heatmap gets the rest
  sidebar_frac <- 0.025
  heat_frac    <- 1 - n_group_cols * sidebar_frac
  widths       <- c(rep(sidebar_frac, n_group_cols), heat_frac)

  # Sidebar legends stay right; heatmap colorbar stays bottom via its guide position
  final_plot <- patchwork::wrap_plots(c(sidebar_plots, list(p_heat)),
                                      ncol   = n_group_cols + 1,
                                      widths = widths) +
    patchwork::plot_layout(guides = "collect")
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
message(sprintf("Saving %s  (%.0f x %.1f in) ...", opt$output, opt$width, plot_height))

ext <- tolower(tools::file_ext(opt$output))
switch(ext,
  pdf = ggsave(opt$output, final_plot, width = opt$width, height = plot_height, device = "pdf", limitsize = FALSE),
  png = ggsave(opt$output, final_plot, width = opt$width, height = plot_height,
               dpi = opt$dpi, device = "png", limitsize = FALSE),
  svg = ggsave(opt$output, final_plot, width = opt$width, height = plot_height, device = "svg", limitsize = FALSE),
  stop("Unsupported format: '", ext, "'. Use .pdf, .png, or .svg")
)

message("Done.")
