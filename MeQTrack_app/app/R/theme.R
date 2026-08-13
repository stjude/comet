# app/R/theme.R
# ---------------------------------------------------------------------------
# Design tokens and shared styling helpers for MeQTrack.
#
# Exports:
#   COLORS              list of semantic color constants
#   APP_THEME           bslib::bs_theme passed to page_navbar(theme = ...)
#   theme_meqtrack_gg() ggplot2 theme matching the palette
#   plotly_defaults()   layers font + muted gridlines onto any plotly figure
#
# Scope: Phase A of the UI design refresh (palette + typography + plot
# harmonization). Layout restructuring and status-pill unification are
# Phase B.
# ---------------------------------------------------------------------------

COLORS <- list(
  surface_0    = "#ffffff",   # page background
  surface_1    = "#f8fafc",   # card interior
  surface_2    = "#e2e8f0",   # dividers / borders
  ink_900      = "#0f172a",   # body text
  ink_500      = "#475569",   # muted text / axis labels
  primary      = "#0d9488",   # teal — actions, active nav
  primary_700  = "#0f766e",   # teal deep — button hover, active underline
  primary_900  = "#134e4a",   # teal-ink — dark accent
  accent_light = "#3de0d5",   # navbar gradient start
  accent_soft  = "#99f6e4",   # navbar gradient end / card tints
  violet       = "#7c3aed",   # secondary accent for clustering
  pass         = "#10b981",   # QC pass (stays saturated)
  warn         = "#f59e0b",   # QC warn / cancellation (stays saturated)
  fail         = "#ef4444",   # QC fail / failed run (stays saturated)
  cnv_gain     = "#d97706",   # heatmap +seg.mean (stays saturated)
  cnv_loss     = "#0d9488"    # heatmap −seg.mean (teal, matches primary)
)

# Qualitative palette for scatter-plot metadata coloring (t-SNE / UMAP /
# 3D MDS). Teal-leading, staying within the theme family. Ordering
# balances perceptual separation for up to 8 groups.
SCATTER_PALETTE <- c(
  "#0d9488",  # teal (primary)
  "#a6611a",  # brown (CNV-gain counterpart)
  "#7c3aed",  # violet
  "#d97706",  # orange
  "#10b981",  # emerald green
  "#f59e0b",  # amber
  "#3de0d5",  # light teal
  "#ef4444"   # red (last — signal reserve)
)

# Two-stop palette for Pass/Fail QC binary coloring.
QC_BINARY_COLORS <- c(`FALSE` = "#ef4444", `TRUE` = "#0d9488")

APP_THEME <- bslib::bs_theme(
  version      = 5,
  primary      = COLORS$primary,
  success      = COLORS$pass,
  warning      = COLORS$warn,
  danger       = COLORS$fail,
  "body-color"   = COLORS$ink_900,
  "body-bg"      = COLORS$surface_0,
  "border-color" = COLORS$surface_2,
  base_font    = bslib::font_google("Inter"),
  heading_font = bslib::font_google("Inter"),
  code_font    = bslib::font_google("JetBrains Mono")
)

# ggplot theme. Uses "sans" rather than "Inter" for base_family because
# Inter isn't guaranteed to be installed as a system font (bslib loads it
# via webfont for the browser). Plotly conversion picks up Inter through
# plotly_defaults() below anyway.
theme_meqtrack_gg <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text              = ggplot2::element_text(color = COLORS$ink_900),
      axis.text         = ggplot2::element_text(color = COLORS$ink_500),
      axis.title        = ggplot2::element_text(color = COLORS$ink_500,
                                                size = base_size - 1),
      panel.grid.major  = ggplot2::element_line(color = COLORS$surface_2,
                                                linewidth = 0.3),
      panel.grid.minor  = ggplot2::element_blank(),
      strip.text        = ggplot2::element_text(color = COLORS$ink_500,
                                                size = base_size - 2,
                                                face = "bold"),
      legend.position   = "top",
      legend.title      = ggplot2::element_text(size = base_size - 2,
                                                color = COLORS$ink_500),
      legend.text       = ggplot2::element_text(size = base_size - 2),
      plot.title        = ggplot2::element_text(size = base_size + 2,
                                                face = "bold",
                                                color = COLORS$ink_900),
      plot.subtitle     = ggplot2::element_text(size = base_size - 1,
                                                color = COLORS$ink_500)
    )
}

# ---------------------------------------------------------------------------
# status_pill — single source of truth for Bootstrap state badges.
#
# Replaces inline tags$span(class = "badge bg-X") calls that used to live
# in run_controller, samplesheet_module, app.R, etc. Same visual output,
# one place to change colors / aliases / typography.
#
# Synonym aliases mean callers can use whichever vocabulary fits the call
# site (idle/neutral/pending all render gray; ok/done/success all render
# green) without per-site mapping.
# ---------------------------------------------------------------------------
STATUS_PILL_CLASSES <- c(
  ok        = "bg-success",
  done      = "bg-success",
  success   = "bg-success",
  completed = "bg-success",
  warn      = "bg-warning text-dark",
  warning   = "bg-warning text-dark",
  cancelled = "bg-warning text-dark",
  fail      = "bg-danger",
  failed    = "bg-danger",
  danger    = "bg-danger",
  running   = "bg-primary",
  active    = "bg-primary",
  info      = "bg-primary",
  pending   = "bg-secondary",
  idle      = "bg-secondary",
  neutral   = "bg-secondary"
)

#' Build a DT headerCallback that attaches native title attributes to
#' column headers based on a name → text dictionary. Browser-native
#' tooltips show on hover after a short delay; no Bootstrap init needed.
#' Pass the result as options = list(headerCallback = dt_header_tooltips(...)).
dt_header_tooltips <- function(tooltips) {
  if (!length(tooltips)) return(NULL)
  json <- jsonlite::toJSON(tooltips, auto_unbox = TRUE)
  DT::JS(sprintf(
    "function(thead) {
      var tt = %s;
      $('th', thead).each(function() {
        var n = $(this).text().trim();
        if (tt[n]) { $(this).attr('title', tt[n]); }
      });
    }", json
  ))
}

#' Render a status badge.
#'
#' @param state    one of names(STATUS_PILL_CLASSES); unknowns fall back
#'                 to neutral gray rather than erroring.
#' @param label    visible text. Defaults to state.
#' @param icon     optional Font Awesome icon name (e.g. "circle-check");
#'                 prepended with a half-em space.
#' @param class    extra utility classes (e.g. "me-2") appended verbatim.
status_pill <- function(state, label = NULL, icon = NULL, class = NULL) {
  cls <- STATUS_PILL_CLASSES[state]
  if (length(cls) == 0L || is.na(cls)) cls <- "bg-secondary"
  text <- if (is.null(label)) state else label
  shiny::tags$span(
    class = paste("badge", cls, class),
    if (!is.null(icon)) shiny::tagList(shiny::icon(icon), " "),
    text
  )
}

# Layer MeQTrack font + muted axis styling onto an existing plotly figure.
# Safe to call as the final step in a |>-chain — plotly::layout does a
# deep merge so per-axis titles set earlier (e.g. "Height") survive.
plotly_defaults <- function(fig) {
  axis_style <- list(
    gridcolor = COLORS$surface_2,
    linecolor = COLORS$surface_2,
    zerolinecolor = COLORS$surface_2,
    tickfont = list(color = COLORS$ink_500, size = 11),
    titlefont = list(color = COLORS$ink_500, size = 12)
  )
  fig |>
    plotly::layout(
      font = list(
        family = "Inter, -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
        color  = COLORS$ink_900,
        size   = 12
      ),
      paper_bgcolor = COLORS$surface_0,
      plot_bgcolor  = COLORS$surface_0,
      xaxis  = axis_style,
      yaxis  = axis_style,
      legend = list(
        bgcolor = "rgba(0,0,0,0)",
        font    = list(size = 11, color = COLORS$ink_500)
      )
    ) |>
    # Trim the modebar to a consistent set of buttons across every plot.
    # Hide the plotly logo (we already cite plotly in the docs) and drop
    # buttons that don't map cleanly to methylation data exploration.
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c(
        "select2d", "lasso2d", "autoScale2d",
        "hoverClosestCartesian", "hoverCompareCartesian",
        "toggleSpikelines"
      ),
      toImageButtonOptions = list(
        format = "png", filename = "meqtrack_plot",
        height = 600, width = 800, scale = 2
      )
    )
}
