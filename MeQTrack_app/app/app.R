# app/app.R — Shiny entrypoint for MeQTrack.
#
# Wave 2 shell:
#   - bslib Bootstrap 5 theme
#   - three-pane layout: sidebar (inputs) | main (tabbed content) |
#     footer (status strip)
#   - Samplesheet tab wired to the real samplesheet_module
#   - QC / Dim-Reduction / CNV / Report tabs are placeholders until their
#     respective waves
#
# Launch:
#   From the project root:
#     R -e 'shiny::runApp("app", launch.browser = TRUE, host = "127.0.0.1")'
#   Or double-click meqtrack.command (macOS) / meqtrack.bat (Windows).

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(shinyFiles)
  library(DT)
})

# ---------------------------------------------------------------------------
# Module sources
# ---------------------------------------------------------------------------
# shiny::runApp("app") sets the working directory to the app/ folder
# before sourcing this file, so module paths are relative to app/.
# project_root() in workspace.R walks up one level for pipeline paths.
source(file.path("R", "theme.R"),              local = FALSE)
source(file.path("R", "workspace.R"),          local = FALSE)
source(file.path("R", "validators.R"),         local = FALSE)
source(file.path("R", "array_detect.R"),       local = FALSE)
source(file.path("R", "samplesheet_module.R"), local = FALSE)
source(file.path("R", "pipeline_bridge.R"),    local = FALSE)
source(file.path("R", "run_controller.R"),     local = FALSE)
source(file.path("R", "results_loader.R"),     local = FALSE)
source(file.path("R", "past_runs_module.R"),   local = FALSE)
source(file.path("R", "settings_module.R"),    local = FALSE)
source(file.path("R", "qc_module.R"),          local = FALSE)
source(file.path("R", "dimred_module.R"),      local = FALSE)
source(file.path("R", "cnv_module.R"),         local = FALSE)
source(file.path("R", "report_module.R"),      local = FALSE)
source(file.path("R", "help_module.R"),        local = FALSE)

# ---------------------------------------------------------------------------
# First-launch workspace creation
# ---------------------------------------------------------------------------
WORKSPACE_PATH <- ensure_workspace()

# Expose run outputs (PDFs, interactive HTML plots, etc.) to the browser
# under a stable URL prefix. Every file under <workspace>/runs/ is reachable
# at /runs/<path>. Result modules build iframe src URLs from here.
addResourcePath("runs", file.path(WORKSPACE_PATH, "runs"))

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
# Build module UI pieces once so we can drop slots into multiple layout
# positions (sidebar vs main vs metadata card).
ss_ui        <- samplesheet_ui("samplesheet")
run_ui       <- run_controller_ui("run")
past_runs_ui <- past_runs_module_ui("past_runs")
settings_ui  <- settings_module_ui("settings")
qc_ui        <- qc_module_ui("qc")
dimred_ui    <- dimred_module_ui("dimred")
cnv_ui       <- cnv_module_ui("cnv")
report_ui    <- report_module_ui("report")
help_ui      <- help_module_ui("help")

ui <- bslib::page_navbar(
  id = "main_nav",
  title = tags$span(
    tags$strong("MeQTrack"),
    tags$span(" — methylation analysis",
              style = "color: #134e4a; font-size: 0.85rem; margin-left: .5rem; font-weight: 400;")
  ),
  theme = APP_THEME,
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "meqtrack.css")
  ),
  window_title = "MeQTrack",
  fillable = TRUE,
  # -----------------------------------------------------------------
  # Samplesheet tab (Wave 2)
  # -----------------------------------------------------------------
  bslib::nav_panel(
    title = "Samplesheet",
    icon = icon("table"),
    bslib::layout_sidebar(
      sidebar = bslib::sidebar(
        width = 320,
        ss_ui$sidebar,
        hr(),
        div(
          class = "text-muted small",
          tags$strong("Workspace:"), br(),
          tags$code(WORKSPACE_PATH)
        )
      ),
      bslib::card(
        bslib::card_header("Samples"),
        bslib::card_body(ss_ui$main, fillable = FALSE)
      ),
      bslib::card(
        bslib::card_header("Metadata"),
        bslib::card_body(ss_ui$metadata, fillable = FALSE)
      )
    )
  ),
  # -----------------------------------------------------------------
  # Run tab (Wave 3)
  # -----------------------------------------------------------------
  bslib::nav_panel(
    title = "Run",
    icon = icon("play"),
    # Settings and the pipeline controls/stages sit side by side so both
    # stay visible without scrolling; the log + actions span full width
    # below. layout_columns stacks the columns on narrow screens.
    bslib::layout_columns(
      col_widths = c(7, 5),
      bslib::card(
        bslib::card_header("Settings"),
        bslib::card_body(settings_ui, fillable = FALSE)
      ),
      bslib::card(
        bslib::card_header("Pipeline run"),
        bslib::card_body(run_ui$controls, fillable = FALSE)
      )
    ),
    run_ui$log,
    run_ui$actions
  ),
  # -----------------------------------------------------------------
  # Past runs library (Wave 6 Theme 6b)
  # -----------------------------------------------------------------
  bslib::nav_panel(
    title = "Past runs",
    icon = icon("clock-rotate-left"),
    bslib::card(
      bslib::card_header("Past runs"),
      bslib::card_body(past_runs_ui, fillable = FALSE)
    )
  ),
  # -----------------------------------------------------------------
  # Result tabs
  # -----------------------------------------------------------------
  bslib::nav_panel(
    title = "QC", icon = icon("chart-column"),
    bslib::card(
      bslib::card_header("QC review"),
      bslib::card_body(qc_ui, fillable = FALSE)
    )
  ),
  bslib::nav_panel(
    title = "Dim. reduction", icon = icon("diagram-project"),
    bslib::card(
      bslib::card_header("t-SNE / UMAP / Clustering / Reference projection"),
      bslib::card_body(dimred_ui, fillable = FALSE)
    )
  ),
  bslib::nav_panel(
    title = "CNV", icon = icon("dna"),
    bslib::card(
      bslib::card_header("CNV review"),
      bslib::card_body(cnv_ui, fillable = FALSE)
    )
  ),
  bslib::nav_panel(
    title = "Report", icon = icon("file-lines"),
    bslib::card(
      bslib::card_header("Report"),
      bslib::card_body(report_ui, fillable = FALSE)
    )
  ),
  bslib::nav_panel(
    title = "Help", icon = icon("circle-question"),
    bslib::card(
      bslib::card_header("Getting started"),
      bslib::card_body(help_ui, fillable = FALSE)
    )
  ),
  # -----------------------------------------------------------------
  # Footer / status strip
  # -----------------------------------------------------------------
  bslib::nav_spacer(),
  bslib::nav_item(
    uiOutput("header_status", inline = TRUE)
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  ss_state <- samplesheet_server(
    "samplesheet",
    workspace     = function() WORKSPACE_PATH,
    project_root_ = project_root
  )

  # Past runs library — emits a run_dir path when the user clicks "Open".
  past_run_selected <- past_runs_module_server(
    "past_runs",
    workspace = function() WORKSPACE_PATH
  )

  # Settings — five tunable QC + dim. reduction parameters. Repopulates
  # from a past run's manifest when the user attaches one.
  parameters <- settings_module_server(
    "settings",
    attach_run = past_run_selected
  )

  run_state <- run_controller_server(
    "run",
    ss_state      = ss_state,
    workspace     = function() WORKSPACE_PATH,
    project_root_ = project_root,
    attach_run    = past_run_selected,
    parameters    = parameters
  )

  # When the user attaches a past run, jump to the QC tab so the result is
  # immediately visible. This is the most informative landing — they can
  # navigate from there to Dim. reduction / CNV / Report as needed.
  observeEvent(past_run_selected(), {
    if (!is.null(past_run_selected())) {
      bslib::nav_select(id = "main_nav", selected = "QC", session = session)
    }
  }, ignoreInit = TRUE)

  # Wave 4: result views. results_loader emits a bundle when run_state flips
  # to COMPLETED; modules render from that bundle and reset when it goes NULL.
  results <- results_loader_server("results", run_state)
  qc_module_server("qc",         results)
  dimred_module_server("dimred", results)
  cnv_module_server("cnv",       results)
  report_module_server("report", results)
  help_module_server("help")

  # Header status (top-right): once a run has started, surface its state;
  # otherwise fall back to samplesheet readiness.
  output$header_status <- renderUI({
    rs <- run_state()
    if (!identical(rs$state, RUN_STATE_IDLE)) {
      return(render_state_badge(rs$state, rs$stage,
                                started_at = NULL, ended_at = NULL))
    }
    st <- ss_state()
    if (isTRUE(st$valid)) {
      status_pill("ok", "Samplesheet OK",
                  icon = "circle-check", class = "me-2")
    } else if (!is.null(st$samplesheet_path)) {
      status_pill("warning", "Samplesheet has issues",
                  icon = "triangle-exclamation", class = "me-2")
    } else {
      status_pill("neutral", "No samplesheet loaded", class = "me-2")
    }
  })

  # Log the samplesheet state to the launcher terminal whenever it changes;
  # cheap observability before we have a proper status panel.
  observe({
    st <- ss_state()
    if (is.null(st$samplesheet_path)) return()
    message(sprintf(
      "[samplesheet] %s | array_type=%s | valid=%s | n_ok=%s | n_bad=%s",
      st$samplesheet_path,
      st$array_type,
      isTRUE(st$valid),
      if (is.null(st$validated_df)) NA else sum(st$validated_df$Status == STATUS_OK),
      if (is.null(st$validated_df)) NA else sum(st$validated_df$Status != STATUS_OK)
    ))
  })
}

shinyApp(ui, server)
