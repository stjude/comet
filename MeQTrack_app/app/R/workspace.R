# app/R/workspace.R
# ---------------------------------------------------------------------------
# Workspace + project path helpers.
#
# The "workspace" is the user-owned folder where samplesheets live, IDATs
# live, and run outputs are written. Default is ~/MeQTrack. It is separate
# from the app install (the MeQTrack_app project folder) so that upgrading
# the app doesn't touch the user's data.
#
# The "project root" is the folder containing pipeline/, app/, renv/, etc.
# We use it to resolve relative Basename entries that target `pipeline/data/
# example/...` (the bundled example ships that way, matching the CLI
# pipeline's behaviour of setwd(script_dir) at startup).
# ---------------------------------------------------------------------------

# Default workspace root. User can override via a settings panel in a later
# wave; for Wave 2 it's a hard default.
DEFAULT_WORKSPACE <- path.expand("~/MeQTrack")

#' Return the absolute path to the project root.
#'
#' app.R is launched from the project root by both the launcher script and
#' `shiny::runApp("app")`. We use getwd() and verify the `pipeline/` folder
#' exists alongside us.
project_root <- function() {
  cwd <- normalizePath(".", winslash = "/", mustWork = TRUE)
  # If we're inside app/ (rare), step up one.
  if (basename(cwd) == "app" && dir.exists(file.path(cwd, "..", "pipeline"))) {
    cwd <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
  }
  cwd
}

#' Return the absolute path to the bundled pipeline folder.
pipeline_dir <- function() {
  file.path(project_root(), "pipeline")
}

#' Return the absolute path to the user's workspace.
#' Creates it (with `samplesheets/`, `idats/`, `runs/` subfolders) if missing.
#'
#' @param path   Override for the default workspace root.
ensure_workspace <- function(path = DEFAULT_WORKSPACE) {
  path <- path.expand(path)
  for (sub in c("", "samplesheets", "idats", "runs")) {
    d <- if (nzchar(sub)) file.path(path, sub) else path
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
