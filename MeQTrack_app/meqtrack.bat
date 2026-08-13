@echo off
REM meqtrack.bat - Windows launcher for MeQTrack.
REM
REM Double-click this file in File Explorer to start the app. A command
REM window opens, this script runs, Shiny starts, and your default browser
REM opens to the app on http://127.0.0.1:<port>.
REM
REM Keep the window open while the app is running. Close it to stop the
REM app.
REM
REM If you hit launcher-specific issues, fall back to launching directly
REM from PowerShell or Command Prompt:
REM     Rscript -e "shiny::runApp('app', launch.browser = TRUE, host = '127.0.0.1')"

setlocal

REM Resolve the script's own directory so the launcher works regardless
REM of where the user placed MeQTrack_app.
cd /d "%~dp0"

REM Silence renv's "project is out-of-sync" startup notice. setup.R installs
REM current package versions and re-snapshots, but renv can't fully pin the
REM GitHub (conumee2) / local-tarball (yamapData) packages, so it always flags
REM a residual lock-vs-library difference. It's cosmetic — every package the
REM app needs is installed — so we suppress the check to avoid alarming users.
set RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

echo ==============================================================
echo   MeQTrack - launching...
echo   Project: %CD%
echo ==============================================================

REM 1. Check Rscript is on PATH.
where Rscript >nul 2>&1
if errorlevel 1 (
  echo.
  echo   ERROR: Rscript not found on your PATH.
  echo.
  echo   MeQTrack needs R ^>= 4.4 installed. Install from:
  echo       https://cran.r-project.org/bin/windows/base/
  echo   and then re-run this launcher.
  echo.
  pause
  exit /b 1
)

REM 2. Check pandoc for the report generator.
where pandoc >nul 2>&1
if errorlevel 1 (
  echo.
  echo   WARNING: pandoc not found on your PATH.
  echo.
  echo   The HTML report generator falls back to a plain-text report when
  echo   pandoc is missing. Install from https://pandoc.org/installing.html
  echo   to enable the full report.
  echo.
  echo   Continuing anyway...
  echo.
)

REM 3. First-launch provisioning. The distribution zip ships renv.lock and
REM    renv\activate.R but NOT the renv\library\ cache (machine-specific,
REM    would balloon the zip). On first launch we populate the library via
REM    setup.R (idempotent — no-op when already provisioned).
REM
REM    We probe by trying to load `shiny`, NOT by checking if renv\library\
REM    exists. R keys its library by minor version (R-4.5 vs R-4.6 are
REM    separate subdirs), so a user who upgrades R between launches has a
REM    populated renv\library\ that still fails to load anything against
REM    the new R. requireNamespace catches that case; setup.R then
REM    populates the library against the current R.
set NEED_SETUP=0
Rscript -e "quit(status = !requireNamespace('shiny', quietly = TRUE))" >nul 2>&1
if errorlevel 1 set NEED_SETUP=1
if "%NEED_SETUP%"=="1" (
  echo.
  echo   First-time setup: installing R packages.
  echo   This takes 5-15 minutes on first run; subsequent launches are instant.
  echo   Watch this window for progress; the app starts automatically when done.
  echo.
  Rscript setup.R
  if errorlevel 1 (
    echo.
    echo   ERROR: setup.R failed. See the messages above.
    pause
    exit /b 1
  )
  echo.
  echo   Setup complete. Starting the app...
  echo.
)

REM 4. Start Shiny.
echo Starting Shiny server on http://127.0.0.1 ...
echo (Close this window to stop the app.)
echo.

Rscript -e "shiny::runApp('app', launch.browser = TRUE, host = '127.0.0.1')"

endlocal
