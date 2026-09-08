# HPC job submission module for methylation array analysis pipeline

#' Generate HPC job submission scripts
#'
#' @param config Configuration list
#' @param opt Command line options
#' @param output_dir Output directory
#' @param scheduler Scheduler to use (slurm, pbs, lsf)
#' @return List of generated script paths
generate_hpc_scripts <- function(config, opt, output_dir, scheduler = "lsf") {
  message("Generating HPC job submission scripts for ", toupper(scheduler), " scheduler...")
  
  # Create directory for scripts
  scripts_dir <- file.path(output_dir, "hpc_scripts")
  dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
  
  script_paths <- list()
  
  # Generate wrapper script
  wrapper_script <- generate_wrapper_script(config, opt, scripts_dir, scheduler)
  script_paths$wrapper <- wrapper_script
  
  # Generate step scripts
  steps <- c("preprocess", "qc", "filtering", "dim_reduction", "cnv", "visualization")
  for (step in steps) {
    step_script <- generate_step_script(step, config, opt, scripts_dir, scheduler)
    script_paths[[step]] <- step_script
  }
  
  # Make scripts executable
  if (.Platform$OS.type == "unix") {
    invisible(sapply(unlist(script_paths), function(script) {
      system(paste("chmod +x", script))
    }))
  }
  
  # Generate submission commands file
  submission_commands <- generate_submission_commands(script_paths, config, scripts_dir)
  script_paths$submission_commands <- submission_commands
  
  # Generate LSF workflow script if using LSF
  if (scheduler == "lsf") {
    workflow_script <- generate_lsf_workflow(script_paths, config, output_dir)
    script_paths$lsf_workflow <- workflow_script
  }
  
  message("HPC job submission scripts generated in: ", scripts_dir)
  return(script_paths)
}

#' Generate wrapper script for running the entire pipeline
#'
#' @param config Configuration list
#' @param opt Command line options
#' @param scripts_dir Directory for scripts
#' @param scheduler Scheduler to use (slurm, pbs, lsf)
#' @return Path to generated script
generate_wrapper_script <- function(config, opt, scripts_dir, scheduler = "lsf") {
  # Determine script extension based on platform
  ext <- ifelse(.Platform$OS.type == "unix", "sh", "bat")
  
  # Set up script content based on platform
  if (.Platform$OS.type == "unix") {
    # Base script content for UNIX/Linux/Mac without scheduler headers
    base_script_content <- c(
      "#!/bin/bash",
      "",
      "# MeQTrack - Wrapper Script",
      paste0("# Generated: ", Sys.time()),
      "",
      "# Pipeline root directory",
      paste0("PIPELINE_DIR=\"", normalizePath(getwd()), "\""),
      "",
      "# Command line arguments",
      paste0("CONFIG=\"", opt$config, "\""),
      paste0("OUTPUT_DIR=\"", opt$output, "\""),
      paste0("ARRAY_TYPE=\"", opt$array_type, "\""),
      paste0("THREADS=", opt$threads),
      paste0("SAMPLE_SHEET=\"", config$sample_sheet, "\""),
      "",
      "# Create output directory",
      "mkdir -p \"$OUTPUT_DIR\"",
      "",
      "# Load modules (uncomment and adjust as needed)",
      "#module load r/4.2.0",
      "",
      "# Run pipeline",
      "Rscript \"$PIPELINE_DIR/methylation_pipeline.R\" \\",
      "  --config \"$CONFIG\" \\",
      "  --array_type \"$ARRAY_TYPE\" \\",
      "  --output \"$OUTPUT_DIR\" \\",
      "  --threads \"$THREADS\" \\",
      "  --input \"$SAMPLE_SHEET\"",
      "",
      "echo \"Pipeline completed.\"",
      ""
    )
    
    # Add scheduler-specific headers
    if (scheduler == "slurm") {
      # SLURM headers
      slurm_headers <- c(
        "#SBATCH --job-name=meqtrack_pipeline",
        "#SBATCH --nodes=1",
        "#SBATCH --mem=16g",
        "#SBATCH --time=24:00:00",
        paste0("#SBATCH --cpus-per-task=", opt$threads),
        "#SBATCH --output=meqtrack_pipeline_%j.out",
        "#SBATCH --error=meqtrack_pipeline_%j.err"
      )
      
      # Insert headers after shebang
      script_content <- c(
        base_script_content[1],
        slurm_headers,
        base_script_content[2:length(base_script_content)]
      )
    } else if (scheduler == "pbs") {
      # Apply PBS headers
      script_content <- add_pbs_headers(base_script_content, "pipeline", "32g", "24:00:00", opt$threads)
    } else if (scheduler == "lsf") {
      # Apply LSF headers
      lsf_headers <- c(
        "#BSUB -J meqtrack_pipeline",
        paste0("#BSUB -n ", opt$threads),
        "#BSUB -R \"rusage[mem=32GB]\"",
        "#BSUB -W 1440", # 24 hours in minutes
        "#BSUB -cwd $(pwd)",
        "#BSUB -o meqtrack_pipeline.%J.out",
        "#BSUB -e meqtrack_pipeline.%J.err"
      )
      
      # Insert headers after shebang
      script_content <- c(
        base_script_content[1],
        lsf_headers,
        base_script_content[2:length(base_script_content)]
      )
    } else {
      # Use base script without scheduler headers
      script_content <- base_script_content
    }
  } else {
    # Windows
    script_content <- c(
      "@echo off",
      "",
      ":: MeQTrack - Wrapper Script",
      paste0(":: Generated: ", Sys.time()),
      "",
      ":: Pipeline root directory",
      paste0("set PIPELINE_DIR=", normalizePath(getwd())),
      "",
      ":: Command line arguments",
      paste0("set CONFIG=\"", opt$config, "\""),
      paste0("set OUTPUT_DIR=\"", opt$output, "\""),
      paste0("set ARRAY_TYPE=\"", opt$array_type, "\""),
      paste0("set THREADS=", opt$threads),
      paste0("set SAMPLE_SHEET=\"", config$sample_sheet, "\""),
      "",
      ":: Create output directory",
      "if not exist %OUTPUT_DIR% mkdir %OUTPUT_DIR%",
      "",
      ":: Run pipeline",
      "Rscript \"%PIPELINE_DIR%\\methylation_pipeline.R\" ^",
      "  --config %CONFIG% ^",
      "  --array_type %ARRAY_TYPE% ^",
      "  --output %OUTPUT_DIR% ^",
      "  --threads %THREADS% ^",
      "  --input %SAMPLE_SHEET%",
      "",
      "echo Pipeline completed.",
      ""
    )
  }
  
  # Write script to file
  script_path <- file.path(scripts_dir, paste0("run_pipeline.", ext))
  writeLines(script_content, script_path)
  
  return(script_path)
}

#' Generate script for running a specific pipeline step
#'
#' @param step Pipeline step
#' @param config Configuration list
#' @param opt Command line options
#' @param scripts_dir Directory for scripts
#' @param scheduler Scheduler to use (slurm, pbs, lsf)
#' @return Path to generated script
generate_step_script <- function(step, config, opt, scripts_dir, scheduler = "lsf") {
  # Determine script extension based on platform
  ext <- ifelse(.Platform$OS.type == "unix", "sh", "bat")
  
  # Determine memory and time requirements based on step
  mem_req <- switch(step,
                   preprocess = "16g",
                   qc = "8g",
                   filtering = "8g",
                   dim_reduction = "16g",
                   cnv = "16g",
                   visualization = "8g",
                   "8g")  # default
  
  time_req <- switch(step,
                    preprocess = "4:00:00",
                    qc = "1:00:00",
                    filtering = "1:00:00",
                    dim_reduction = "2:00:00",
                    cnv = "8:00:00",
                    visualization = "2:00:00",
                    "2:00:00")  # default
  
  # Set up script content based on platform
  if (.Platform$OS.type == "unix") {
    # Base script content without scheduler headers
    base_script_content <- c(
      "#!/bin/bash",
      "",
      paste0("# MeQTrack - ", step, " Step"),
      paste0("# Generated: ", Sys.time()),
      "",
      "# Pipeline root directory",
      paste0("PIPELINE_DIR=\"", normalizePath(getwd()), "\""),
      "",
      "# Command line arguments",
      paste0("CONFIG=\"", opt$config, "\""),
      paste0("OUTPUT_DIR=\"", opt$output, "\""),
      paste0("ARRAY_TYPE=\"", opt$array_type, "\""),
      paste0("THREADS=", opt$threads),
      paste0("STEP=\"", step, "\""),
      paste0("SAMPLE_SHEET=\"", config$sample_sheet, "\""),
      "",
      "# Create output directory",
      "mkdir -p \"$OUTPUT_DIR\"",
      "",
      "# Load modules (uncomment and adjust as needed)",
      "#module load r/4.2.0",
      "",
      "# Run pipeline step",
      "Rscript \"$PIPELINE_DIR/methylation_pipeline.R\" \\",
      "  --config \"$CONFIG\" \\",
      "  --step \"$STEP\" \\",
      "  --array_type \"$ARRAY_TYPE\" \\",
      "  --output \"$OUTPUT_DIR\" \\",
      "  --threads \"$THREADS\" \\",
      "  --input \"$SAMPLE_SHEET\"",
      "",
      paste0("echo \"", step, " step completed.\""),
      ""
    )
    
    # Add scheduler-specific headers
    if (scheduler == "slurm") {
      # SLURM headers
      slurm_headers <- c(
        paste0("#SBATCH --job-name=meqtrack_", step),
        "#SBATCH --nodes=1",
        paste0("#SBATCH --mem=", mem_req),
        paste0("#SBATCH --time=", time_req),
        paste0("#SBATCH --cpus-per-task=", opt$threads),
        "#SBATCH --output=meqtrack_%j.out",
        "#SBATCH --error=meqtrack_%j.err"
      )
      
      # Insert headers after shebang
      script_content <- c(
        base_script_content[1],
        slurm_headers,
        base_script_content[2:length(base_script_content)]
      )
    } else if (scheduler == "pbs") {
      # Apply PBS headers
      script_content <- add_pbs_headers(base_script_content, step, mem_req, time_req, opt$threads)
    } else if (scheduler == "lsf") {
      # Apply LSF headers
      script_content <- add_lsf_headers(base_script_content, step, mem_req, time_req, opt$threads)
    } else {
      # Use base script without scheduler headers
      script_content <- base_script_content
    }
  } else {
    # Windows
    script_content <- c(
      "@echo off",
      "",
      paste0(":: MeQTrack - ", step, " Step"),
      paste0(":: Generated: ", Sys.time()),
      "",
      ":: Pipeline root directory",
      paste0("set PIPELINE_DIR=", normalizePath(getwd())),
      "",
      ":: Command line arguments",
      paste0("set CONFIG=\"", opt$config, "\""),
      paste0("set OUTPUT_DIR=\"", opt$output, "\""),
      paste0("set ARRAY_TYPE=\"", opt$array_type, "\""),
      paste0("set THREADS=", opt$threads),
      paste0("set STEP=\"", step, "\""),
      paste0("set SAMPLE_SHEET=\"", config$sample_sheet, "\""),
      "",
      ":: Create output directory",
      "if not exist %OUTPUT_DIR% mkdir %OUTPUT_DIR%",
      "",
      ":: Run pipeline step",
      "Rscript \"%PIPELINE_DIR%\\methylation_pipeline.R\" ^",
      "  --config %CONFIG% ^",
      "  --step %STEP% ^",
      "  --array_type %ARRAY_TYPE% ^",
      "  --output %OUTPUT_DIR% ^",
      "  --threads %THREADS% ^",
      "  --input %SAMPLE_SHEET%",
      "",
      paste0("echo ", step, " step completed."),
      ""
    )
  }
  
  # Write script to file
  script_path <- file.path(scripts_dir, paste0(step, "_step.", ext))
  writeLines(script_content, script_path)
  
  return(script_path)
}

#' Generate job submission commands for different HPC systems
#'
#' @param script_paths List of script paths
#' @param config Configuration list
#' @param scripts_dir Directory for scripts
#' @return Path to submission commands file
generate_submission_commands <- function(script_paths, config, scripts_dir) {
  # Create submission commands for SLURM
  slurm_commands <- c(
    "# SLURM submission commands",
    "# -------------------------",
    "",
    "# Run entire pipeline",
    paste0("sbatch ", basename(script_paths$wrapper)),
    "",
    "# Run individual steps",
    paste0("sbatch ", basename(script_paths$preprocess)),
    paste0("sbatch ", basename(script_paths$qc)),
    paste0("sbatch ", basename(script_paths$filtering)),
    paste0("sbatch ", basename(script_paths$dim_reduction)),
    paste0("sbatch ", basename(script_paths$cnv)),
    paste0("sbatch ", basename(script_paths$visualization)),
    "",
    "# Run steps with dependencies (run in sequence)",
    paste0("PREPROCESS_ID=$(sbatch ", basename(script_paths$preprocess), " | awk '{print $4}')"),
    paste0("QC_ID=$(sbatch --dependency=afterok:$PREPROCESS_ID ", basename(script_paths$qc), " | awk '{print $4}')"),
    paste0("FILTERING_ID=$(sbatch --dependency=afterok:$QC_ID ", basename(script_paths$filtering), " | awk '{print $4}')"),
    paste0("DIM_REDUCTION_ID=$(sbatch --dependency=afterok:$FILTERING_ID ", basename(script_paths$dim_reduction), " | awk '{print $4}')"),
    paste0("CNV_ID=$(sbatch --dependency=afterok:$FILTERING_ID ", basename(script_paths$cnv), " | awk '{print $4}')"),
    paste0("sbatch --dependency=afterok:$DIM_REDUCTION_ID,afterok:$CNV_ID ", basename(script_paths$visualization)),
    ""
  )
  
  # Create submission commands for PBS/Torque
  pbs_commands <- c(
    "# PBS/Torque submission commands",
    "# ------------------------------",
    "",
    "# Run entire pipeline",
    paste0("qsub ", basename(script_paths$wrapper)),
    "",
    "# Run individual steps",
    paste0("qsub ", basename(script_paths$preprocess)),
    paste0("qsub ", basename(script_paths$qc)),
    paste0("qsub ", basename(script_paths$filtering)),
    paste0("qsub ", basename(script_paths$dim_reduction)),
    paste0("qsub ", basename(script_paths$cnv)),
    paste0("qsub ", basename(script_paths$visualization)),
    "",
    "# Run steps with dependencies (run in sequence)",
    paste0("PREPROCESS_ID=$(qsub ", basename(script_paths$preprocess), ")"),
    paste0("QC_ID=$(qsub -W depend=afterok:$PREPROCESS_ID ", basename(script_paths$qc), ")"),
    paste0("FILTERING_ID=$(qsub -W depend=afterok:$QC_ID ", basename(script_paths$filtering), ")"),
    paste0("DIM_REDUCTION_ID=$(qsub -W depend=afterok:$FILTERING_ID ", basename(script_paths$dim_reduction), ")"),
    paste0("CNV_ID=$(qsub -W depend=afterok:$FILTERING_ID ", basename(script_paths$cnv), ")"),
    paste0("qsub -W depend=afterok:$DIM_REDUCTION_ID,afterok:$CNV_ID ", basename(script_paths$visualization)),
    ""
  )
  
  # Create submission commands for LSF
  lsf_commands <- c(
    "# LSF (bsub) submission commands",
    "# -----------------------------",
    "",
    "# Run entire pipeline",
    paste0("bsub < ", basename(script_paths$wrapper)),
    "",
    "# Run individual steps",
    paste0("bsub < ", basename(script_paths$preprocess)),
    paste0("bsub < ", basename(script_paths$qc)),
    paste0("bsub < ", basename(script_paths$filtering)),
    paste0("bsub < ", basename(script_paths$dim_reduction)),
    paste0("bsub < ", basename(script_paths$cnv)),
    paste0("bsub < ", basename(script_paths$visualization)),
    "",
    "# Run steps with dependencies (run in sequence)",
    paste0("PREPROCESS_ID=$(bsub < ", basename(script_paths$preprocess), " | grep -o \"<[0-9]*>\" | tr -d \"<>\")"),
    paste0("bsub -w \"done($PREPROCESS_ID)\" < ", basename(script_paths$qc)),
    paste0("QC_ID=$(bsub -w \"done($PREPROCESS_ID)\" < ", basename(script_paths$qc), " | grep -o \"<[0-9]*>\" | tr -d \"<>\")"),
    paste0("FILTERING_ID=$(bsub -w \"done($QC_ID)\" < ", basename(script_paths$filtering), " | grep -o \"<[0-9]*>\" | tr -d \"<>\")"),
    paste0("DIM_REDUCTION_ID=$(bsub -w \"done($FILTERING_ID)\" < ", basename(script_paths$dim_reduction), " | grep -o \"<[0-9]*>\" | tr -d \"<>\")"),
    paste0("CNV_ID=$(bsub -w \"done($FILTERING_ID)\" < ", basename(script_paths$cnv), " | grep -o \"<[0-9]*>\" | tr -d \"<>\")"),
    paste0("bsub -w \"done($DIM_REDUCTION_ID) && done($CNV_ID)\" < ", basename(script_paths$visualization)),
    ""
  )
  
  # Combine all commands
  all_commands <- c(
    "# MeQTrack - HPC Submission Commands",
    paste0("# Generated: ", Sys.time()),
    "# ---------------------------------------------------------",
    "",
    "# This file contains examples of submission commands for different HPC systems.",
    "# Adjust as needed for your specific environment.",
    "",
    slurm_commands,
    pbs_commands,
    lsf_commands
  )
  
  # Write commands to file
  commands_file <- file.path(scripts_dir, "submission_commands.txt")
  writeLines(all_commands, commands_file)
  
  return(commands_file)
}

#' Generate LSF workflow submission script
#' 
#' @param script_paths List of script paths
#' @param config Configuration list
#' @param output_dir Output directory
#' @param queue LSF queue name
#' @return Path to generated workflow script
generate_lsf_workflow <- function(script_paths, config, output_dir, queue = "normal") {
  # Get script directory
  scripts_dir <- file.path(output_dir, "hpc_scripts")
  
  # Create LSF workflow script
  lsf_workflow <- c(
    "#!/bin/bash",
    "# MeQTrack - LSF Workflow Submission Script",
    paste0("# Generated: ", Sys.time()),
    "",
    "# Step 1: Run preprocessing",
    paste0("PREP_ID=$(bsub -q ", queue, " < ", basename(script_paths$preprocess), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"Preprocessing job submitted with ID: $PREP_ID\"",
    "",
    "# Step 2: Run QC (depends on preprocessing)",
    paste0("QC_ID=$(bsub -q ", queue, " -w \"done($PREP_ID)\" < ", basename(script_paths$qc), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"QC job submitted with ID: $QC_ID\"",
    "",
    "# Step 3: Run filtering (depends on QC)",
    paste0("FILTER_ID=$(bsub -q ", queue, " -w \"done($QC_ID)\" < ", basename(script_paths$filtering), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"Filtering job submitted with ID: $FILTER_ID\"",
    "",
    "# Step 4: Run dimensionality reduction (depends on filtering)",
    paste0("DIM_ID=$(bsub -q ", queue, " -w \"done($FILTER_ID)\" < ", basename(script_paths$dim_reduction), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"Dimensionality reduction job submitted with ID: $DIM_ID\"",
    "",
    "# Step 5: Run CNV analysis (depends on filtering)",
    paste0("CNV_ID=$(bsub -q ", queue, " -w \"done($FILTER_ID)\" < ", basename(script_paths$cnv), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"CNV analysis job submitted with ID: $CNV_ID\"",
    "",
    "# Step 6: Run visualization (depends on dim reduction and CNV)",
    paste0("VIZ_ID=$(bsub -q ", queue, " -w \"done($DIM_ID) && done($CNV_ID)\" < ", basename(script_paths$visualization), " | grep -o \"<[0-9]*>\" | tr -d \"<>\" )"),
    "echo \"Visualization job submitted with ID: $VIZ_ID\"",
    "",
    "echo \"All jobs submitted successfully\""
  )
  
  # Write workflow script
  workflow_file <- file.path(scripts_dir, "submit_lsf_workflow.sh")
  writeLines(lsf_workflow, workflow_file)
  
  # Make executable
  if (.Platform$OS.type == "unix") {
    system(paste("chmod +x", workflow_file))
  }
  
  return(workflow_file)
}

#' Add PBS/Torque specific headers to a script
#'
#' @param script_content Script content
#' @param step Pipeline step
#' @param mem_req Memory requirement
#' @param time_req Time requirement
#' @param threads Number of threads
#' @return Modified script content
add_pbs_headers <- function(script_content, step, mem_req, time_req, threads) {
  # Replace SLURM headers with PBS headers
  headers <- c(
    paste0("#PBS -N meqtrack_", step),
    "#PBS -l nodes=1:ppn=1",
    paste0("#PBS -l mem=", mem_req),
    paste0("#PBS -l walltime=", time_req),
    paste0("#PBS -l ncpus=", threads),
    "#PBS -j oe",
    paste0("#PBS -o meqtrack_", step, ".out")
  )
  
  # Find where to insert headers (after shebang)
  insert_pos <- 2
  
  # Insert headers
  script_content <- c(
    script_content[1],
    headers,
    script_content[(insert_pos+1):length(script_content)]
  )
  
  return(script_content)
}

#' Add LSF (bsub) specific headers to a script
#'
#' @param script_content Script content
#' @param step Pipeline step
#' @param mem_req Memory requirement
#' @param time_req Time requirement
#' @param threads Number of threads
#' @return Modified script content
add_lsf_headers <- function(script_content, step, mem_req, time_req, threads) {
  # Convert memory requirement to LSF format (e.g., 16g to 16GB)
  mem_req_lsf <- gsub("g$", "GB", mem_req)
  
  # Convert time requirement to minutes for LSF (HH:MM:SS to minutes)
  time_parts <- as.numeric(strsplit(time_req, ":")[[1]])
  time_minutes <- time_parts[1] * 60 + time_parts[2]
  
  # Replace SLURM headers with LSF headers
  headers <- c(
    paste0("#BSUB -J meqtrack_", step),
    paste0("#BSUB -n ", threads),
    paste0("#BSUB -R \"rusage[mem=", mem_req_lsf, "]\""),
    paste0("#BSUB -W ", time_minutes),
    "#BSUB -cwd $(pwd)",
    paste0("#BSUB -o meqtrack_", step, ".%J.out"),
    paste0("#BSUB -e meqtrack_", step, ".%J.err")
  )
  
  # Find where to insert headers (after shebang)
  insert_pos <- 2
  
  # Insert headers
  script_content <- c(
    script_content[1],
    headers,
    script_content[(insert_pos+1):length(script_content)]
  )
  
  return(script_content)
}
