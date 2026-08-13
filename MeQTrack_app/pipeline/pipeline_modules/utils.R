# Utility functions for methylation array analysis pipeline

#' Load R package dynamically
#'
#' @param package_name Name of the package to load
#' @param quiet Whether to suppress messages
#' @return TRUE if package was loaded successfully, FALSE otherwise
load_package <- function(package_name, quiet = FALSE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!quiet) {
      message(paste("Package", package_name, "is not installed. Trying to install..."))
    }
    
    # Try to install package
    install_result <- try(install.packages(package_name, repos = "https://cran.r-project.org"), 
                        silent = quiet)
    
    if (inherits(install_result, "try-error")) {
      if (!quiet) {
        message(paste("Failed to install package", package_name))
      }
      return(FALSE)
    }
  }
  
  # Try to load package
  load_result <- try(library(package_name, character.only = TRUE), silent = quiet)
  
  if (inherits(load_result, "try-error")) {
    if (!quiet) {
      message(paste("Failed to load package", package_name))
    }
    return(FALSE)
  }
  
  return(TRUE)
}

#' Load Bioconductor package dynamically
#'
#' @param package_name Name of the Bioconductor package to load
#' @param quiet Whether to suppress messages
#' @return TRUE if package was loaded successfully, FALSE otherwise
load_bioc_package <- function(package_name, quiet = FALSE) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!quiet) {
      message(paste("Bioconductor package", package_name, "is not installed. Trying to install..."))
    }
    
    # Check if BiocManager is installed
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install_result <- try(install.packages("BiocManager", repos = "https://cran.r-project.org"),
                          silent = quiet)
      
      if (inherits(install_result, "try-error")) {
        if (!quiet) {
          message("Failed to install BiocManager")
        }
        return(FALSE)
      }
    }
    
    # Try to install Bioconductor package
    install_result <- try(BiocManager::install(package_name, update = FALSE, ask = FALSE),
                        silent = quiet)
    
    if (inherits(install_result, "try-error")) {
      if (!quiet) {
        message(paste("Failed to install Bioconductor package", package_name))
      }
      return(FALSE)
    }
  }
  
  # Try to load package
  load_result <- try(library(package_name, character.only = TRUE), silent = quiet)
  
  if (inherits(load_result, "try-error")) {
    if (!quiet) {
      message(paste("Failed to load Bioconductor package", package_name))
    }
    return(FALSE)
  }
  
  return(TRUE)
}

#' Check for required dependencies
#'
#' @param dependencies List of required package names
#' @param bioc_dependencies List of required Bioconductor package names
#' @param quiet Whether to suppress messages
#' @return TRUE if all dependencies are available, FALSE otherwise
check_dependencies <- function(dependencies = NULL, bioc_dependencies = NULL, quiet = FALSE) {
  all_ok <- TRUE
  
  # Check R packages
  if (!is.null(dependencies)) {
    for (pkg in dependencies) {
      if (!load_package(pkg, quiet)) {
        all_ok <- FALSE
        if (!quiet) {
          message(paste("Missing required dependency:", pkg))
        }
      }
    }
  }
  
  # Check Bioconductor packages
  if (!is.null(bioc_dependencies)) {
    for (pkg in bioc_dependencies) {
      if (!load_bioc_package(pkg, quiet)) {
        all_ok <- FALSE
        if (!quiet) {
          message(paste("Missing required Bioconductor dependency:", pkg))
        }
      }
    }
  }
  
  return(all_ok)
}

#' Install package from source file
#'
#' @param package_name Name of the package
#' @param source_path Path to the source tarball (.tar.gz)
#' @param quiet Whether to suppress messages
#' @return TRUE if package was installed successfully, FALSE otherwise
install_source_package <- function(package_name, source_path, quiet = FALSE) {
  # Check if package is already installed
  if (requireNamespace(package_name, quietly = TRUE)) {
    if (!quiet) {
      message(paste("Package", package_name, "is already installed"))
    }
    return(TRUE)
  }
  
  # Check if source file exists
  if (!file.exists(source_path)) {
    if (!quiet) {
      message(paste("Source file not found:", source_path))
    }
    return(FALSE)
  }
  
  if (!quiet) {
    message(paste("Installing", package_name, "from source:", source_path))
  }
  
  # Install package from source
  install_result <- try(
    install.packages(source_path, repos = NULL, type = "source", quiet = quiet),
    silent = quiet
  )
  
  if (inherits(install_result, "try-error")) {
    if (!quiet) {
      message(paste("Failed to install", package_name, "from source"))
    }
    return(FALSE)
  }
  
  # Verify installation
  if (requireNamespace(package_name, quietly = TRUE)) {
    if (!quiet) {
      message(paste("Successfully installed", package_name, "from source"))
    }
    return(TRUE)
  } else {
    if (!quiet) {
      message(paste("Package", package_name, "was installed but cannot be loaded"))
    }
    return(FALSE)
  }
}

#' Install all required dependencies for the pipeline
#'
#' @param quiet Whether to suppress messages
#' @return TRUE if all dependencies were installed successfully, FALSE otherwise
install_dependencies <- function(quiet = FALSE) {
  # CRAN packages
  cran_packages <- c(
    "optparse", "data.table", "ggplot2", "plotly", "Rtsne", "umap",
    "dendextend", "circlize", "htmlwidgets", "rmarkdown", "knitr",
    "DT", "parallel", "yaml", "ggrepel"
  )
  
  # Bioconductor packages.
  # Note: the CNV backend uses `conumee2`, not the older `conumee` package.
  # conumee2 is handled by the project's setup.R (Bioc with GitHub fallback)
  # and is intentionally not listed here.
  bioc_packages <- c(
    "minfi", "limma", "missMethyl", "RColorBrewer", "matrixStats", "snifter",
    "DMRcate", "GenomicRanges", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylationEPICv2manifest",
    "IlluminaHumanMethylationEPICv2anno.20a1.hg38", "sesame", "Gviz"
  )
  
  # Install dependencies
  cran_ok <- check_dependencies(cran_packages, NULL, quiet)
  bioc_ok <- check_dependencies(NULL, bioc_packages, quiet)
  
  # Install yamapData from source
  yamap_ok <- TRUE
  yamap_source_path <- file.path("data", "yamapData_0.0.3.tar.gz")
  
  # Try different possible locations for the source file
  possible_paths <- c(
    yamap_source_path,
    file.path(".", "data", "yamapData_0.0.3.tar.gz"),
    file.path("..", "data", "yamapData_0.0.3.tar.gz"),
    "yamapData_0.0.3.tar.gz"
  )
  
  # Check if yamapData is already installed
  if (!requireNamespace("yamapData", quietly = TRUE)) {
    if (!quiet) {
      message("Installing yamapData from source...")
    }
    
    # Find the source file
    source_found <- FALSE
    for (path in possible_paths) {
      if (file.exists(path)) {
        yamap_ok <- install_source_package("yamapData", path, quiet)
        source_found <- TRUE
        break
      }
    }
    
    if (!source_found) {
      if (!quiet) {
        message("yamapData source file not found in any of these locations:")
        for (path in possible_paths) {
          message(paste("  -", path))
        }
        message("Please ensure yamapData_0.0.3.tar.gz is in the data/ directory")
      }
      yamap_ok <- FALSE
    }
  } else {
    if (!quiet) {
      message("yamapData is already installed")
    }
  }
  
  return(cran_ok && bioc_ok && yamap_ok)
}

#' Install yamapData package specifically
#'
#' @param source_path Path to yamapData source file (optional)
#' @param quiet Whether to suppress messages
#' @return TRUE if yamapData was installed successfully, FALSE otherwise
install_yamapData <- function(source_path = NULL, quiet = FALSE) {
  # Check if yamapData is already installed
  if (requireNamespace("yamapData", quietly = TRUE)) {
    if (!quiet) {
      message("yamapData is already installed")
    }
    return(TRUE)
  }
  
  # Default source path if not provided
  if (is.null(source_path)) {
    possible_paths <- c(
      file.path("data", "yamapData_0.0.3.tar.gz"),
      file.path(".", "data", "yamapData_0.0.3.tar.gz"),
      file.path("..", "data", "yamapData_0.0.3.tar.gz"),
      "yamapData_0.0.3.tar.gz"
    )
    
    # Find the source file
    source_found <- FALSE
    for (path in possible_paths) {
      if (file.exists(path)) {
        source_path <- path
        source_found <- TRUE
        break
      }
    }
    
    if (!source_found) {
      if (!quiet) {
        message("yamapData source file not found in any of these locations:")
        for (path in possible_paths) {
          message(paste("  -", path))
        }
        message("Please ensure yamapData_0.0.3.tar.gz is in the data/ directory")
      }
      return(FALSE)
    }
  }
  
  # Install yamapData from source
  return(install_source_package("yamapData", source_path, quiet))
}

#' Check if file exists and is readable
#'
#' @param file_path Path to file
#' @return TRUE if file exists and is readable, FALSE otherwise
check_file <- function(file_path) {
  return(file.exists(file_path) && file.access(file_path, 4) == 0)
}

#' Generate timestamp string
#'
#' @param format Timestamp format (default: %Y%m%d_%H%M%S)
#' @return Timestamp string
get_timestamp <- function(format = "%Y%m%d_%H%M%S") {
  return(format(Sys.time(), format))
}
