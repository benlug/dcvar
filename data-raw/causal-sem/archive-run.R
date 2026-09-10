# Keep small reports in the repository. Keep full fits at the run location.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Supply the run directory and a new archive directory.")
input <- arguments[1]
output <- arguments[2]
if (dir.exists(output)) stop("Use a new archive directory.")
if (!file.exists(file.path(input, "run.rds"))) stop("The run manifest is missing.")
dir.create(output, recursive = TRUE)
files <- list.files(input, pattern = "[.]csv$", full.names = TRUE)
if (length(files) && !all(file.copy(files, output))) stop("A report file cannot be copied.")
source_directory <- file.path(input, "stan-sources")
if (dir.exists(source_directory) && !file.copy(source_directory, output, recursive = TRUE)) {
  stop("The Stan source snapshots cannot be copied.")
}
manifest <- readRDS(file.path(input, "run.rds"))
connection <- file(file.path(output, "manifest.txt"), open = "wt")
manifest_lines <- c("Causal SEM validation run", paste("Full results:", normalizePath(input)),
                    "", capture.output(dput(manifest$config)), "",
                    capture.output(print(manifest$session)))
writeLines(sub("[[:blank:]]+$", "", manifest_lines), connection)
close(connection)
cat("Reports saved to", output, "\n")
