source("data-raw/causal-sem/check-scales.R")

check_paper_measurement <- function(source_file = Sys.getenv("DCVAR_EFFECTLITER_SOURCE", "")) {
  if (nzchar(source_file)) {
    environment <- new.env(parent = baseenv())
    environment$na.omit <- stats::na.omit
    sys.source(source_file, envir = environment)
    generator <- environment$generateMeasurementModel
    provenance <- list(path = normalizePath(source_file), md5 = tools::md5sum(source_file))
  } else if (requireNamespace("EffectLiteR", quietly = TRUE)) {
    generator <- EffectLiteR::generateMeasurementModel
    provenance <- list(package = "EffectLiteR",
                         version = as.character(utils::packageVersion("EffectLiteR")))
  } else {
    stop("Install EffectLiteR or set DCVAR_EFFECTLITER_SOURCE to its measurement source file.")
  }
  data <- data.frame(z1 = ordered(1:3), z2 = ordered(1:3), z3 = ordered(1:3))
  syntax <- generator(names = "L", indicators = list(L = names(data)), ncells = 2,
                       model = "tau-cong-categorical", data = data)
  stopifnot(grepl("L =~ c(1,1)*z1", syntax, fixed = TRUE),
            grepl("z1 | c(0,0)*t1", syntax, fixed = TRUE))
  for (item in names(data)) {
    stopifnot(grepl(paste0(item, " ~*~ c(1,NA)*", item), syntax, fixed = TRUE),
              grepl(paste0(item, " ~ c(0,0)*1"), syntax, fixed = TRUE))
  }
  list(passed = TRUE, syntax = syntax, provenance = provenance, scales = check_scales())
}

if (sys.nframe() == 0L) {
  result <- check_paper_measurement()
  print(result[c("passed", "provenance")])
  cat(result$syntax, "\n")
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)) saveRDS(result, args[1])
}
