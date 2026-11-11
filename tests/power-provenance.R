library(inferguidR)
# Every procedure has one provenance row: its original(s) and the supported
# combinations, or the reason it has the base engine only.
provenance <- utils::read.csv(system.file("validation", "power-provenance.csv", package = "inferguidR", mustWork = TRUE),
  colClasses = "character")
stopifnot(identical(names(provenance), c("method_id", "engines", "original", "supported", "reason_base_only")))
listing <- list_methods()
stopifnot(setequal(provenance$method_id, listing$id), !anyDuplicated(provenance$method_id),
  identical(provenance$engines[match(listing$id, provenance$method_id)], listing$engines))
reference <- provenance$engines == "base+reference"
stopifnot(all(nzchar(provenance$original[reference])), all(nzchar(provenance$supported[reference])),
  all(nzchar(provenance$reason_base_only[!reference])), !any(nzchar(provenance$reason_base_only[reference])))
cat("PASS power provenance of", nrow(provenance), "procedures,", sum(reference), "with a reference engine\n")
