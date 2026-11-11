# Engine registry for power procedures.
#
# A procedure is computed by inferguidR's own forward model ("base") or by the
# forward model of an established R package that computes the same quantity
# ("reference"). Entries are explicit and keyed by method id; an entry lists
# option branches, and the first branch whose `when()` accepts the options and
# inputs is used. A branch names the original (package, functions, minimum
# version, call), supplies `power(p)` with the base forward model's interface,
# and states whether the generic solvers may search over it (`generic_search`).
# Originals are optional (Suggests); a missing, older or inapplicable original
# gives an "unsupported" result and never falls back to the base engine.
# Files named engine-reference-*.R are collated before this file.
.ig_engine_info <- c(.ig_engine_info_t(), .ig_engine_info_f(), .ig_engine_info_z())

.ig_engines <- c("base", "reference")

# Resolve the reference engine of a specification: the record kept in the
# result's provenance, or an ig_status explaining why it is unavailable.
.ig_engine_resolve <- function(spec, registry = .ig_engine_info) {
  id <- spec$method_id
  entry <- registry[[id]]
  if (is.null(entry))
    return(ig_status("unsupported", sprintf("no reference engine is registered for %s; use engine = \"base\"", id)))
  branch <- NULL
  for (candidate in entry$branches) {
    if (isTRUE(candidate$when(spec$options, spec$inputs))) {
      branch <- candidate
      break
    }
  }
  if (is.null(branch))
    return(ig_status("unsupported", sprintf(
      "the original for %s does not provide this option combination; use engine = \"base\"", id)))
  required <- sprintf("%s (>= %s)", branch$package, branch$min_version)
  if (!requireNamespace(branch$package, quietly = TRUE))
    return(ig_status("unsupported", sprintf(
      "engine = \"reference\" for %s requires the optional package %s, which is not installed", id, required)))
  version <- utils::packageVersion(branch$package)
  if (version < branch$min_version)
    return(ig_status("unsupported", sprintf(
      "engine = \"reference\" for %s requires %s; version %s is installed", id, required, version)))
  absent <- branch$functions[!vapply(branch$functions, exists, NA, envir = asNamespace(branch$package),
    inherits = FALSE)]
  if (length(absent))
    return(ig_status("unsupported", sprintf(
      "engine = \"reference\" for %s requires %s::%s, which version %s does not provide", id, branch$package,
      absent[[1L]], version)))
  list(engine = "reference", package = branch$package, version = as.character(version), call = branch$call,
    branch = branch)
}

# Signal that a design lies outside the original's domain (for example fewer
# observations than it accepts); the runner reports it as unsupported.
.ig_ref_outside <- function(message)
  stop(structure(class = c("inferguidR_reference_domain", "error", "condition"),
    list(message = message, call = NULL)))

# The reference model of a method: the original's forward model replaces the
# base model, the base search hooks are dropped, and the smallest total sample
# size the original can evaluate becomes the lower bound of the search.
.ig_ref_method <- function(def, resolved, spec) {
  branch <- resolved$branch
  def$power <- branch$power
  # A procedure's design rule (for example equal group sizes) is kept by
  # rebuilding its a priori search around the original's forward model.
  def$a_priori_search <- if (is.function(branch$search)) branch$search(branch$power) else NULL
  def$sensitivity_search <- NULL
  def$engine <- list(package = resolved$package, fun = resolved$call, version = resolved$version,
                     engine = "reference")
  n_min <- if (is.function(branch$n_min)) branch$n_min(spec$inputs) else branch$n_min
  if (!is.null(n_min)) {
    at <- which(vapply(def$params, function(x) identical(x$name, "n_total"), NA))
    if (length(at)) def$params[[at]]$domain$lower <- max(def$params[[at]]$domain$lower, n_min)
  }
  def
}

.ig_engine_listing <- function(ids) {
  engines <- ifelse(ids %in% names(.ig_engine_info), "base+reference", "base")
  reference <- vapply(ids, function(id) {
    entry <- .ig_engine_info[[id]]
    if (is.null(entry)) return(NA_character_)
    paste(unique(vapply(entry$branches, function(b) paste0(b$package, "::", b$functions[[1L]]), "")),
      collapse = "; ")
  }, "", USE.NAMES = FALSE)
  list(engines = engines, reference = reference)
}
