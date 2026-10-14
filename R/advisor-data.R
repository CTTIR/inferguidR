# Advisor data: loading, indexing and versioning the rule set -----------------
#
# The advisor is driven by versioned data, never by code generated from it:
#
# * inst/advisor/rules.json             questions, transitions, outcomes
# * inst/advisor/analysis-variants.json analysis variants named by outcomes
# * inst/advisor/uzh-tree.json          UZH methods tree (source evidence)
# * inst/advisor/descrtab2-tree.json    DescrTab2 test-choice tree (source evidence)
# * inst/advisor/source-crosswalk.csv   crosswalk between both trees and variants
# * inst/metadata/learning-resources.json  CTTIR learning-resource registry
#
# The files in inst/advisor/ are byte-identical copies of the audited
# specifications in specs/advisor/. A loaded rule set is an `advisor_rules`
# object: the parsed documents plus lookup indexes, compiled conditions and
# version/hash provenance. The bundled set is parsed once per session and
# cached.

.ig_advisor_rules_schema <- "inferguidR.advisor-rules/0.1.0"
.ig_advisor_variants_schema <- "inferguidR.analysis-variants/0.1.0"

# Package-level cache for the bundled rule set.
.ig_advisor_cache <- new.env(parent = emptyenv())

# Installed path of a bundled advisor or metadata file ("" if absent).
.ig_advisor_bundled_path <- function(which) {
  switch(which,
         rules = system.file("advisor", "rules.json", package = "inferguidR"),
         variants = system.file("advisor", "analysis-variants.json", package = "inferguidR"),
         uzh_tree = system.file("advisor", "uzh-tree.json", package = "inferguidR"),
         descrtab2_tree = system.file("advisor", "descrtab2-tree.json", package = "inferguidR"),
         crosswalk = system.file("advisor", "source-crosswalk.csv", package = "inferguidR"),
         resources = system.file("metadata", "learning-resources.json", package = "inferguidR"),
         scenarios = system.file("advisor", "acceptance-scenarios.json", package = "inferguidR"),
         stop("unknown advisor file: ", which, call. = FALSE))
}

# Text of a local file (UTF-8) and the SHA-256 of that text with line endings
# normalised to LF, so that the content hash does not depend on the platform
# a copy was checked out on.
.ig_advisor_read_file <- function(path, what) {
  if (!.ig_is_string(path) || grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)) {
    .ig_abort("invalid_input", sprintf("%s must be the path of a local file", what))
  }
  txt <- .ig_advisor_read_local(path, what)
  txt <- gsub("\r\n?", "\n", txt, useBytes = TRUE)
  list(text = txt, sha256 = .ig_sha256(txt))
}

# Whole text of a local regular file (UTF-8). A missing file, a directory, a
# file over the size limit or a read failure (including R warnings such as
# permission problems) is a typed invalid_input condition, never an untyped
# error.
.ig_advisor_read_local <- function(path, what) {
  if (!file.exists(path)) {
    .ig_abort("invalid_input", sprintf("%s not found: %s", what, path))
  }
  if (dir.exists(path)) {
    .ig_abort("invalid_input", sprintf("%s is a directory, not a file: %s", what, path))
  }
  size <- file.size(path)
  if (is.na(size) || size > .ig_json_max_bytes) {
    .ig_abort("invalid_input", sprintf("%s is larger than %d bytes", what, .ig_json_max_bytes))
  }
  fail <- function(e) {
    .ig_abort("invalid_input", sprintf("%s could not be read (%s): %s", what, conditionMessage(e), path))
  }
  txt <- tryCatch(readChar(path, size, useBytes = TRUE), error = fail, warning = fail)
  if (!length(txt)) txt <- "" # empty file
  Encoding(txt) <- "UTF-8"
  txt
}

# Read one JSON document safely (no evaluation, duplicate keys rejected).
.ig_advisor_read_json <- function(path, what) {
  f <- .ig_advisor_read_file(path, what)
  structure(list(doc = .ig_json_parse(f$text)), sha256 = f$sha256)
}

.ig_advisor_read_crosswalk <- function(path) {
  f <- .ig_advisor_read_file(path, "source crosswalk")
  cw <- utils::read.csv(text = f$text, colClasses = "character", na.strings = character(),
                        check.names = FALSE, encoding = "UTF-8")
  structure(list(doc = cw), sha256 = f$sha256)
}

#' Load an advisor rule set
#'
#' Reads the versioned advisor data (rule graph, analysis variants, both
#' source trees, the source crosswalk and the learning-resource registry),
#' checks it and compiles every condition with an allowlisted interpreter.
#' Conditions are JSON data (`q`/`in`, `all`, `any`, `not`, `always`);
#' expression strings are never parsed or evaluated.
#'
#' Without arguments the rule set bundled with the package is returned (it is
#' read once and cached). Its files are byte-identical copies of the audited
#' specification and pass [advisor_validate_rules()] in the package tests, so
#' at load time only their structure is checked. Paths can be given to load a
#' different copy, for example a revised draft; such a copy is always fully
#' validated first (including the exhaustive walk, about two seconds) and
#' refused if any check fails. Use [advisor_validate_rules()] to inspect the
#' problems of a draft without loading it.
#'
#' @param rules,variants,uzh_tree,descrtab2_tree,crosswalk,resources Paths of
#'   local files, or `NULL` for the bundled copy: the rule graph
#'   (`rules.json`), analysis variants, UZH tree, DescrTab2 tree, source
#'   crosswalk (CSV) and learning-resource registry.
#' @return An object of class `advisor_rules`. Its `versions` element holds
#'   the rules id and version, the source-tree, variant and resource versions
#'   and SHA-256 content hashes of every input. A rule set that fails a check
#'   raises a condition of class `inferguidR_invalid_input` whose
#'   `$status$details$problems` lists every problem as `"[check] message"`.
#' @examples
#' rs <- advisor_rules()
#' rs$versions$rules_version
#' @seealso [advisor_validate_rules()], [advisor_session()]
#' @family advisor
#' @export
advisor_rules <- function(rules = NULL, variants = NULL, uzh_tree = NULL,
                          descrtab2_tree = NULL, crosswalk = NULL, resources = NULL) {
  paths <- list(rules = rules, variants = variants, uzh_tree = uzh_tree,
                descrtab2_tree = descrtab2_tree, crosswalk = crosswalk,
                resources = resources)
  bundled <- all(vapply(paths, is.null, logical(1)))
  if (bundled && !is.null(.ig_advisor_cache$bundled)) return(.ig_advisor_cache$bundled)
  if (!bundled) {
    v <- advisor_validate_rules(rules, variants, uzh_tree, descrtab2_tree, crosswalk, resources,
                                walk = TRUE)
    if (!isTRUE(v$valid)) {
      .ig_abort_invalid(sprintf("[%s] %s", v$problems$check, v$problems$message),
                        "the advisor rule set failed validation")
    }
  }
  rs <- .ig_advisor_load_checked(paths)
  problems <- .ig_advisor_guarded(.ig_advisor_structure_problems(rs), "schema",
                                  "the structure could not be checked")
  if (nrow(problems)) {
    .ig_abort_invalid(sprintf("[%s] %s", problems$check, problems$message),
                      "the advisor rule set is not structurally valid")
  }
  rs <- .ig_advisor_compile(rs)
  if (bundled) .ig_advisor_cache$bundled <- rs
  rs
}

# Evaluate a problem-table producing expression; an unexpected R error (a
# malformed document whose shape the checks did not anticipate) becomes one
# problem row of class `check` instead of escaping as an untyped error.
.ig_advisor_guarded <- function(expr, check, context) {
  tryCatch(expr, error = function(e) {
    .ig_advisor_problem_table(check, sprintf("%s: %s", context, conditionMessage(e)))
  })
}

# .ig_advisor_load() with every failure turned into a typed invalid_input
# (or unsupported, for a damaged installation) condition.
.ig_advisor_load_checked <- function(paths) {
  tryCatch(.ig_advisor_load(paths),
           inferguidR_error = function(e) stop(e),
           error = function(e) {
             .ig_abort("invalid_input", paste("the advisor data could not be read:", conditionMessage(e)))
           })
}

# Parse and index all inputs (no validation beyond what indexing needs).
.ig_advisor_load <- function(paths) {
  for (nm in names(paths)) {
    if (is.null(paths[[nm]])) {
      p <- .ig_advisor_bundled_path(nm)
      if (!nzchar(p)) {
        .ig_abort("unsupported", sprintf("bundled advisor file '%s' is missing from the installation", nm))
      }
      paths[[nm]] <- p
    }
  }
  read <- list(
    rules = .ig_advisor_read_json(paths$rules, "rules file"),
    variants = .ig_advisor_read_json(paths$variants, "analysis-variants file"),
    uzh_tree = .ig_advisor_read_json(paths$uzh_tree, "UZH tree file"),
    descrtab2_tree = .ig_advisor_read_json(paths$descrtab2_tree, "DescrTab2 tree file"),
    resources = .ig_advisor_read_json(paths$resources, "learning-resource registry"),
    crosswalk = .ig_advisor_read_crosswalk(paths$crosswalk)
  )
  docs <- lapply(read, `[[`, "doc")
  hashes <- lapply(read, attr, "sha256")
  .ig_advisor_index(docs, hashes)
}

# Build lookup tables. Malformed entries are skipped here and reported by the
# structural checks. `hashes` are the SHA-256 content hashes of the inputs.
.ig_advisor_index <- function(docs, hashes) {
  r <- docs$rules
  by_id <- function(items, key) {
    items <- Filter(function(z) is.list(z) && .ig_is_string(z[[key]]), items %||% list())
    stats::setNames(items, vapply(items, `[[`, "", key))
  }
  questions <- by_id(r$questions, "question_id")
  outcomes <- by_id(r$outcomes, "outcome_id")
  transitions <- Filter(function(t) is.list(t) && .ig_is_string(t$from), r$transitions %||% list())
  trans_from <- split(transitions, vapply(transitions, `[[`, "", "from"))
  trans_from <- lapply(trans_from, function(ts) {
    ord <- vapply(ts, function(t) if (is.numeric(t$order) && length(t$order) == 1L) as.numeric(t$order) else NA_real_,
                  numeric(1))
    ts[order(ord, na.last = TRUE)]
  })
  variants <- by_id(docs$variants$variants, "analysis_variant_id")
  refs <- by_id(r$references, "ref_id")
  uzh_nodes <- vapply(Filter(function(n) .ig_is_string(n$node_id), docs$uzh_tree$nodes %||% list()),
                      `[[`, "", "node_id")
  dt2_nodes <- c(docs$descrtab2_tree$nodes %||% list(), docs$descrtab2_tree$code_only_nodes %||% list())
  dt2_nodes <- vapply(Filter(function(n) .ig_is_string(n$id), dt2_nodes), `[[`, "", "id")
  res <- Filter(function(z) is.list(z) && .ig_is_string(z$resource_id),
                if (is.null(names(docs$resources))) docs$resources else list())
  hashes <- list(
    rules = hashes$rules,
    analysis_variants = hashes$variants,
    uzh_tree = hashes$uzh_tree,
    descrtab2_tree = hashes$descrtab2_tree,
    source_crosswalk = hashes$crosswalk,
    learning_resources = hashes$resources
  )
  dt2 <- docs$descrtab2_tree
  versions <- list(
    rules_id = r$rules_id %||% NA_character_,
    rules_version = r$rules_version %||% NA_character_,
    rules_schema_version = r$schema_version %||% NA_character_,
    rules_status = r$status %||% NA_character_,
    rules_sha256 = hashes$rules,
    uzh_tree_version = docs$uzh_tree$tree_version %||% NA_character_,
    descrtab2_tree_version = if (.ig_is_string(dt2$audit_date)) {
      paste0("audit ", dt2$audit_date,
             if (.ig_is_string(dt2$version_code$audited_version))
               paste0(" (DescrTab2 ", dt2$version_code$audited_version, ")") else "")
    } else NA_character_,
    analysis_variants_version = docs$variants$variants_version %||% NA_character_,
    # The resource registry carries no version field; its content hash is its version.
    learning_resources_version = paste0("sha256:", hashes$learning_resources),
    data_sha256 = hashes
  )
  # Option lookup per question: option records by id.
  options <- lapply(questions, function(q) {
    opts <- Filter(function(o) is.list(o) && .ig_is_string(o$option_id),
                   if (is.list(q$options)) q$options else list())
    stats::setNames(opts, vapply(opts, `[[`, "", "option_id"))
  })
  structure(list(
    docs = docs,
    questions = questions,
    options = options,
    outcomes = outcomes,
    transitions = trans_from,
    variants = variants,
    references = refs,
    uzh_nodes = unique(uzh_nodes),
    descrtab2_nodes = unique(dt2_nodes),
    resources = res,
    versions = versions
  ), class = "advisor_rules")
}

# Attach compiled conditions (functions of the effective-answer vector).
.ig_advisor_compile <- function(rs) {
  rs$compiled <- list(
    applies = lapply(rs$questions, function(q) .ig_advisor_compile_condition(q$applies_when)),
    transitions = lapply(rs$transitions, function(ts) {
      lapply(ts, function(t) .ig_advisor_compile_condition(t$when))
    }),
    modifiers = lapply(rs$outcomes, function(o) {
      lapply(o$status_modifiers %||% list(), function(m) .ig_advisor_compile_condition(m$when))
    })
  )
  rs
}

# Is a node id an outcome (rather than a question)?
.ig_advisor_is_outcome <- function(rs, id) !is.null(rs$outcomes[[id]])

#' @export
format.advisor_rules <- function(x, ...) {
  v <- x$versions
  c(sprintf("<advisor_rules> %s %s (%s)", v$rules_id, v$rules_version, v$rules_status),
    sprintf("  questions: %d, transitions: %d, outcomes: %d, analysis variants: %d",
            length(x$questions), sum(lengths(x$transitions)), length(x$outcomes),
            length(x$variants)),
    sprintf("  UZH tree %s; DescrTab2 tree %s; variants %s", v$uzh_tree_version,
            v$descrtab2_tree_version, v$analysis_variants_version),
    sprintf("  learning resources: %d (%s)", length(x$resources),
            substr(v$learning_resources_version, 1, 23)))
}

#' @export
print.advisor_rules <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}
