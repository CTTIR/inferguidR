# Copyright (c) 2026 Raban Heller. MIT licence; see ../../LICENSE.md.

write_numerical_inventory <- function(package_path, cases, output) {
  registry <- list.files(file.path(package_path, "registry"), "[.]json$", full.names = TRUE)
  methods <- unlist(lapply(registry, function(f) jsonlite::read_json(f, simplifyVector = FALSE)$methods), recursive = FALSE)
  stopifnot(length(methods) == 50L, !anyDuplicated(vapply(methods, `[[`, "", "method_id")))
  modes <- options <- combinations <- list()
  status <- function(rows) {
    if (!nrow(rows)) return("uncovered")
    if (any(!rows$passed)) return("sample_failed")
    if (all(rows$evidence_kind == "guard_contract")) return("guard_checked_no_numeric_solution")
    "sample_checked_not_certified"
  }
  encode <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null"))
  applicable <- function(rule, values, mode) {
    if (is.null(rule)) return(TRUE)
    if (!is.null(rule$condition)) return(applicable(rule$condition, values, mode))
    if (!is.null(rule$all)) return(all(vapply(rule$all, applicable, FALSE, values = values, mode = mode)))
    if (!is.null(rule$any)) return(any(vapply(rule$any, applicable, FALSE, values = values, mode = mode)))
    if (!is.null(rule$option) && !is.null(rule[["in"]])) return(values[[rule$option]] %in% unlist(rule[["in"]]))
    if (!is.null(rule$mode)) return(mode %in% unlist(rule$mode))
    # Input-dependent and unknown predicates cannot establish evidence.
    FALSE
  }
  for (method in methods) for (mode in unlist(method$modes)) {
    id <- method$method_id
    evidence <- cases[cases$method == id & cases$mode == mode, , drop = FALSE]
    modes[[length(modes) + 1L]] <- data.frame(method = id, mode = mode, evidence_status = status(evidence),
      case_count = nrow(evidence), evidence_cases = paste(evidence$case, collapse = ";"))
    declared <- method$options
    if (!length(declared)) declared <- list(list(option_id = "(none)", value_type = "none"))
    defaults <- stats::setNames(lapply(declared, `[[`, "default"), vapply(declared, `[[`, "", "option_id"))
    resolve <- function(supplied) {
      values <- defaults
      for (name in names(supplied)) values[[name]] <- supplied[[name]]
      values
    }
    case_options <- if (nrow(evidence)) lapply(seq_len(nrow(evidence)), function(i) resolve(jsonlite::fromJSON(evidence$options_json[i], simplifyVector = FALSE))) else list()
    for (option in declared) {
      values <- vapply(option$values, function(v) as.character(v$code), "")
      if (!length(values)) values <- "(continuous or metadata; see declaration_json)"
      for (value in values) {
        matches <- if (nrow(evidence)) vapply(seq_len(nrow(evidence)), function(i) {
          supplied <- case_options[[i]][[option$option_id]]
          applicable(option$applies_when, case_options[[i]], mode) &&
            length(supplied) == 1L && identical(as.character(supplied), value)
        }, FALSE) else logical()
        ev <- evidence[matches, , drop = FALSE]
        options[[length(options) + 1L]] <- data.frame(method = id, mode = mode, option = option$option_id,
          value = value, default_json = encode(option$default), applicability_json = encode(option$applies_when),
          declaration_json = encode(option), evidence_status = status(ev), case_count = nrow(ev),
          evidence_cases = paste(ev$case, collapse = ";"))
      }
    }
    # Enumerate the entire finite declaration Cartesian product. This includes
    # inactive and input-constrained combinations; applicability is explicit.
    # Continuous domains stay as domain markers, never falsely discretised.
    choices <- stats::setNames(lapply(declared, function(o) {
      values <- vapply(o$values, function(v) as.character(v$code), "")
      if (length(values)) values else "(domain)"
    }), vapply(declared, `[[`, "", "option_id"))
    grid <- expand.grid(choices, stringsAsFactors = FALSE)
    active_key <- function(values) {
      keep <- vapply(declared, function(o) applicable(o$applies_when, values, mode), FALSE)
      encode(values[names(choices)[keep]])
    }
    evidence_keys <- vapply(case_options, active_key, "")
    for (i in seq_len(nrow(grid))) {
      values <- as.list(grid[i, , drop = FALSE])
      key <- active_key(values)
      ev <- evidence[evidence_keys == key, , drop = FALSE]
      combinations[[length(combinations) + 1L]] <- data.frame(method = id, mode = mode,
        options_json = encode(values), active_options_json = key,
        applicability_status = "registry_candidate_not_input_validated",
        evidence_status = status(ev), evidence_cases = paste(ev$case, collapse = ";"))
    }
  }
  mode_rows <- do.call(rbind, modes); option_rows <- do.call(rbind, options); combination_rows <- do.call(rbind, combinations)
  utils::write.csv(mode_rows, file.path(output, "method-mode-inventory.csv"), row.names = FALSE)
  utils::write.csv(option_rows, file.path(output, "method-mode-option-inventory.csv"), row.names = FALSE)
  utils::write.csv(combination_rows, file.path(output, "method-mode-option-combinations.csv"), row.names = FALSE)
  jsonlite::write_json(list(schema_version = 1, scope = "complete declarations, marginal values and finite Cartesian option candidates; input feasibility is not inferred",
    evidence_scope = "this run only; historical or check-public evidence is not silently inferred",
    methods = methods), file.path(output, "registry-declarations.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  list(methods = length(methods), method_modes = nrow(mode_rows), method_mode_option_values = nrow(option_rows),
    finite_option_candidates = nrow(combination_rows),
    uncovered_method_modes = sum(mode_rows$evidence_status == "uncovered"),
    uncovered_option_values = sum(option_rows$evidence_status == "uncovered"))
}
