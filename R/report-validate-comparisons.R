# Check retained pairwise aggregates against their models and declared inference.
.ig_report_validate_comparisons <- function(r, s) {
  bad <- function(message) .ig_abort("invalid_input", paste("Invalid classifier comparison:", message))
  scalar <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  near <- function(a, b) scalar(a) && scalar(b) && abs(a - b) <= 1e-10 * max(1, abs(a), abs(b))
  cell <- function(x, j) if (is.null(x)) NULL else if (is.list(x)) x[[j]] else x[j]
  x <- r$comparisons
  ids <- vapply(r$models, `[[`, "", "model_id")
  expected <- if (s$input_kind == "scored_data") choose(length(ids), 2) else 0
  if (expected == 0) {
    if (!is.null(x)) bad("comparisons are not available for this input")
    return(invisible(TRUE))
  }
  types <- .ig_cls_table_columns$comparisons
  if (!is.list(x) || anyDuplicated(names(x)) || !setequal(names(x), names(types))) bad("missing or unexpected columns")
  if (any(vapply(x, function(v) max(1L, length(v)), 0L) != expected)) bad("incomplete or duplicated pair inventory")
  seen <- character()
  for (j in seq_len(expected)) {
    row <- lapply(x, cell, j = j)
    for (k in names(types)) {
      v <- row[[k]]
      if (!is.null(v) && !(if (types[[k]] == "character") .ig_is_string(v) else scalar(v))) bad("invalid cell type")
    }
    a <- row$model_a; b <- row$model_b
    if (!.ig_is_string(a) || !.ig_is_string(b) || !a %in% ids || !b %in% ids || a == b) bad("unknown or repeated model")
    key <- paste(sort(match(c(a, b), ids)), collapse = ":")
    if (key %in% seen) bad("duplicate model pair")
    seen <- c(seen, key)
    if (!identical(row$method, "delong_paired")) bad("unknown paired method")
    if (!scalar(row$level) || row$level != s$uncertainty$level) bad("confidence level differs from specification")
    if (!.ig_is_string(row$inference_status) || !row$inference_status %in% c("ok", "undefined", "unsupported")) bad("unknown inference status")
    for (side in c("a", "b")) {
      m <- r$models[[match(row[[paste0("model_", side)]], ids)]]
      ix <- which(unlist(m$summaries$metric_id, use.names = FALSE) == "roc_auc")
      if (length(ix) != 1L) bad("model AUC summary is missing or duplicated")
      value <- cell(m$summaries$value, ix)
      actual <- row[[paste0("auc_", side)]]
      if (is.null(value)) {
        if (!is.null(actual)) bad("unavailable AUC must remain unavailable")
      } else if (!near(actual, value)) bad("AUC differs from model summary")
    }
    if (scalar(row$auc_a) && scalar(row$auc_b)) {
      if (!near(row$difference, row$auc_a - row$auc_b)) bad("difference disagrees with AUCs")
    } else if (!is.null(row$difference)) bad("difference requires both AUCs")
    infer <- c("se", "lower", "upper", "z", "p_value")
    if (row$inference_status == "ok") {
      if (!all(vapply(row[infer], scalar, FALSE)) || row$se <= 0 || row$p_value < 0 || row$p_value > 1 || row$lower > row$upper) bad("invalid inference values")
      zcrit <- stats::qnorm((1 - row$level) / 2, lower.tail = FALSE)
      if (!near(row$z, row$difference / row$se) || !near(row$p_value, 2 * stats::pnorm(-abs(row$z))) ||
          !near(row$lower, row$difference - zcrit * row$se) || !near(row$upper, row$difference + zcrit * row$se)) bad("inference fields disagree")
    } else {
      if (any(!vapply(row[infer], is.null, FALSE)) || !.ig_is_string(row$reason_code)) bad("unavailable inference needs null values and a reason")
    }
  }
  invisible(TRUE)
}
