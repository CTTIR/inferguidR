# Deterministic example models -------------------------------------------------------------
#
# Neutral, synthetic examples for documentation, tests and the application.
# No random numbers are drawn: scores are normal quantiles at ppoints(),
# rounded so that ties occur, and units are interleaved by a fixed
# permutation. No source data, model or population of any real project is
# reproduced, and no policy target is attached.

# Synthetic scored data: 60 events and 140 non-events with two probability
# scores (model_a, model_b) rounded to two decimals.
.ig_cls_synthetic_data <- function() {
  n1 <- 60L
  n0 <- 140L
  z1 <- 1.1 + 1.2 * stats::qnorm(stats::ppoints(n1))
  z0 <- stats::qnorm(stats::ppoints(n0))
  # deterministic interleaving (fractional parts of an irrational multiple)
  ord1 <- order((seq_len(n1) * 0.6180339887) %% 1)
  ord0 <- order((seq_len(n0) * 0.4142135624) %% 1)
  z1 <- z1[ord1]
  z0 <- z0[ord0]
  latent <- c(z1, z0)
  y <- c(rep(1L, n1), rep(0L, n0))
  wobble <- 0.9 * sin(seq_along(latent) * 2.399963)
  d <- data.frame(
    id = sprintf("u%03d", seq_along(y)),
    outcome = ifelse(y == 1L, "event", "no_event"),
    model_a = round(stats::plogis(-1.2 + 1.1 * latent), 2),
    model_b = round(stats::plogis(-1.0 + 0.8 * latent + wobble), 2),
    stringsAsFactors = FALSE
  )
  # mix events and non-events in the row order
  perm <- order((seq_len(nrow(d)) * 0.7548776662) %% 1)
  d <- d[perm, ]
  rownames(d) <- NULL
  d$id <- sprintf("u%03d", seq_len(nrow(d)))
  d
}

#' Deterministic example classification analyses
#'
#' Neutral synthetic examples (no real data, models or policy targets):
#'
#' * `"scored_synthetic"`: 200 units (60 events) with two predicted
#'   probabilities `model_a` and `model_b` rounded to two decimals (ties
#'   occur); `part = "tuning"` and `part = "test"` give two disjoint halves
#'   for held-out evaluation.
#' * `"binormal_equal_variance"`: binormal model with AUC 0.8 and
#'   \eqn{\sigma = 1}, prevalence 0.1, N = 1000, cutoff at the Youden optimum
#'   \eqn{d/2}.
#' * `"binormal_unequal_variance"`: AUC 0.8, \eqn{\sigma = 1.5}, prevalence
#'   0.1, N = 1000, cutoff 1.
#' * `"counts_fixed"`: one fixed 2x2 table (TP 40, FN 10, FP 20, TN 130).
#'
#' @param name Example name.
#' @param part For `"scored_synthetic"`: `"all"` (default), `"tuning"` or
#'   `"test"`.
#' @param models For `"scored_synthetic"`: score columns to use.
#' @return A list with `spec` (a `classification_spec`), `data` (the data
#'   frame for scored data, otherwise `NULL`) and `description`.
#' @examples
#' ex <- classification_examples("binormal_unequal_variance")
#' evaluate_cutoff(ex$spec)
#' @family classification
#' @export
classification_examples <- function(name = c("scored_synthetic", "binormal_equal_variance",
                                             "binormal_unequal_variance", "counts_fixed"),
                                    part = c("all", "tuning", "test"), models = c("model_a", "model_b")) {
  name <- match.arg(name)
  part <- match.arg(part)
  if (name != "scored_synthetic" && part != "all") {
    .ig_abort("invalid_input", "`part` applies to the scored example only")
  }
  if (name == "scored_synthetic") {
    models <- match.arg(models, several.ok = TRUE)
    d <- .ig_cls_synthetic_data()
    if (part != "all") {
      keep <- if (part == "tuning") seq_len(nrow(d)) %% 2L == 1L else seq_len(nrow(d)) %% 2L == 0L
      d <- d[keep, ]
      rownames(d) <- NULL
    }
    meta <- lapply(models, function(m) {
      classifier_model(m, name = if (m == "model_a") "Example model A" else "Example model B",
                       version = "synthetic-1",
                       score_meaning = "predicted probability of the event (synthetic)",
                       unit = "unit")
    })
    spec <- classification_spec(d, truth = "outcome", positive_label = "event", negative_label = "no_event",
                                scores = models, models = meta, cutoff = 0.3,
                                score_scale = "probability", id = "id",
                                evaluation_split = if (part == "test") "held_out_test" else
                                  if (part == "tuning") "held_out_tuning" else "same_data")
    return(list(spec = spec, data = d, description = sprintf(
      "Synthetic scored data (%s part): %d units, %d events, %d model(s); probabilities rounded to 2 decimals.",
      part, nrow(d), sum(d$outcome == "event"), length(models))))
  }
  if (name == "counts_fixed") {
    spec <- classification_spec(counts = c(tp = 40, fn = 10, fp = 20, tn = 130), cutoff = 0.5,
                                rule = NULL)
    return(list(spec = spec, data = NULL, description = "One fixed 2x2 table: TP 40, FN 10, FP 20, TN 130."))
  }
  sigma <- if (name == "binormal_equal_variance") 1 else 1.5
  d <- stats::qnorm(0.8) * sqrt(1 + sigma^2)
  cutoff <- if (sigma == 1) d / 2 else 1
  spec <- classification_spec(theory = list(auc = 0.8, sigma = sigma, n_total = 1000),
                              cutoff = cutoff, prevalence = 0.1)
  list(spec = spec, data = NULL, description = sprintf(
    "Binormal model: AUC 0.8, sigma %s, prevalence 0.1 (assumed), N 1000, cutoff %s on the z scale of the non-events.",
    format(sigma), format(signif(cutoff, 6))))
}
