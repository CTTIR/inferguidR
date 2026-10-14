# Mutation tests for advisor_validate_rules() and advisor_rules(), part 2:
# defects that only the exhaustive walk of every answer path can find, the
# loader's refusal of invalid copies, and restoring under another rules
# version. Part 1 (test-advisor-mutation.R) holds the static checks.
#
# Every case writes a deliberately broken copy of rules.json to a temporary
# file (the other inputs stay the bundled copies) and must be rejected with
# the check that targets the defect. Part 1 shows that an unmodified copy
# written the same way passes. The bundled files are never modified.

tmw_doc_cache <- new.env(parent = emptyenv())

tm_doc <- function() {
  if (is.null(tmw_doc_cache$rules)) {
    path <- system.file("advisor", "rules.json", package = "inferguidR")
    tmw_doc_cache$rules <- jsonlite::parse_json(
      paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n"),
      simplifyVector = FALSE)
  }
  tmw_doc_cache$rules
}

tm_write_json <- function(doc) {
  path <- tempfile("advisor-mutant-", fileext = ".json")
  txt <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  writeLines(enc2utf8(as.character(txt)), path, useBytes = TRUE)
  path
}

# Position of the one record whose `key` equals `id`.
tm_pos <- function(items, key, id) {
  i <- which(vapply(items, function(z) if (is.list(z) && is.character(z[[key]])) z[[key]] else "", "") == id)
  if (length(i) != 1L) stop(sprintf("test setup: %s '%s' must occur exactly once", key, id))
  i
}
tm_q <- function(d, id) tm_pos(d$questions, "question_id", id)
tm_t <- function(d, id) tm_pos(d$transitions, "transition_id", id)
tm_o <- function(d, id) tm_pos(d$outcomes, "outcome_id", id)

# Contiguous orders 1..n again for the transitions from one question.
tm_renumber <- function(d, from) {
  idx <- which(vapply(d$transitions, function(t) t$from, "") == from)
  ord <- order(vapply(d$transitions[idx], function(t) as.numeric(t$order), numeric(1)))
  for (k in seq_along(ord)) d$transitions[[idx[ord[k]]]]$order <- k
  d
}

# Validate a mutated copy (with the walk) and check that it is rejected with
# every expected check.
tm_expect_rejected <- function(case) {
  path <- tm_write_json(case$mutate(tm_doc()))
  on.exit(unlink(path), add = TRUE)
  v <- advisor_validate_rules(rules = path, walk = TRUE)
  expect_s3_class(v, "advisor_rules_validation")
  checks <- unique(v$problems$check)
  missing <- setdiff(case$expect, checks)
  expect(!isTRUE(v$valid) && !length(missing),
         sprintf("%s (%s) was not rejected as expected: missing check(s) [%s]; reported [%s]",
                 case$id, case$what, paste(missing, collapse = ", "), paste(checks, collapse = ", ")))
  invisible(v)
}

# Defects that only an exhaustive walk can find (blueprint mutation set of the
# 01A review plus cycles, applies_when and reachability).
tm_walk_cases <- list(
  list(id = "W01", what = "removed transition of q.records_per_individual (guard reached)",
       expect = "guard_reached", mutate = function(d) {
         d$transitions[[tm_t(d, "t.records_per_individual.05")]] <- NULL
         tm_renumber(d, "q.records_per_individual")
       }),
  list(id = "W02", what = "groups_over_time routed to out.paired.t",
       expect = c("design_consistency", "unreachable_outcome"), mutate = function(d) {
         d$transitions[[tm_t(d, "t.structure.06")]]$to <- "out.paired.t"
         d
       }),
  list(id = "W03", what = "paired binary routed to Boschloo",
       expect = "design_consistency", mutate = function(d) {
         d$transitions[[tm_t(d, "t.missingness.46")]]$to <- "out.two_group.boschloo_central"
         d
       }),
  list(id = "W04", what = "several records per individual routed to Pearson",
       expect = "design_consistency", mutate = function(d) {
         d$transitions[[tm_t(d, "t.records_per_individual.04")]]$to <- "out.assoc.pearson"
         d
       }),
  list(id = "W05", what = "Welch replaced by the pooled t-test when equal variances are not justified",
       expect = c("design_consistency", "unreachable_outcome"), mutate = function(d) {
         d$transitions[[tm_t(d, "t.missingness.30")]]$to <- "out.two_group.student"
         d
       }),
  list(id = "W06", what = "Boschloo replaced by Fisher for fixed group sizes",
       expect = "design_consistency", mutate = function(d) {
         d$transitions[[tm_t(d, "t.missingness.37")]]$to <- "out.two_group.fisher_2x2"
         d
       }),
  list(id = "W07", what = "transition back to a question already on the path (cycle)",
       expect = "cycle", mutate = function(d) {
         d$transitions[[tm_t(d, "t.missingness.30")]]$to <- "q.direction"
         d
       }),
  list(id = "W08", what = "question reached although its applies_when is false",
       expect = "applies_when", mutate = function(d) {
         d$questions[[tm_q(d, "q.direction")]]$applies_when <- list(q = "q.aim", "in" = list("describe"))
         d
       }),
  list(id = "W09", what = "outcome that no answer path reaches",
       expect = "unreachable_outcome", mutate = function(d) {
         o <- d$outcomes[[tm_o(d, "out.two_group.welch")]]
         o$outcome_id <- "out.two_group.welch_copy"
         d$outcomes[[length(d$outcomes) + 1L]] <- o
         d
       }),
  list(id = "W10", what = "question that no answer path reaches",
       expect = "unreachable_question", mutate = function(d) {
         q <- d$questions[[tm_q(d, "q.variance")]]
         q$question_id <- "q.orphan"
         d$questions[[length(d$questions) + 1L]] <- q
         g <- d$transitions[[tm_t(d, "t.variance.03")]]
         g$transition_id <- "t.orphan.01"
         g$from <- "q.orphan"
         g$order <- 1L
         d$transitions[[length(d$transitions) + 1L]] <- g
         d
       }),
  list(id = "W11", what = "transition shadowed by an earlier one (never fires)",
       expect = "unused_transition", mutate = function(d) {
         t <- d$transitions[[tm_t(d, "t.aim.01")]]
         t$transition_id <- "t.aim.shadowed"
         t$order <- 5L
         d$transitions[[tm_t(d, "t.aim.05")]]$order <- 6L
         d$transitions[[length(d$transitions) + 1L]] <- t
         d
       }),
  list(id = "W12", what = "documented answer-path count does not match the walk",
       expect = "summary_counts", mutate = function(d) {
         d$summary$draft_check$answer_paths_enumerated <- 5075L
         d
       })
)

test_that("walk-detected defects are rejected", {
  for (case in tm_walk_cases) tm_expect_rejected(case)
})

test_that("advisor_rules() refuses invalid copies and loads valid ones fully checked", {
  d <- tm_doc()
  bad <- d
  bad$transitions[[tm_t(bad, "t.missingness.46")]]$to <- "out.two_group.boschloo_central"
  miscounted <- d
  miscounted$summary$n_outcomes <- miscounted$summary$n_outcomes + 1L
  pb <- tm_write_json(bad)
  ps <- tm_write_json(miscounted)
  pc <- tm_write_json(d)
  on.exit(unlink(c(pb, ps, pc)), add = TRUE)
  e <- tryCatch(advisor_rules(rules = pb), inferguidR_error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_true(any(grepl("^\\[design_consistency\\]", unlist(e$status$details$problems))))
  # content defects are refused too, not only structural ones
  e <- tryCatch(advisor_rules(rules = ps), inferguidR_error = function(e) e)
  expect_true(any(grepl("^\\[summary_counts\\]", unlist(e$status$details$problems))))
  # a valid copy loads; its content hash differs from the bundled files
  rc <- advisor_rules(rules = pc)
  expect_s3_class(rc, "advisor_rules")
  expect_identical(rc$versions$rules_version, advisor_rules()$versions$rules_version)
  expect_false(identical(rc$versions$rules_sha256, advisor_rules()$versions$rules_sha256))
  # sessions are bound to the exact rule content they were created with
  s <- advisor_answer(advisor_session(rules = rc), "compare", rules = rc)
  expect_identical(s$state$pending_question_id, "q.stage")
  e <- tryCatch(advisor_answer(s, "analysis"), inferguidR_error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_match(conditionMessage(e), "content hash")
})

test_that("restoring under a different rules version requires an explicit migration", {
  d <- tm_doc()
  d$rules_version <- "2026-10-06.1"
  d$change_log[[length(d$change_log) + 1L]] <- list(
    rules_version = "2026-10-06.1", date = "2026-10-06",
    summary_en = "Test copy with a new version label.",
    migration_en = "No answers change; this entry exists only for the migration test.")
  p <- tm_write_json(d)
  on.exit(unlink(p), add = TRUE)
  newer <- advisor_rules(rules = p)
  expect_identical(newer$versions$rules_version, "2026-10-06.1")
  s <- advisor_replay(c(q.aim = "compare", q.stage = "analysis", q.outcome_type = "continuous"))
  js <- advisor_save(s)
  e <- tryCatch(advisor_restore(js, rules = newer), inferguidR_error = function(e) e)
  expect_s3_class(e, "inferguidR_migration_required")
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_identical(e$status$details$reason_code, "migration_required")
  expect_identical(e$status$details$saved_rules_version, "2026-10-05.1")
  expect_identical(e$status$details$current_rules_version, "2026-10-06.1")
  versions <- vapply(e$status$details$change_log, `[[`, "", "rules_version")
  expect_true("2026-10-06.1" %in% versions)
  # and the other way round: a state saved under the newer version is not re-interpreted
  s2 <- advisor_replay(c(q.aim = "compare", q.stage = "analysis"), rules = newer)
  e2 <- tryCatch(advisor_restore(advisor_save(s2, rules = newer)), inferguidR_error = function(e) e)
  expect_s3_class(e2, "inferguidR_migration_required")
  expect_identical(e2$status$details$saved_rules_version, "2026-10-06.1")
  # the same version restores exactly
  expect_identical(advisor_restore(advisor_save(s2, rules = newer), rules = newer), s2)
})
