# Mutation tests for advisor_validate_rules(), part 1: defects found by the
# static checks, defective companion files and unreadable inputs. Part 2
# (test-advisor-mutation-walk.R) holds the defects that need the exhaustive
# walk, and the advisor_rules() loader and migration tests.
#
# Every case writes a deliberately broken copy of one advisor data file to a
# temporary file (the other inputs stay the bundled copies) and must be
# rejected with the check that targets the defect. A control copy written the
# same way passes, so a rejection is caused by the mutation and not by the
# rewrite itself. The bundled files are never modified.

tm_cache <- new.env(parent = emptyenv())

tm_doc <- function(file = "rules.json") {
  if (is.null(tm_cache[[file]])) {
    path <- system.file("advisor", file, package = "inferguidR")
    tm_cache[[file]] <- jsonlite::parse_json(
      paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n"),
      simplifyVector = FALSE)
  }
  tm_cache[[file]]
}

tm_crosswalk <- function() {
  utils::read.csv(system.file("advisor", "source-crosswalk.csv", package = "inferguidR"),
                  colClasses = "character", na.strings = character(), check.names = FALSE,
                  encoding = "UTF-8")
}

tm_write_json <- function(doc) {
  path <- tempfile("advisor-mutant-", fileext = ".json")
  txt <- jsonlite::toJSON(doc, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  writeLines(enc2utf8(as.character(txt)), path, useBytes = TRUE)
  path
}

tm_write_text <- function(txt, ext = ".json") {
  path <- tempfile("advisor-mutant-", fileext = ext)
  writeLines(enc2utf8(txt), path, useBytes = TRUE)
  path
}

tm_write_csv <- function(df) {
  path <- tempfile("advisor-mutant-", fileext = ".csv")
  utils::write.csv(df, path, row.names = FALSE, fileEncoding = "UTF-8")
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
tm_opt <- function(q, id) tm_pos(q$options, "option_id", id)

# Validate a mutated copy without the walk (static checks only) and check
# that it is rejected with every expected check. `arg` names the input that
# `path` replaces; without `path` the case mutates rules.json.
tm_expect_rejected <- function(case, arg = "rules", path = NULL) {
  if (is.null(path)) path <- tm_write_json(case$mutate(tm_doc()))
  on.exit(unlink(path), add = TRUE)
  args <- stats::setNames(list(path), arg)
  v <- do.call(advisor_validate_rules, c(args, list(walk = FALSE)))
  expect_s3_class(v, "advisor_rules_validation")
  checks <- unique(v$problems$check)
  missing <- setdiff(case$expect, checks)
  expect(!isTRUE(v$valid) && !length(missing),
         sprintf("%s (%s) was not rejected as expected: missing check(s) [%s]; reported [%s]",
                 case$id, case$what, paste(missing, collapse = ", "), paste(checks, collapse = ", ")))
  invisible(v)
}

test_that("control copies rewritten like the mutants pass validation", {
  p <- tm_write_json(tm_doc())
  on.exit(unlink(p), add = TRUE)
  v <- advisor_validate_rules(rules = p)
  if (nrow(v$problems)) print(v$problems)
  expect_true(v$valid)
  expect_identical(v$walk$answer_paths, 5076L)
  pv <- tm_write_json(tm_doc("analysis-variants.json"))
  pc <- tm_write_csv(tm_crosswalk())
  on.exit(unlink(c(pv, pc)), add = TRUE)
  expect_true(advisor_validate_rules(variants = pv, crosswalk = pc, walk = FALSE)$valid)
})

# Defects found by the static checks (no walk needed).
tm_static_cases <- list(
  list(id = "S01", what = "dont_know option removed", expect = "dont_know_option", mutate = function(d) {
    qi <- tm_q(d, "q.variance")
    d$questions[[qi]]$options[[tm_opt(d$questions[[qi]], "dont_know")]] <- NULL
    d
  }),
  list(id = "S02", what = "rank-test null assumption removed", expect = "required_assumptions",
       mutate = function(d) {
         oi <- tm_o(d, "out.two_group.wmw")
         d$outcomes[[oi]]$assumptions <- Filter(function(a) a$assumption_id != "a.rank_null_identical",
                                                d$outcomes[[oi]]$assumptions)
         d
       }),
  list(id = "S03", what = "dead-end outcome without a next action", expect = "next_actions",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.contradictory_design")]]$next_actions <- list()
         d
       }),
  list(id = "S04", what = "pre-test wording in an option label", expect = "forbidden_terms",
       mutate = function(d) {
         qi <- tm_q(d, "q.variance")
         oi <- tm_opt(d$questions[[qi]], "justified")
         d$questions[[qi]]$options[[oi]]$label_en <- "Equal variances (Levene test not significant)"
         d
       }),
  list(id = "S05", what = "Wilcoxon modifier pointing to an unlisted mixed-model alternative",
       expect = "modifiers", mutate = function(d) {
         m <- d$outcomes[[tm_o(d, "out.paired.t")]]$status_modifiers[[1]]
         d$outcomes[[tm_o(d, "out.paired.wilcoxon")]]$status_modifiers[[1]] <- m
         d
       }),
  list(id = "S06", what = "operator outside the allowlist", expect = "condition", mutate = function(d) {
    d$transitions[[tm_t(d, "t.aim.01")]]$when <- list(eval = "TRUE")
    d
  }),
  list(id = "S07", what = "condition with an extra key", expect = "condition", mutate = function(d) {
    d$transitions[[tm_t(d, "t.aim.01")]]$when <- list(q = "q.aim", "in" = list("describe"), call = "x")
    d
  }),
  list(id = "S08", what = "condition naming an unknown option", expect = "condition", mutate = function(d) {
    ti <- tm_t(d, "t.missingness.30")
    parts <- d$transitions[[ti]]$when$all
    k <- which(vapply(parts, function(p) identical(p$q, "q.variance"), logical(1)))
    d$transitions[[ti]]$when$all[[k]][["in"]] <- list("unequal")
    d
  }),
  list(id = "S09", what = "applies_when naming an unknown question", expect = "condition",
       mutate = function(d) {
         d$questions[[tm_q(d, "q.direction")]]$applies_when <- list(q = "q.nonexistent", "in" = list("yes"))
         d
       }),
  list(id = "S10", what = "transition to a target that does not exist", expect = "transitions",
       mutate = function(d) {
         d$transitions[[tm_t(d, "t.structure.02")]]$to <- "q.nowhere"
         d
       }),
  list(id = "S11", what = "final guard removed", expect = "transitions", mutate = function(d) {
    d$transitions[[tm_t(d, "t.structure.05")]] <- NULL
    d
  }),
  list(id = "S12", what = "guard not last", expect = "transitions", mutate = function(d) {
    d$transitions[[tm_t(d, "t.structure.05")]]$order <- 5L
    d$transitions[[tm_t(d, "t.structure.06")]]$order <- 6L
    d
  }),
  list(id = "S13", what = "gap in transition orders", expect = "transitions", mutate = function(d) {
    d$transitions[[tm_t(d, "t.structure.04")]]$order <- 99L
    d
  }),
  list(id = "S14", what = "duplicate question id", expect = "unique_ids", mutate = function(d) {
    d$questions[[length(d$questions) + 1L]] <- d$questions[[tm_q(d, "q.aim")]]
    d
  }),
  list(id = "S15", what = "duplicate transition id", expect = "unique_ids", mutate = function(d) {
    d$transitions[[tm_t(d, "t.structure.03")]]$transition_id <- "t.structure.02"
    d
  }),
  list(id = "S16", what = "duplicate option id", expect = "options", mutate = function(d) {
    qi <- tm_q(d, "q.stage")
    d$questions[[qi]]$options[[length(d$questions[[qi]]$options) + 1L]] <- d$questions[[qi]]$options[[1]]
    d
  }),
  list(id = "S17", what = "planning claims exact_gpower_match", expect = "planning", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.student")]]$planning$planning_compatibility <- "exact_gpower_match"
    d
  }),
  list(id = "S18", what = "planning claims an independently validated extension", expect = "planning",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.two_group.welch")]]$planning$planning_compatibility <-
           "independently_validated_extension"
         d
       }),
  list(id = "S19", what = "classifier outcome with test-power planning", expect = "planning",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.classifier.heldout")]]$planning$planning_availability <-
           "planned_independent_validation"
         d
       }),
  list(id = "S20", what = "outside-scope outcome with planning", expect = "planning", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.survival_outside")]]$planning$planning_compatibility <- "not_yet_validated"
    d
  }),
  list(id = "S21", what = "German question text missing", expect = "bilingual", mutate = function(d) {
    d$questions[[tm_q(d, "q.unit")]]$text_de <- NULL
    d
  }),
  list(id = "S22", what = "German rationale missing", expect = "bilingual", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.ancova")]]$rationale_de <- NULL
    d
  }),
  list(id = "S23", what = "German option label missing", expect = "bilingual", mutate = function(d) {
    qi <- tm_q(d, "q.stage")
    d$questions[[qi]]$options[[tm_opt(d$questions[[qi]], "planning")]]$label_de <- NULL
    d
  }),
  list(id = "S24", what = "recommended outcome without an analysis variant",
       expect = c("outcome_status", "template"), mutate = function(d) {
         d$outcomes[[tm_o(d, "out.ancova")]]$analysis_variant_id <- NULL
         d
       }),
  list(id = "S25", what = "needs_information outcome naming a variant", expect = "outcome_status",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.needs_information")]]$analysis_variant_id <- "t_means_two_independent_welch"
         d
       }),
  list(id = "S26", what = "status outside the vocabulary", expect = "outcome_status", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.welch")]]$status <- "validated"
    d
  }),
  list(id = "S27", what = "code template of a different test", expect = "template", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.welch")]]$r_template_id <- "tpl.t_means_two_independent_student.v1"
    d
  }),
  list(id = "S28", what = "unknown UZH and DescrTab2 source nodes", expect = "source_nodes",
       mutate = function(d) {
         oi <- tm_o(d, "out.two_group.welch")
         d$outcomes[[oi]]$source_nodes <- list(uzh = list("uzh.l.does_not_exist"),
                                               descrtab2 = list("dt2_leaf_does_not_exist"))
         d
       }),
  list(id = "S29", what = "unknown literature reference", expect = "references", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.welch")]]$references <- list("nobody2099")
    d
  }),
  list(id = "S30", what = "unknown crosswalk id", expect = "crosswalk_refs", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.welch")]]$crosswalk_ids <- list("cw999")
    d
  }),
  list(id = "S31", what = "unknown analysis variant as alternative", expect = "variant_refs",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.two_group.welch")]]$alternatives[[1]]$analysis_variant_id <- "t_test_magic"
         d
       }),
  list(id = "S32", what = "override_allowed that is not logical", expect = "variant_refs",
       mutate = function(d) {
         d$outcomes[[tm_o(d, "out.two_group.welch")]]$alternatives[[1]]$override_allowed <- "yes"
         d
       }),
  list(id = "S33", what = "routing on a diagnostic result (input_kind)", expect = "input_kind",
       mutate = function(d) {
         d$questions[[tm_q(d, "q.variance")]]$input_kind <- "diagnostic_result"
         d
       }),
  list(id = "S34", what = "forbidden routing inputs no longer declared", expect = "input_kind",
       mutate = function(d) {
         d$engine$forbidden_routing_inputs <- list()
         d
       }),
  list(id = "S35", what = "route_with_assumption without an assumption text", expect = "dont_know_option",
       mutate = function(d) {
         qi <- tm_q(d, "q.covariates")
         oi <- tm_opt(d$questions[[qi]], "dont_know")
         d$questions[[qi]]$options[[oi]]$assumption_en <- NULL
         d
       }),
  list(id = "S36", what = "dont_know routed as itself", expect = "dont_know_option", mutate = function(d) {
    qi <- tm_q(d, "q.variance")
    d$questions[[qi]]$options[[tm_opt(d$questions[[qi]], "dont_know")]]$routes_as <- "dont_know"
    d
  }),
  list(id = "S37", what = "dont_know policy that invents a default", expect = "dont_know_option",
       mutate = function(d) {
         qi <- tm_q(d, "q.variance")
         d$questions[[qi]]$options[[tm_opt(d$questions[[qi]], "dont_know")]]$dont_know_policy <- "assume_equal"
         d
       }),
  list(id = "S38", what = "dont_know without a German explanation",
       expect = c("dont_know_option", "bilingual"), mutate = function(d) {
         qi <- tm_q(d, "q.aim")
         d$questions[[qi]]$options[[tm_opt(d$questions[[qi]], "dont_know")]]$explanation_de <- NULL
         d
       }),
  list(id = "S39", what = "paired McNemar no longer rejects Boschloo",
       expect = c("rejected_routes", "variant_usage"), mutate = function(d) {
         oi <- tm_o(d, "out.paired.mcnemar")
         d$outcomes[[oi]]$rejected_routes <- Filter(
           function(r) r$analysis_variant_id != "boschloo_two_sided_central", d$outcomes[[oi]]$rejected_routes)
         d
       }),
  list(id = "S40", what = "status modifier that upgrades", expect = "modifiers", mutate = function(d) {
    d$outcomes[[tm_o(d, "out.two_group.welch")]]$status_modifiers[[1]]$effect <- "recommended"
    d
  }),
  list(id = "S41", what = "documented outcome count does not match", expect = "summary_counts",
       mutate = function(d) {
         d$summary$n_outcomes <- d$summary$n_outcomes + 1L
         d
       }),
  list(id = "S42", what = "unknown rules schema version", expect = "schema", mutate = function(d) {
    d$schema_version <- "inferguidR.advisor-rules/9.0.0"
    d
  }),
  list(id = "S43", what = "start question that does not exist", expect = "schema", mutate = function(d) {
    d$engine$start_question <- "q.start"
    d
  }),
  list(id = "S44", what = "question entry that is not an object", expect = "schema", mutate = function(d) {
    d$questions[[length(d$questions) + 1L]] <- "q.extra"
    d
  })
)

test_that("statically detectable defects are rejected", {
  for (case in tm_static_cases) tm_expect_rejected(case)
})

test_that("an expression string is rejected without being evaluated", {
  marker <- tempfile("advisor-eval-marker-")
  d <- tm_doc()
  d$transitions[[tm_t(d, "t.aim.01")]]$when <- sprintf("file.create('%s')", marker)
  d$questions[[tm_q(d, "q.stage")]]$applies_when <- sprintf("isTRUE(file.create('%s'))", marker)
  p <- tm_write_json(d)
  on.exit(unlink(c(p, marker)), add = TRUE)
  v <- tm_expect_rejected(list(id = "E01", what = "expression strings as conditions", expect = "condition"),
                          path = p)
  expect_true(any(grepl("never evaluated", v$problems$message)))
  e <- tryCatch(advisor_rules(rules = p), inferguidR_error = function(e) e)
  expect_s3_class(e, "inferguidR_invalid_input")
  expect_false(file.exists(marker))
})

test_that("defective companion files and unreadable inputs are rejected without errors", {
  # analysis variants
  d <- tm_doc("analysis-variants.json")
  i <- tm_pos(d$variants, "analysis_variant_id", "t_means_two_independent_welch")
  d1 <- d
  d1$variants[[i]]$learning_method_ids <- list("no_such_learning_method")
  tm_expect_rejected(list(id = "F01", what = "learning method id not in the resource registry",
                          expect = "learning_refs"), arg = "variants", path = tm_write_json(d1))
  d2 <- d
  d2$variants[[i]]$used_by <- d2$variants[[i]]$used_by[-1]
  tm_expect_rejected(list(id = "F02", what = "used_by out of step with rules.json", expect = "variant_usage"),
                     arg = "variants", path = tm_write_json(d2))
  d3 <- d
  d3$schema_version <- "inferguidR.analysis-variants/9.0.0"
  tm_expect_rejected(list(id = "F03", what = "unknown variants schema version", expect = "schema"),
                     arg = "variants", path = tm_write_json(d3))
  # source crosswalk
  cw <- tm_crosswalk()
  c1 <- cw
  c1$difference_reason[c1$crosswalk_id == "cw002"] <- ""
  tm_expect_rejected(list(id = "F04", what = "crosswalk row without a difference reason",
                          expect = "crosswalk_reason"), arg = "crosswalk", path = tm_write_csv(c1))
  tm_expect_rejected(list(id = "F05", what = "crosswalk row removed", expect = c("crosswalk_refs", "summary_counts")),
                     arg = "crosswalk", path = tm_write_csv(cw[cw$crosswalk_id != "cw002", ]))
  # unreadable rules files: always valid = FALSE with a schema problem, never an R error
  raw <- paste(readLines(system.file("advisor", "rules.json", package = "inferguidR"), encoding = "UTF-8",
                         warn = FALSE), collapse = "\n")
  unreadable <- list(
    F06 = tm_write_text(sub("{", "{\"rules_id\": \"shadow\",", raw, fixed = TRUE)), # duplicate key
    F07 = tm_write_text(substr(raw, 1, 5000)),                                      # truncated
    F08 = tm_write_text("[]"),
    F09 = tm_write_text("5"),
    F10 = tm_write_text("{\"schema_version\": \"inferguidR.advisor-rules/0.1.0\", \"questions\": [1, 2]}"),
    F11 = "https://example.org/rules.json",
    F12 = file.path(tempdir(), "no-such-advisor-rules.json")
  )
  for (id in names(unreadable)) {
    v <- advisor_validate_rules(rules = unreadable[[id]])
    expect_false(v$valid, info = id)
    expect_true("schema" %in% v$problems$check, info = id)
    e <- tryCatch(advisor_rules(rules = unreadable[[id]]), error = function(e) e)
    expect_s3_class(e, "inferguidR_invalid_input")
  }
  unlink(unlist(unreadable[!names(unreadable) %in% c("F11", "F12")]))
})

